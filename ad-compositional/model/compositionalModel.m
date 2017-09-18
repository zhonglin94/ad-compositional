classdef compositionalModel < ReservoirModel
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% DESCRIPTION
% Full implicit compositional model
%}
 properties
     % Equation of state
     Eos    
     % Compositional variable names
     compsVarNames
    % Maximum relative compositional variable increment
    dcMaxRel
     % Maximum relative compositional variable increment
    dcMaxAbs   
end   
methods
    function model = compositionalModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid, varargin{:});
        % Equation of state
        model.Eos = [];
        model.compsVarNames =[];

        % Max increments of components
        model.dcMaxAbs = inf;
        model.dcMaxRel = inf;
        % All phases are present 
        model.oil = true;
        model.gas = true;
        model.water = true;
        model.saturationVarNames = {'sw', 'so', 'sg'};       
        
        % Compositional -> use CNV style convergence 
        model.useCNVConvergence = true;         
        
        model = merge_options(model, varargin{:});
    end
    
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsCompositional(state0, state, model, dt, ...
                        drivingForces, varargin{:}); 

    end
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        compsVars = lower(model.compsVarNames);
        satVars = model.saturationVarNames;
        switch (lower(name))
            case compsVars
                    fn = 'zc';
                    index = model.getCompsVarsIndex(name);      
            case {'zc','xc','yc', 'ki'}
                fn = lower(name);
                index = 1:numel(compsVars);
            case {'alpha','rhoo','rhog', 'mwo','mwg'}
                fn = lower(name);
                index = 1;         
            case{'rho', 'mw', 'mu'}
                fn = lower(name);
                index = 1:numel(satVars);
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name);
        end
   end
   % --------------------------------------------------------------------%        
    function compsIndex = getCompsVarsIndex(model, name)
        compsVars = lower(model.compsVarNames);
        compsIndex = find(strcmpi(compsVars, name), true);
    end
    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@ReservoirModel(model, state);
        ng = model.G.cells.num; % Number of grid
        nc = numel(model.compsVarNames); % Numer of components
            % zc must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'zc', [ng, nc], [1, 2]);
    end   
    % --------------------------------------------------------------------%
    function [sO, sG] = getOGSat(model, st, alpha, Z, p, t, sW)
        % Get oil and gas saturation
        R  = 8.314472; %Pa.m3/mol.K
        vOmolar = (R.*t.*Z(:,1))./p; % m3/mol
        vGmolar = (R.*t.*Z(:,2))./p;
        sG = (st == 2).*(1-sW)+(st == 3).*alpha.*vGmolar.*(1-sW)./(alpha.*vGmolar+(1-alpha).*vOmolar); %OK
        sO = (st ==1).*(1-sW)+(st ==3).*(1-sG-sW);      
    end
    % --------------------------------------------------------------------%
    function [muO, muG] = getOGViscosity(model,MW,rhoO,rhoG,rhocO,rhocG, xc,yc, t)         
           Tc = model.inputdata.mix.Tc;
           Pc = model.inputdata.mix.Pc;
          % Reduced density
         rhorO = rhoO./rhocO;
         rhorG = rhoG./rhocG;
         muO = LohrenzBrayClark(rhorO, xc, t, Tc, Pc, MW); 
         muG = LohrenzBrayClark(rhorG, yc, t, Tc, Pc, MW);
    end      
   % --------------------------------------------------------------------%
    function [rhoO, rhoG, rhocO, rhocG, mwO, mwG] = getOGDensity(model, MW, xc, yc, Z, p, t)
        R  = 8.314472; %Pa.m3/mol.K
        MW = MW.*milli;% kg/mol
        Vc = model.inputdata.mix.Vc.*milli; % l/mole-> m3/mol
        [mwO, mwG, VcO, VcG] = deal(0);
        nc  = numel(model.compsVarNames);
        for  i = 1: nc
            %Sum of molar weight
            mwO = xc{i}.*MW(i)+mwO;  
            mwG = yc{i}.*MW(i)+mwG;
            % Sum of critical volume
            VcO = xc{i}.*Vc(i)+VcO;       
            VcG = yc{i}.*Vc(i)+VcG;             
        end       
        if isnumeric(Z)
            % oil density             
            rhoO = mwO./((R.*t.*Z(:,1))./p); rhocO = mwO./VcO;
            % gas density
            rhoG = mwG./((R.*t.*Z(:,2))./p); rhocG = mwG./VcG;
        else
            rhoO = mwO./((R.*t.*Z{1})./p); rhocO = mwO./VcO;
            rhoG = mwG./((R.*t.*Z{2})./p); rhocG = mwG./VcG;    
        end
    end    
    % --------------------------------------------------------------------%              
       function [fug, Z, rho, mw, mu] = initADIProps(model, p, xc, yc)
        t =model.inputdata.t;    
        MW = model.inputdata.mix.MW;   
        [fugL, ZL] = model.fugacity(xc, p, t, 1); 
        [fugV, ZV] = model.fugacity(yc, p, t, 0); 
        % Oil and gas fugacity [ADI]
        fug = {fugL, fugV};
        Z = {ZL, ZV};
        % Water, oil and gas density [ADI]
        [rhoO, rhoG, rhocO, rhocG, mwO, mwG] = model.getOGDensity(MW, xc, yc, Z, p, t);          
         f = model.fluid;
         pW =p; % Ignore capillary pressure
         rhoW = f.bW(pW).*f.rhoWS;
         rho ={rhoW,rhoO,rhoG};       
          % Water, oil and gas weight [ADI] 
         mwW = 18.*milli.*ones(model.G.cells.num,1);
         mw = {mwW, mwO, mwG}; 
         % Water, oil and gas viscosity [ADI] 
         [muO, muG] = model.getOGViscosity(MW,rhoO,rhoG,rhocO,rhocG, xc,yc, t); 
         muW = f.muW(pW);
         mu = {muW,muO,muG}; 
       end
        % --------------------------------------------------------------------%         
        function   [rho, mw, mu, s] = propsAfterFlash(model, p, t, st, alpha, Z, xc, yc, sW)
            % Get the density and viscosity for oil and gas  
            MW = model.inputdata.mix.MW;   
            if ~iscell(xc) && ~iscell(yc)
            xc = num2cell(xc, 1); 
            yc = num2cell(yc,1);
            end
             % Oil and gas densities, Critical densities, mole weights
             [rhoO, rhoG, rhocO, rhocG, mwO, mwG] = model.getOGDensity(MW, xc, yc, Z, p, t); 
             f = model.fluid;
             pW =p; % Ignore capillary pressure
             rhoW = f.bW(pW).*f.rhoWS;
             rho =[rhoW,rhoO,rhoG];
             mwW = 18.*milli.*ones(numel(mwO),1); 
             % Water, oil and gas weight 
             mw = [mwW, mwO, mwG];
             % Water, oil and gas viscosities
             [muO, muG] = model.getOGViscosity(MW,rhoO,rhoG,rhocO,rhocG, xc,yc, t);   
             muW = f.muW(pW);
             mu = [muW,muO,muG];
             % Oil and gas saturations
             [sO, sG] = model.getOGSat(st, alpha, Z, p, t, sW); % Oil and gas saturation         
             s = [sO, sG];
        end
        % --------------------------------------------------------------------%      
        function [xc, yc] = splitPhase(model, zc, ki, alpha)
            nc  = numel(model.compsVarNames);
            ng = model.G.cells.num;
            [xc,yc] = deal(zeros(size(zc,1), nc));
            for i =1:nc
            xc(:,i) = zc(:,i)./(1+(ki(:,i)-1).*alpha);
            yc(:,i) = ki(:,i).*xc(:,i);            
            end   
        end    
    
        % --------------------------------------------------------------------%    
        function state = flashState(model, state)
           % Implement flash for each block
            zc = num2cell(state.zc,2); 
            p = num2cell(state.pressure);
           [ ki, alpha, st, Z] =cellfun(@(zc, p)...
                                            flasher(model, zc, p), zc, p,  'UniformOutput', false);
            state.ki =cell2mat(ki);
            state.alpha = cell2mat(alpha);
            state.status =cell2mat(st);
            state.Z =cell2mat(Z);
            [xc, yc] = model.splitPhase(state.zc, state.ki, state.alpha);
            state.xc = xc;
            state.yc = yc;                
            % Update properties for state
            t =model.inputdata.t;      
            sW = state.s(:,1);
            [rho, mw, mu, s] = model.propsAfterFlash(state.pressure, t, state.status,...
                                               state.alpha, state.Z, state.xc, state.yc, sW);
            state.rho = rho;
            state.mw = mw;
            state.mu =mu;
             s = [sW, s];
             state.s = s;
        end
        % --------------------------------------------------------------------%
        function [qWs, qOs, qGs] = flashWell(model,qWm, qOm, qGm)
            % Overall hydrocarbon composition    
            qWm = double(qWm);
            % Convert ADI to double
            qOm = cellfun(@(x) double(x), qOm);   
            qGm =cellfun(@(x) double(x), qGm);    
            zc = qOm+qGm;
            zt = sum(zc);
            % Normalize to 1 and convert to cell
            zc = zc./zt;        
            % Water, oil and gas volume rate
            p = 1.*atm; t = (25+273.15).*Kelvin; % Surface condition
            [ ki, alpha, st, Z] = flasher(model, zc, p, t);
            [xc, yc] = model.splitPhase(zc, ki, alpha);
            % Get fluid properties at surface condition
            [rho, mw, mu, ~] = propsAfterFlash(model, p, t, st, alpha, Z, xc, yc, []);
            rhom = rho./mw;
            qWs = qWm./rhom(1);
            qOs= zt.*(1-alpha)./rhom(2);
            qGs = zt.*alpha./rhom(3);    
        end
        % --------------------------------------------------------------------%
        function [restVars, xcVars, ycVars,satVars, wellVars] = splitPrimaryVariables(model, vars)
            % Split a set of primary variables into three groups:
            % Well variables, saturation variables, gas composition variables and 
            % oil compOsition variables and the rest. This is useful because the
            % saturation variables oil and gas composition variables usually are updated
            % together, and the well variables are a special case.
            gvar = model.getCompsVarNames;
            xcVars = cellfun(@(x) strcat('x',x), gvar,'UniformOutput' , false);
            ycVars = cellfun(@(x) strcat('y',x), gvar,'UniformOutput' , false);

            isXc = cellfun(@(x) any(strcmpi(xcVars, x)), vars);
            isYc = cellfun(@(x) any(strcmpi(ycVars, x)), vars);
            xcVars = vars(isXc);
            ycVars = vars(isYc);
            for i = 1:numel(xcVars)
                vars = model.stripVars(vars,xcVars{i});
                vars = model.stripVars(vars,ycVars{i});
            end
            [restVars, satVars, wellVars] = splitPrimaryVariables@ReservoirModel(model, vars);
        end   
        % --------------------------------------------------------------------%
            function vars = getCompsVarNames(model)
                    vars = model.compsVarNames;
            end
        % --------------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
             % Get the density and viscosity for oil and gas
            p =state.pressure;
            ng = model.G.cells.num;
            xc =state.xc;
            yc =state.yc;
            state.ki = xc./yc;
            if ~iscell(xc) && ~iscell(yc)
                xc = num2cell(xc, 1); % 
                yc = num2cell(yc,1);
            end        
            [fug, Z, rho, mw, mu] = model.initADIProps(p, xc, yc);
            state.rho = cell2mat(rho);
            state.mw =cell2mat(mw);
            state.mu =cell2mat(mu);
            state.Z = cell2mat(Z);
            rhom = state.rho./state.mw;
            s =state.s;
            state.alpha = rhom(:,3).*s(:,3)./( rhom(:,2).*s(:,2)+ rhom(:,3).*s(:,3));
            % Update overall composition
            for i = 1:ng
                state.zc(i,:) =  (rhom(i,3).*s(i,3).*state.yc(i,:)+rhom(i,2).*s(i,2).*state.xc(i,:))./( rhom(i,2).*s(i,2)+ rhom(i,3).*s(i,3));
            end
            report =[];
        end
        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
         [restVars, xcVars, ycVars,satVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables);
         % Update mole fraction in one go
        state  = model.updateXc(state, dx, problem, xcVars);   
        state  = model.updateYc(state, dx, problem, ycVars);       
        state = model.updateSat(state, dx, problem, satVars);
        if ~isempty(restVars)
            % Handle pressure seperately
            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);
            restVars = model.stripVars(restVars, 'pressure');

            % Update remaining variables (tracers, temperature etc)
            for i = 1:numel(restVars);
                 state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
            end
        end
        report = [];
        end
        % --------------------------------------------------------------------%        
        function state = updateSat(model, state, dx, problem, satVars)
                ds = zeros(model.G.cells.num, numel(satVars)+1);
                for i =2:numel(satVars)+1
                    ds(:,i) = model.getIncrement(dx, problem, satVars{i-1}); % Water, oil and gas
                end
                ds(:,1) = -sum(ds(:,2:end),2);
                state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMaxRel, model.dsMaxAbs);
        end
    % --------------------------------------------------------------------%        
    function state = updateXc(model, state, dx, problem, xcVars)
        if nargin < 5
            % Get the xc names directly from the problem
            [~, xcVars] = ...
                splitPrimaryVariables(model, problem.primaryVariables);
        end        
        if isempty(xcVars)
            % No saturations passed, nothing to do here.
            return
        end        
        dxc = zeros(model.G.cells.num, numel(xcVars)+1);
        for i = 1:numel(xcVars)
             v = model.getIncrement(dx, problem, xcVars{i});
             dxc(:,i) =v;
        end
         dxc(:,end) = -sum(dxc,2);
         state = model.updateStateFromIncrement(state, dxc, problem, 'xc', model.dcMaxRel, model.dcMaxAbs);
        % We update all mole fractions simultanously    
    end
    % --------------------------------------------------------------------%   
    function state = updateYc(model, state, dx, problem, ycVars)
        if nargin < 5
            % Get the yc names directly from the problem
            [~, ~, ycVars] = ...
                splitPrimaryVariables(model, problem.primaryVariables);
        end        
        if isempty(ycVars)
            % No saturations passed, nothing to do here.
            return
        end        
        dyc = zeros(model.G.cells.num, numel(ycVars)+1);
        for i = 1:numel(ycVars)
             v = model.getIncrement(dx, problem, ycVars{i});
             dyc(:,i) =v;
        end
         dyc(:,end) = -sum(dyc,2);
         state = model.updateStateFromIncrement(state, dyc, problem, 'yc', model.dcMaxRel, model.dcMaxAbs);
        % We update all mole fractions simultanously    
    end
    % --------------------------------------------------------------------%   
    function state = setProp(model, state, name, value)
        [fn, index] = model.getVariableField(name);
        if strcmpi(fn,'zc')
            state.(fn)(:, index) = num2cell(value,1);
        else
            state = setProp@PhysicalModel(model, state, name, value);
        end
    end
    % --------------------------------------------------------------------%         
    function [fug, Z] = fugacity(model, zi, p, t, phase)
        % Return the fugacities and compressibility factor of a given mixture. 
        % Inputted zi should be in cell
        Eos = model.Eos;
        mix = model.inputdata.props.mix;
        fn = str2func(['fugacity', Eos]);
        [fug, Z]  = fn(zi, mix, p, t, phase);
    end
end
end

       