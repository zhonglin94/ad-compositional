function [problem, state] = equationsCompositional(state0, state, model, dt, drivingForces, varargin)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% DESCRIPTION
% Generate linearized problem for the compositional equations which 
% include the conservation equations, fugacity equations and well equations      
 %}

opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});
% Shorter names for some commonly used parts of the model and forces.
o = model.operators;
f = model.fluid;
W = drivingForces.W;
nc = numel(model.compsVarNames);
ns = numel(model.saturationVarNames);
ng =model.G.cells.num;
gvar = model.compsVarNames;
xcVars = cellfun(@(gvar) strcat('x',gvar), gvar,'UniformOutput' , false);
ycVars = cellfun(@(gvar) strcat('y',gvar), gvar,'UniformOutput' , false);
% Properties at current timestep
[p, sO, sG, xc, yc] = model.getProps(state, ...
    'pressure', 'oil', 'gas', 'xc', 'yc');
% Properties at previous timestep
[p0, sO0, sG0, xc0, yc0, rho0, mw0] = model.getProps(state0, ...
    'pressure', 'oil', 'gas', 'xc', 'yc', 'rho', 'mw');
xc0 = num2cell(xc0,1);
yc0 = num2cell(yc0,1);
s0 = {1-sO0-sG0, sO0,sG0};

xc = num2cell(xc,1);
yc = num2cell(yc,1);
% Initialize primary variabes, include fluid compositions
 [p, sO, sG, xc{1:end-1}, yc{1:end-1}] = initVariablesADI(p, sO, sG, xc{1:end-1}, yc{1:end-1}); 
tempx =0; tempy =0;
for i = 1:nc-1
    tempx =tempx+xc{:,i};
    tempy =tempy+yc{:,i};
end
    xc{:,end} = 1-tempx;
    yc{:,end} = 1-tempy;
    s = {1-sO-sG, sO, sG};
 % Get the gas and oil properties and derivatives for ADI variables
 [fug, ~, rho, mw, mu] = model.initADIProps(p, xc, yc); %#
 rhom = cell(ns,1);
 for  i = 1:ns
     rhom{i} = rho{i}./mw{i};
 end
rhom0 = rho0./mw0;
rhom0 = num2cell(rhom0,1);
% Evaluate relative permeability
[krW, krO, krG] = model.evaluateRelPerm(s);
%
T = o.T;
%Water
mobW   = krW./mu{1};
pW =p;
dpW    = o.Grad(pW);
% water upstream-index
upcw  = (double(dpW)<=0);
vW = -o.faceUpstr(upcw, mobW).*T.*dpW;
% Oil
dpO    = o.Grad(p);
mobO   = krO./mu{2};
% Oil upstream-index
upco = (double(dpO)<=0);
vO = cell(1, nc);
for i = 1:nc
    mobOcomps = mobO.*xc{i};
    vO{i}    = - o.faceUpstr(upco, mobOcomps).*T.*dpO;
end   
% Gas
pG = p;
mobG   = krG./mu{3};
dpG    = o.Grad(pG);
    % Gas upstream-index
upcg    = (double(dpG)<=0);
vG = cell(1,nc);
for i =1:nc
   mobGcomps = mobG.*yc{i};
   vG{i} = - o.faceUpstr(upcg, mobGcomps).*T.*dpG;
end
mob = {mobW, mobO, mobG};
%% Equations 
% The index 1 stands for water
water = (o.pv/dt).*( rhom{1}.*s{1} - rhom0{1}.*s0{1}) + o.Div(o.faceUpstr(upcw, rhom{1}).*vW);
[hydrocarbon, fugacity] = deal(cell(nc,1));
for i =1:nc
% Molar balance equations
% The index 2 and 3 stand for oil and gas respectively
hydrocarbon{i} = (o.pv/dt).*( rhom{2}.*s{2}.*xc{i} - rhom0{2}.*s0{2}.*xc0{i}) + o.Div(o.faceUpstr(upcw, rhom{2}).*vO{i})+...
                          (o.pv/dt).*( rhom{3}.*s{3}.*yc{i} - rhom0{3}.*s0{3}.*yc0{i} ) + o.Div(o.faceUpstr(upcw, rhom{3}).*vG{i});
% Fugacity equations
fugacity{i} = fug{1}{i}-fug{2}{i};
end
%% Insert a simple well
wc = W(1).cells; % Perforation grid cells
WI = W(1).WI;    % Well indices
bhp =W.val;
    qwm = WI .* (mob{1}(wc)).*rhom{1}(wc).* (bhp - p(wc)); % Water molar rate
    water(wc) = water(wc)-qwm;
    [qom, qgm] = deal(cell(1,nc)); % Oil and gas molar rate
for i =1:nc
    qom{i} = WI .* mob{2}(wc).*rhom{2}(wc).*xc{i}(wc) .* (bhp - p(wc));
    qgm{i} =WI .* (mob{3}(wc)).*rhom{3}(wc).*yc{i}(wc) .* (bhp - p(wc));
    hydrocarbon{i}(wc) = hydrocarbon{i}(wc) -qom{i}-qgm{i};
end
%% Flash at surface condition
[qws, qos, qgs] = flashWell(model, qwm, qom, qgm);
state.wellSol.qWs = qws;
state.wellSol.qOs =qos;
state.wellSol.qGs =qgs;
state.wellSol.qs = [qws, qos, qgs];
%% Put the set of equations into cell arrays along with their names/types.
eqs = {water, hydrocarbon{:}, fugacity{:}};
names = {'pressure', 'so','sg', xcVars{1:end-1}, ycVars{1:end-1}};
primaryVars = names;
types = cell(1, 2.*numel(gvar)+1);
% 2nc+1 variables
for i =1:(2*numel(gvar)+1) 
types{i} = 'cell'; 
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

