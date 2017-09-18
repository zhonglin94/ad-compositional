function [ state, model, schedule ] = compsExampleSetUp
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% We define a 6x6x3 grid, spanning a 120*120*12 domain with homogenous
% permeability and porosity. The effect of gravity is ignored. For
% illustrative purpose, we added a simple well with a bhp control of 70
% bars and no well equations are needed. The reservoir
% fluids is always in two phase region during the simulation. 
%
% SYNOPSIS:
%    [state, model, schedule] = compsExampleSetUp
%
% RETURNS:
%     state - The state of reservoir at a time step
%   model - Compositional model
% schedule - Include well and control information
%}
    cartDims = [6,6,3];
    physDims = [120, 120, 12]*meter;
    G = cartGrid(cartDims, physDims);
    G = computeGeometry(G);
    gravity reset off
    
    poro = 0.2;
    perm = 1.2*milli*darcy;
    rock = makeRock(G, perm, poro);

    %% Set up wells and schedule
    T = 1000*day;
     % Producer
     W = [];
     W = verticalWell(W, G, rock, 1, 1, 3, 'Name', 'P1', 'radius', 5*inch, ...
        'Type', 'bhp', 'Val', 70*barsa, 'Sign', -1);   % Producer

    dT_target = 50*day; 
    dt = rampupTimesteps(T, dT_target, 10);
    
    % Set up a schedule
     schedule = struct();     
     schedule.step.val = dt; 
     schedule.step.control = (cumsum(dt) > T) + 1;
     schedule.control = struct('W', W);

     % Input properties of each component
     %          C1      n-C4   n-C10
     zi =   [0.5532  0.3630 0.0838]; %Mole Fraction
     % Id represent the index of components in 'Components_db'
     id = [5 9 16];
     [mix, names] = getCompsProps(zi, id);
     %Relative permeability
     props = getRelPermTable();     
     props.DENSITY = [0 1.0378e+03  0];
     props.PVTW = [2.768e+07,1.029,4.539e-10,3.100e-04,0];
     props.ROCK = [7.865e+02,1.038e+03,0.9698];
      fluid = initCompsADIFluid(props); 
      EoS = 'PR';
     model = compositionalModel(G, rock, fluid, 'EOS', EoS, 'compsVarNames', names, 'Verbose', true); 
     % Temperature
     t =60+273.15;  %Kelvin
     model.inputdata.t = t;
     % Water and hydrocarbon properties
     props.mix= mix;
     model.inputdata.props = props;
     model.inputdata.mix = mix;
     % Initial the pressure and water saturation
     state = initResSol(G, 120*barsa, [0.2, 0, 0]);
     % Overall hydrocarbon composition
     state.zc = repmat(zi, model.G.cells.num,1); 
     % Implement flash for state and then get water, oil and gas properties
     state = flashState(model, state);
end