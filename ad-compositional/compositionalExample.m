%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Set up and run a fully implicit compositonal example.
% This ad-compositional module is developed base on the frame of MRST 2017a
% and should be run with the corresponding version. This work can be
% divided into four main parts: (1) compsDataBase contains data of common components
% used in petroleum engineering; (2) compositionalModel includes necessary information and functions;
% (3) equationsCompositional carries the conservation equations, fugacity
% equations and well equations and is responsible for constructing jacobian matrix;
% (4) flashPakage is equipped with flasher and fugacity which
 % can calculate the phase equilibrium and fugacity.
%
% SHORTCOMINGS
% (1) This module cannot handle the flow problem that across the phase
%      boundary.
% (2) Only a simple well with bhp constraint is added
%}
clear;
clc;
close all;
% We begin by loading the required modules
mrstModule clear
mrstModule add ad-compositional ad-core ad-props mrst-gui  
[ state, model, schedule] = compsExampleSetUp;
% Solve
nls = NonLinearSolver('useLinesearch', true);
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', nls);
% Plot after simulation
time = zeros(numel(states),1);
qs = zeros(numel(states),1);
sumqs = zeros(numel(states),1);
for i = 1:numel(states)
    figure(1)
    p = plotCellData(model.G,convertTo(states{i}.pressure, barsa));
     plotWell(model.G, schedule.control(1).W);  
    title('Pressure, bar');
    xlabel('Pressure, bar');
    axis tight, view(35, 40), colorbar('SouthOutside');
   
    figure(2) % Plot for gas rate and cumulative gas versus time    
    time(i) = convertTo(report.ReservoirTime(i), day);
    qs(i) = convertTo(-wellSols{i}.qGs , meter^3/day);
    % Line specification
    lineSpec = {'LineWidth', 2}; 
    markerSpec = {'Marker','o',  'MarkerFaceColor','green', 'MarkerEdgeColor', 'blue'};
    plot(time(1:i), qs(1:i), lineSpec{:}, markerSpec{:});
    % Axis
    xmax = convertTo(report.ReservoirTime(end), day);
    ymax = 1.1*qs(1); 
    axis([0 xmax 0 ymax]);
    title('Surface gas rate');
    xlabel('Time, days');
    ylabel('gas rate, m3/d');
    hold on
   pause(0.3);
end
