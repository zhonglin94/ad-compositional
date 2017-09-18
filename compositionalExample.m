%% Set up and run a fully implicit compositonal example

clear;
clc;
close all;

% We begin by loading the required modules
mrstModule clear
mrstModule add ad-compositional ad-core ad-props mrst-gui  
[ state, model, schedule] = compsExampleSetUp;
% solve
nls = NonLinearSolver('useLinesearch', true);
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', nls);
% plot after simulation
time = zeros(numel(states),1);
qs = zeros(numel(states),1);
xc = zeros(numel(states),numel(state.xc(1,:)));
yc = zeros(numel(states),numel(state.yc(1,:)));
ki = zeros(numel(states),numel(state.ki(1,:)));
Z = zeros(numel(states),numel(state.Z(1,:)));
sumqs = zeros(numel(states),1);
for i = 1:numel(states)
    xc(i,:) = states{i}.xc(1,:);
    yc(i,:) = states{i}.yc(1,:);
    ki(i,:) =states{i}.ki(1,:);
    Z(i,:) = states{i}.Z(1,:);
    figure(1)
    p = plotCellData(model.G,convertTo(states{i}.pressure, barsa));
     plotWell(model.G, schedule.control(1).W);    % Pick the only well control present
    title('Pressure, bar');
    xlabel('Pressure, bar');
    axis tight, view(35, 40), colorbar('SouthOutside');
   
    figure(2) % Plot for gas rate and cumulative gas versus time    
    time(i) = convertTo(report.ReservoirTime(i), day);
    qs(i) = convertTo(-wellSols{i}.qGs , meter^3/day);
    
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
   pause(0.2);
end
