function [modelType,nbrEfields,PwrLimit,goal_function,iteration,hstreshold, SavePath1,particle_settings,freq] = InputData()
% Opens window to enter input data

prompt = {'Model type:',...
    'Number of E-fields:',...
    'Antenna power limit (% of 150 W):',...
    'Goal function (M1-M1, M1-HTQ, M2):',...
    'Iteration number',...
    'Hotspot treshold',...
    'Path to folder for data storage', ...
    'Particle swarm size:', ...
    'Max iterations:', ...
    'Max stall iterations:',...
    'Frequency(ies), MHz, one per row:'};
title = 'Inputs';
num_lines = [1,1,1,1,1,1,1,1,1,1,5];
defaultans = {'child','14','100','M1-M1','1','2.0','Path','20','10','10',['450']};
options.Resize = 'on';
[input] = inputdlg(prompt,title,num_lines,defaultans,options);

modelType = input{1};
nbrEfields = str2num(input{2});
PwrLimit = str2num(input{3})/100;
goal_function = input{4};
iteration=str2num(input{5});
hstreshold=str2num(input{6});
SavePath1=input{7};
particle_settings = [str2num(input{8}),str2num(input{9}),str2num(input{10})];
frequencies = input{11};
f = size(frequencies);
for j = 1:f(1)
    freq(j) = str2num(frequencies(j,:));
end
end
