% MAIN
% This script contains all necessary steps to complete a hyperthermia
% treatment plan. 

% Run all or run in sections
% Current directory should be MAIN !
% ------------------------------------------------------------------------

% Inputs and compilation
clear all
addpath ..
olddir = pwd;
filename = which('Main');
[mainpath,~,~] = fileparts(filename);

cd(mainpath)
addpath([mainpath filesep 'Scripts' filesep 'Optimization'])

[modelType,nbrEfields,PwrLimit,goal_function,iteration,hstreshold,particle_settings,freq] = InputData;
hyp_compile
hyp_init

%% Optimization
%Ändrat
fclose('all');
save ('iteration','iteration');
cd Results\P_and_unscaled_settings
delete settings_child_450MHz.txt P_child_450MHz.mat E_child_450MHz.oct
cd(mainpath);
%Färdigändrat
if length(freq) ==1 
    EF_optimization_single(freq, nbrEfields, modelType, goal_function, particle_settings)
elseif length(freq) ==2
    EF_optimization_double(freq, nbrEfields, modelType, goal_function, particle_settings)
elseif length(freq) >2
    error('Optimization does not currently work for more than two frequencies. Prehaps combine_single can be of use?')
end

%% Generate FEniCS Parameters
fclose('all');
cd Results\Input_to_FEniCS\child;
delete *.mat *.txt *.xml *.obj;
cd(mainpath)
message1 = msgbox('Optimization finished! Generating FEniCS parameters. ','Success');
t = timer('ExecutionMode', 'singleShot', 'StartDelay',3,'TimerFcn',@(~,~)close(message1));
start(t);
isopath = [mainpath filesep '..' filesep 'Libs' filesep 'iso2mesh'];
fenicspath = [mainpath filesep 'Scripts' filesep 'Prep_FEniCS']; 
preppath = [fenicspath filesep 'Mesh_generation'];
addpath(isopath, fenicspath, preppath) 

generate_fenics_parameters(modelType, freq, true)
generate_amp_files(modelType, freq, nbrEfields, PwrLimit)
cd(olddir)
%% FEniCS
message2 = msgbox(['Time to move on to FEniCS! Remember to input the files you just generated. '...
    'Then press any key to continue (In command window).'],'Time for FEniCS');
pause 

%% Temperature
temppath = [mainpath filesep 'Scripts' filesep 'Temperature'];
addpath(temppath);
evaluate_temp(modelType, freq, true);
save_scaled_settings(modelType, freq, nbrEfields);
disp('Finished, storing temperature files and settings...')
HotSpotExtractor_child(hstreshold, iteration, freq);

   
message4 = msgbox('You have now finished the optimization and temperature transformation. You rock!', 'Hot stuff!');
u = timer('ExecutionMode', 'singleShot', 'StartDelay',5,'TimerFcn',@(~,~)close(message4));
start(u);
LocateLeo = [mainpath filesep 'Scripts'];
GoodJob = [LocateLeo filesep 'Leonardo-DiCaprio-Clap.gif'];
fig = figure;
gifplayer(GoodJob,0.1);
r = timer('ExecutionMode', 'singleShot', 'StartDelay',1.5,'TimerFcn',@(~,~)close(fig));
start(r);0
cd(olddir)
fclose('all');
