function create_initial_perf_nonlin(tissue_mat, perfusion_mat, modelType)

% Function that generates a matrix of the perfusion used for calculating
% the very first temperature in CompletePennes. Here, the intial
% temperature is assumed to be 37 deg C.

% model for modifying the perfusion for muscle, fat and tumor is obtained
% from "Impact of Nonlinear Heat Transfer on Temperature Control in
% Regional Hyperthermia" by Lang, Erdmann, Seebass

T=37; %assume initial temp

initial_perf_mat=perfusion_mat;

perf_muscle= 0.45+3.55*exp(-(T-45)^2/12);
perf_fat= 0.36+0.36*exp(-(T-45)^2/12);
perf_tumor= 0.833-(T-37)^4.8/(5.438*10^3);

if startsWith(modelType, 'duke')
    index_muscle = find(tissue_mat==48);
    index_fat = find(tissue_mat==27);
    index_tumor = find(tissue_mat==80);
    
    initial_perf_mat(index_muscle)=perf_muscle;
    initial_perf_mat(index_fat)=perf_fat;
    initial_perf_mat(index_tumor)=perf_tumor;
    
elseif startsWith(modelType, 'child')
    index_muscle = find(tissue_mat==3);
    index_fat = find(tissue_mat==8); % Obs this is WM, cant find fat 
    index_tumor = find(tissue_mat==9);
    
    initial_perf_mat(index_muscle)=perf_muscle;
    initial_perf_mat(index_fat)=perf_fat;
    initial_perf_mat(index_tumor)=perf_tumor;  
end

save(get_path('initial_perf_mat'), 'initial_perf_mat', '-v7.3');

end
