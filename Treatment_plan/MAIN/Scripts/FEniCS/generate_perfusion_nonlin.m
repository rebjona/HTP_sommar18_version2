function perfusion_mat=generate_perfusion_nonlin(perf_mat, temp_mat, tissue_mat, modelType)
% ha bas W för första T beräkning, ta in nya T i matlab för att beräkna
% nytt W och sen använda det för beräkna tredje T, repeat
% spara W för varje temp-steg så det går läsa in i CompletePennes

% All values and formula for estimating the perfusion are obtained by the
% report "Impact of Nonlinear Heat Transfer on Temperature Control in
% Regional Hyperthermia" by Lang, Erdmann, Seebass

if startsWith(modelType, 'duke')
    index_muscle = 48;
    index_fat = 27;
    index_tumor = 80;  
    
elseif startsWith(modelType, 'child')
    index_muscle = 3;
    index_fat = 8; % Obs this is WM, cant find fat in tissue_file
    index_tumor = 9;
end

perfusion_mat=perf_mat;

% Change perfusion for the tumor -------------------
tumor=(tissue_mat==index_tumor);

% If temp in tumor < 37C
temp_lower_37=(temp_mat<37);
index_tumor_cold= find(temp_lower_37==tumor & temp_lower_37==1);
perfusion_mat(index_tumor_cold)=0.833; 

% If 37C <= temp in tumor <= 42C
temp_leq_42=(temp_mat<=42);
temp_geq_37=(temp_mat>=37);
index_tumor_middle= find(temp_leq_42==tumor & temp_geq_37==temp_leq_42 & temp_leq_42==1);
perfusion_mat(index_tumor_middle)=0.833-(((temp_mat(index_tumor_middle)-37).^4.8)/(5.438*10^3));% assumed 10^3 from report 

% If temp in tumor > 42C
temp_higher_42=(temp_mat>42);
index_tumor_hot= find(temp_higher_42==tumor & temp_higher_42==1);
perfusion_mat(index_tumor_hot)=0.416; 


% Change perfusion for muscle -----------------------
muscle=(tissue_mat==index_muscle);

% If temp in muscle <= 45C 
temp_leq_45=(temp_mat<=45);
index_muscle_less= find(temp_leq_45==muscle & temp_leq_45==1);
perfusion_mat(index_muscle_less)=0.45+3.55*exp(-((temp_mat(index_muscle_less)-45).^2)/12); 


% If temp in muscle > 45C
temp_higher_45=(temp_mat>45);
index_muscle_hot=find(temp_higher_45==muscle & temp_higher_45==1);
perfusion_mat(index_muscle_hot)=4;

% Change perfusion for fat --------------------------
fat=(tissue_mat==index_fat);

% If temp in fat <= 45C
index_fat_less=find(temp_leq_45==fat & temp_leq_45==1);
perfusion_mat(index_fat_less)=0.36+0.36*exp(-((temp_mat(index_fat_less)-45).^2)/12);

% Is temp in fat > 45C
index_fat_hot=find(temp_higher_45==fat & temp_higher_45==1);
perfusion_mat(index_fat_hot)=0.72;

%-----------------------------------------------------
current_folder=pwd;
savepath=  [pwd filesep '..' filesep 'Input_to_FEniCS' filesep 'perfusion_mat_nonlin'];
save(savepath, 'perfusion_mat', '-v7.3');
% borde man spara in m olika index så man kan kolla på dom senare?

end


