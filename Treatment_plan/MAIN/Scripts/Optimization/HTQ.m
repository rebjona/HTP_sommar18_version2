function q = HTQ(P, tumor_tissue,healthy_tissue)
% Calculates HTQ for a PLD field in either mat- or oct-form. 
% ------INPUTS----------------------------------------------
% P:              mat or oct PLD field
% tumor_tissue:   mat or oct with 1 in tumour and 0 otherwise
% healthy_tissue: mat or oct with 1 in healthy tissue and 0 otherwise
% ------OUTPUTS---------------------------------------------
%q:               HTQ value of P
% ----------------------------------------------------------

% Oct/mat conversion
if isa(P, 'Yggdrasil.Octree')
    P = P.to_mat();
end
if isa(tumor_tissue, 'Yggdrasil.Octree')
    tumor_tissue = tumor_tissue.to_mat();
end
if isa(healthy_tissue, 'Yggdrasil.Octree')
    healthy_tissue = healthy_tissue.to_mat();
end
%ÄNDRAT - TILLAGT
load 'tissue_mat_child';
tissue_mat=data;
water_ind = 30;
    ext_air_ind = 1;
    int_air_ind = 5;
    tumor_ind = 9;
healthy_tissue = tissue_mat~=water_ind & ...
    tissue_mat~=ext_air_ind & ...
    tissue_mat~=tumor_ind & ...
    tissue_mat~=int_air_ind;% & ...
healthy_tissue;
%FÄRDIGÄNDRAT - TILLAGT
tumor_tissue = logical(tumor_tissue);
healthy_tissue = logical(healthy_tissue);

% Denomenator in HTQ: Mean in tumor
PLDtumor = P(tumor_tissue);
meanPLDtumor = mean(PLDtumor);

% Nomenator in HTQ: highest PLD percentage in healthy tissue
PLDhealthy = P(healthy_tissue);
PLDv1 = mean(PLDhealthy(PLDhealthy > prctile(PLDhealthy,99)));

q=double(PLDv1/meanPLDtumor);

end