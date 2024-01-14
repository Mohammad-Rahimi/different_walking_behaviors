function [URcond, URmus, UR, syncon]=funur(data,W,C)
%[URcond, URmus, UR, syncon]=funur(data,w,c) 
% This function calculates uncentered correlation coefficients of data and
% reconstructed data = WC
% W and C are used to generate the reconstructed data (recdat= W*C)
% It determines the mean error in the overall reconstruction (UR)
% It determines the error in the reconstruction of each muscle tuning 
% curve (URmus) and each muscle activation pattern for every single
% perturbation direction (URcond)
% Input:
%       data    matrix of observed data  (e.g., data=[mus pert_dir])
%       W       matrix of synergy vectors 
%       C       matrix of coefficiens 
% Output:
%       URcond   matrix with error % for each condition(e.g., error= [pert_dir error])
%       URmus    matrix with error % for each muscle (e.g., error= [mus error])
%       UR       matrix with overall error
% called functions:
%       rsqr_uncentered.m 
%       
% 
% this function is called by:
%       plot_syn.m
%       plot_syn3D.m
%
% Written by: GTO May 24th, 2006
% Last modified:
%
%

[nmuscles ncond]=size(data);
[nsyn ndum]=size(C);

%Calculate reconstructed values
ReconData=W*C;

Wstr=mat2cell(W,[nmuscles],[ones(1,nsyn)]); % <1x7>
Cstr=mat2cell(C,[ones(1,nsyn)],[ncond]);
Recon_parts=cellfun(@(N1,N2) N1*N2,Wstr,Cstr','UniformOutput',0);
% coef_det_uncentered = @(recon,orig) 1-((sum((recon(:).*orig(:)).^2))/(sum(orig(:).^2)*sum(recon(:).^2)));
% synconstr=cellfun(@(D1) coef_det_uncentered(D1,data),Recon_parts,'UniformOutput',0);
% syncon=cell2mat(synconstr);

% Calculate contribution of each synergy to overall variability
synconstr=cellfun(@(D1) 100*sum(sum(D1))/sum(sum(ReconData)),Recon_parts,'UniformOutput',0);
% synconstr=cellfun(@(D1) 100*sum(sum(D1))/sum(sum(data)),Recon_parts,'UniformOutput',0);
syncon=cell2mat(synconstr);


% synconstr=cellfun(@(D1) overallr2uncentered(data,D1),Recon_parts,'UniformOutput',0);
% syncon=cell2mat(synconstr);


%Calculate overall variability(1x1)
UR=overallr2uncentered(data,ReconData);

%Make fake reconstructed data with 70% error in the prediction
%ReconData=data.*1.7;
    
%Calculate error in the reconstruction of each direction
%URcond(1 x nconditions)
[URcond]=rsqr_uncentered(data',ReconData');
URcond=100*(URcond);   

%Calculate error in the reconstruction of each muscle activity level
%URmus(nmus x 1)
[URmus]=rsqr_uncentered(data,ReconData);    
URmus=100*(URmus);





