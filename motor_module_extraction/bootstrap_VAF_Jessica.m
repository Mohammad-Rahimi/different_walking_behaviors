function [synergyDataOut] = bootstrap_VAF_Jessica(synergyDataIn, nsyn, nboot,VAFtype)
%
% [synergyDataOut] = bootstrap_VAF_Jessica(synergyDataIn, nsyn, nboot)
%
% This function is based off of code from SAC (bootstrap_VAF.m and
% bootstrap_VAF_spinalcat.m). 
%
% Put what it does here!
%
% INPUTS:   synergyDataIn - previously identified synergy data structure
%           nsyn - number of synergies (max, e.g. 16 for 16 muscles)
%           nboot - number of bootstrap samples
%           VAFtype - when to calculate VAF (1=after removing UV scaling (default)
%                     2 = before removing UV scaling)
%
% OUTPUTS:  synergyDataOut - synergyDataIn modified to have bootstrapped
%               VAF sorted for both unshuffled and shuffled data
%
%
% Created:  3/20/13, JLA
% Modified: 5/19/15, JLA, changed VAF to be calculated on UV scaled data
%           3/3/16, JLA - added option to calculate VAF before/after UV
%                         scaling
    


% nboot = 500;

VAFtot = [];
VAFmus_tot = [];
sortedVAF = [];
sortedVAFmus = [];

if nargin < 4, VAFtype=1; end

global stdev;

for s = 1:nsyn
    
    DATA = synergyDataIn(s).emgOrig;
    [nmus, nconds] = size(DATA);
    
    orig.DATA = synergyDataIn(s).emgOrig;
    stdev = std(synergyDataIn(s).emgOrig,0,2);
    orig.W = synergyDataIn(s).W;
    orig.C = synergyDataIn(s).C;
    orig.stdev = stdev;
    
    % added JA, 5/19/15
    orig.DATA_uv = diag(1./orig.stdev)*orig.DATA;
    orig.W_uv = diag(1./orig.stdev)*orig.W;
    
        
    % generate CIs for original data
    VAFall = [];
    VAFmusall = [];
    
    % generate random index vector
    bootsam = ceil(nconds*rand(nconds,nboot));
    
    parfor iter=1:nboot
        
        % sample the data and use the same sample from C
        Csamp = orig.C(:,bootsam(:,iter));
        if VAFtype == 1
            DATAsamp = orig.DATA(:,bootsam(:,iter));
            [~, VAFmus, VAF] = funur(DATAsamp,orig.W,Csamp);
        elseif VAFtype == 2 % JA, 5/19/15
            DATAsamp = orig.DATA_uv(:,bootsam(:,iter)); 
            [~, VAFmus, VAF] = funur(DATAsamp, orig.W_uv,Csamp); 
        end
                 
        % store parameters
        VAFall(iter) = VAF;
        VAFmusall(:,iter) = VAFmus;
        
        % clear temporary variables
%         clear DATAsamp Csamp VAF VAFmus
        DATAsamp = [];
        Csamp = [];
        VAF = [];
        VAFmus = [];
        
    end
    
    orig.VAFtot(:,s) = VAFall;
    orig.sortedVAF(:,s) = sort(VAFall);
    synergyDataIn(s).sortedVAF_original = sort(VAFall);
    
    orig.VAFmus_total(:,:,s) = VAFmusall;
    orig.sortedVAFmus(:,:,s) = sort(VAFmusall,2,'ascend');
    synergyDataIn(s).sortedVAFmus_original = sort(VAFmusall,2,'ascend');
    
    
    % shuffle the source data    
    DATAshuff = zeros(size(DATA));
    ind = zeros(size(DATA));
    for i=1:nmus
        ind(i,:) = randperm(nconds);
        DATAshuff(i,:) = DATA(i,ind(i,:));
    end

    % added JA, 5/19/15
    % shuffle the source data    
%     DATAshuff = zeros(size(orig.DATA_uv));
%     ind = zeros(size(orig.DATA_uv));
%     for i=1:nmus
%         ind(i,:) = randperm(nconds);
%         DATAshuff(i,:) = orig.DATA_uv(i,ind(i,:));
%     end
%     
    
    % extract muscle synergies from shuffled data
    [Wshuff, Cshuff, errshuff] = nnmf_Jessica(DATAshuff, s, 2);
    
    shuff.DATA = DATAshuff;
    shuff.W = Wshuff;
    shuff.C = Cshuff;
    shuff.stdev = stdev;
    
    % generate CIs for shuffled data
    VAFall = [];
    VAFmusall = [];
    
    % generate random index vector
    bootsam = ceil(nconds*rand(nconds,nboot));
    
    for iter = 1:nboot
%         iter
%         for i=1:nmus
%             ind(i,:) = randperm(nconds);
%             DATAshuff(i,:) = DATA(i,ind(i,:));
%         end        
%         
%         % extract muscle synergies from shuffled data
%         [Wshuff, Cshuff, ~] = nnmf_Jessica(DATAshuff, s, 2);
%         
%         shuff.DATA = DATAshuff;
%         shuff.W = Wshuff;
%         shuff.C = Cshuff;
%         shuff.stdev = stdev;
        
        % sample the data and use the same sample from C
        DATAsamp = DATAshuff(:,bootsam(:,iter));
        Csamp = Cshuff(:,bootsam(:,iter));
        
%         % reconstruct data and calculate VAF and VAF muscle

        if VAFtype == 1
            [~, VAFmus, VAF] = funur(DATAsamp, Wshuff, Csamp);
        elseif VAFtype == 2 % added JA, 5/19/15
            stdev = std(DATAsamp'); 
            [~, VAFmus, VAF]=funur(diag(1./stdev)*DATAsamp,diag(1./stdev)*Wshuff,Csamp); 
        end

        % store params
        VAFall(iter) = VAF;
        VAFmusall(:,iter) = VAFmus;
        
        % clear temporary variables
%         clear DATAsamp Csamp VAF VAFmus
        DATAsamp = []; Csamp = []; VAF = []; VAFmus = [];
        
    end
    
    shuff.VAFtot(:,s) = VAFall;
    shuff.sortedVAF(:,s) = sort(VAFall);
    synergyDataIn(s).sortedVAF_shuffled = sort(VAFall);
    
    shuff.VAFmus_tot(:,:,s) = VAFmusall;
    shuff.sortedVAFmus(:,:,s) = sort(VAFmusall,2,'ascend');
    synergyDataIn(s).sortedVAFmus_shuffled = sort(VAFmusall,2,'ascend');
    
%     clear Wnew C err bootsam
    Wnew = []; C = []; err = []; bootsmp = [];
end


% % 95% CI
% %sampint=[1,10];
% %sampint=[3,98];
sampint=[13,488];

% % 99% CI
% sampint=[6,495];


synergyDataOut = synergyDataIn;


% figure
% errorbar(1:16,[synergyOutRight(:).VAF],[synergyOutRight(:).VAF]-orig.sortedVAF(14,:),orig.sortedVAF(488,:)-[synergyOutRight(:).VAF])
