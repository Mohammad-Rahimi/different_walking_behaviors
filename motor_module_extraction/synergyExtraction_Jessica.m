function [synergyOut] = synergyExtraction_Jessica(V, numSyn, VAFtype)
%
% [synergyOut] = synergyExtraction_Jessica(V,numSyn)
%
% INPUTS:   V - emg matrix resulting from calcSynergyEMGMatrix
%           numSyn - total # of synergies (will loop through and extract
%                    synergies for 1 - numSyn)
%           VAFtype - when to calculate VAF (1=after removing UV scaling (default)
%                     2 = before removing UV scaling)
%
% OUTPUTS:  synergyOut - structure containing synergy information for each
%                       number of synergies
%
% Modified and compliled from various sources. Requires function
% nnmf_Jessica to run the actual NNMF algorithm on the data. 
%
% Calls the following functions:
%   nnmf_Jessica - for running the NNMF algorithm
%   funur - to calculate VAFs
%   funR - to calcualted r^2 values
%
% Created: 3/13/13, JLA
% Modified: 3/14/13, JLA - added VAF muscle per bin calculation. added
%                       inputs to function numtrials, numconds and numbins
%           ?, JLA - removed VAF muscle per bin calculations
%           3/3/15, JLA - added condition to calculate VAF either
%                         before/after unit variance scaling

% stdev = std(V,0,2);
% 
% nmus = size(V,1);
% 
% % make unit variance
% Vnew = diag(1./stdev)*V;

if nargin < 3, VAFtype = 1; end

parfor r=1:numSyn
    
    [W,C,err]=nnmf_Jessica(V,r,1);

    emgRecon=W*C;

    synergyOut(r).C = C;
    synergyOut(r).W = W;
    synergyOut(r).emgRecon = emgRecon;
    synergyOut(r).emgOrig = V;
    synergyOut(r).err = err;
    
    % Calculate and save VAF
    stdev = std(V'); % added JA, 5/19/15    
    % modified JA, 3/3/16
    if VAFtype == 2
        [VAFcond, VAFmus, VAF]=funur(diag(1./stdev)*V,diag(1./stdev)*W,C); % added JA, 5/19/15
    elseif VAFtype == 1
        [VAFcond, VAFmus, VAF]=funur(V,W,C);
    end
    
    synergyOut(r).VAFcond = VAFcond;
    synergyOut(r).VAFmus = VAFmus;
    synergyOut(r).VAF = VAF;
    

    % Calculate and save r2 values
    [Rcond, Rmus, Rsqr]=funR(V,W,C);
    
    synergyOut(r).Rcond = Rcond;
    synergyOut(r).Rmus = Rmus;
    synergyOut(r).Rsqr = Rsqr;
    
    
%     C=reshape(C',numtrials,numconds,numbins,r);    
%     Data=reshape(V',numtrials,numconds,numbins,nmus);
%     
%     % VAFs and Rs for each direction forming muscle tuning curves
%     % VAF for one given direction and time bin, all are stored in
%     % VAFconds2,VAFmus2, etc.
%     % This information goes in the pannels with direction,bins,muscles
%     % Size of VAFconds2 = [ntrials x nconds x nbins] and VAFmus2 = [nmus x nconds x nbins]
%     temp=[];
%     temp1=[];
%     for b=1:numbins
%         for n=1:numconds
%             temp=squeeze(C(:,n,b,:));
%             temp1=squeeze(Data(:,n,b,:));
%             Co=temp;
%             Da=temp1;
%             if numtrials==1
%                 [VAFconds2(:,n,b), VAFmus2(:,n,b), VAFs2(n,b)]=funur(Da,W,Co);
%                 [Rconds2(:,n,b), Rmus2(:,n,b), Rsqrs2(n,b)]=funR(Da,W,Co);
%             else
%                 [VAFconds2(:,n,b), VAFmus2(:,n,b), VAFs2(n,b)]=funur(Da',W,Co');
%                 [Rconds2(:,n,b), Rmus2(:,n,b), Rsqrs2(n,b)]=funR(Da',W,Co');
%             end
%         end
%     end
%     
%     synergyOut(r).VAFmusByCond = VAFmus2;
%     
%     
%     %VAFs and Rs for the tuning curves formed by multiple trials
%     %VAF for one given time bin, all are stored in VAFconds1,VAFmus1, etc.
%     %VAFconds1 = [60 x nbins] (not useful) and VAFmus1 = [nmus x nbins]
%     
%     temp=[];
%     temp1=[];
%     for b=1:numbins
%         temp=squeeze(C(:,:,b,:));
%         temp1=squeeze(Data(:,:,b,:));
%         Co=reshape(temp,numtrials*numconds,r);
%         Da=reshape(temp1,numtrials*numconds,nmus);
%         [~, VAFmus1(:,b), VAFs1(b)]=funur(Da',W,Co');
% %         [Rconds1(:,b), Rmus1(:,b), Rsqrs1(b)]=funR(Da',W,Co');
%     end
%     
%     synergyOut(r).VAFmusByBin = VAFmus1;
%     synergyOut(r).VAFbin = VAFs1;
    C = [];
    W = [];
    emgRecon = [];
    err = [];
    VAFcond = [];
    VAFmus = [];
    VAF = [];
    Rcond = [];
    Rmus = [];
    Rsqr = []; 
%     clear C W emgRecon err VAFcond VAFmus VAF Rcond Rmus Rsqr
    
end