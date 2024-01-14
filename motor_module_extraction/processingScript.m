clear all
close all
Folder = 'E:\wvu\Allen works\papers to read\muscle synergy\synergyFiles\synergyFiles\';

subjects = {'YASub01','YASub02','YASub03'};
condition = {'VR_0_125','VR_20_125','VR_35_125','VR_50_125'};
indices = 1:60000;

% general information about the data
SampleRate = 1000;
HP = 35; % highpass filter of 35 Hz
LP = 10; % lowpas filter of 10 Hz
% EMGnames = {'TA','MGAS','PERO','SOL','MH','VLAT','GMED','EXOB','INOB','ERSP'};
EMGnames = {'TA','MGAS','SOL','PERO','VLAT','MH','GMED'};

for k=1:length(subjects)
    subj = subjects{k};
    
%     try
        % load data
        load([Folder subj '.mat'])
        
        disp(['processing ' subj])
        
        for i = 1:length(condition)
            
            EMG = Output.(condition{i}).EMG{1}(indices,1:length(EMGnames));
            RHS = Output.(condition{i}).RHS{1}*10; % multiply by 10 because RHS based on markers, with 100Hz sample rate
            RHS = RHS(RHS>indices(1)) - (indices(1)-1);
            %--------------------
           
            %-----------------------------------
            % FILTER THE EMG
            nyquist_frequency = SampleRate/2;
            % High pass filter at HP Hz
            [filt_high_B,filt_high_A] = butter(3,HP/nyquist_frequency,'high');
            emg = filtfilt(filt_high_B, filt_high_A,EMG);
            % Demean and rectify
            emg_mean = repmat(mean(emg(1:SampleRate,:),1),size(emg,1),1);
            emg = abs(emg-emg_mean);
            % Low pass filter at LP Hz
            [filt_low_B,filt_low_A] = butter(3,LP/nyquist_frequency,'low');
            emg = filtfilt(filt_low_B, filt_low_A,emg);
            
            emgMatrix = emg;
            time = (1/SampleRate):(1/SampleRate):(length(emgMatrix)/SampleRate);
            
            %----------------------------
             for ff=1:length(EMGnames)
              subplot(2,1,1)
           
              plot(EMG(:,ff))
              title([subj EMGnames(ff)])
             
              subplot(2,1,2)
              plot(emg(:,ff))
              title(['Filtered' subj EMGnames(ff)])
             
            
                set(gca,'Units','normalized','FontWeight','normal'...
                    ,'FontSize',12,'FontName','Times');

                print([Folder 'muscleSynergies' '-' subj '-' condition{i} '-' EMGnames{ff} '_'],'-dpng')
             end
            %-----------------------------
            emgMatrixGC = [];
            % parse into gait cycles
            for gc = 1:length(RHS)-1
                temp = emgMatrix(RHS(gc):RHS(gc+1),:);
                oldTime = time(RHS(gc):RHS(gc+1));
                newTime = time(RHS(gc)):(time(RHS(gc+1))-time(RHS(gc)))/200:time(RHS(gc+1));
                temp = spline(oldTime,temp',newTime)';
                emgMatrixGC = [emgMatrixGC; temp];
                
                clear temp oldTime newTime
            end
            
            % store all matrices so that they can be normalized to max later
            allEMGmatrix(i).emgMatrixGC = emgMatrixGC;
            
            % find maximum value across all conditions for each muscle
            if i==1
                maxVal = max(allEMGmatrix(i).emgMatrixGC);
            else
                temp = max(allEMGmatrix(i).emgMatrixGC);
                maxVal = max([temp;maxVal]);
            end
            
            clear emgMatrix emgMatrixGC temp
            
        end
        
        
        % normalize to maximum value
        for i = 1:length(condition)
            for m=1:length(EMGnames)
                temp(:,m) = allEMGmatrix(i).emgMatrixGC(:,m)./maxVal(m);
            end
            allEMGmatrix(i).emgMatrixGC_norm = temp;
            clear temp
        end
        
        
        % extract synergies from each condition
        NBOOT = 100; CI = 95;
        chooseFlag = 2; % choose number of synergies based on 1:CI, 2:max of all
       
        for i = 1:length(condition)
            data = allEMGmatrix(i).emgMatrixGC_norm';
            disp('extracting muscle synergies')
            synergyOut = synergyExtraction_Jessica(data,length(EMGnames));
            disp('bootstrapping muscle synergies')
            synergyOut = bootstrap_VAF_Jessica(synergyOut, length(EMGnames), NBOOT);
            disp('choosing number of synergies')
            nsyn = choose_synergies_walking([Folder subj '-' condition{i}],synergyOut,CI,chooseFlag,EMGnames');
            save([Folder 'muscleSynergies' '-' subj '-' condition{i} '_' date],'nsyn','synergyOut','allEMGmatrix','maxVal','EMGnames');
            clear data synergyOut nsyn
            close all
        end
                       
%     catch
%         disp([subj ' not processed'])
%     end
%     
    close all
end
