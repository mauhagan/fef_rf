function LatResults = LatencyAnalysis_utah(monkeyname, stimsize)

% parallize slurm to run all monkeys

% Latecny Analysis
% for all monkeys

% only run for white stimuli - black stimuli may have longer latecy? or run
% black separately

% only run for trials inside the receptive field

% need to make a spike matrix ntrials x 400 ms, spike count in each ms

%
% addpath(genpath('/projects/zb33/analyze'))
% addpath(genpath('/scratch/zb33/Data/'))

Task = 'RFMap';
nch = 96;
tract = '1';
array = 'FEF';

if strcmp(monkeyname,'CJ197')
    array = 'V1';
end

% stimsizes = [2 4 6]; % focus on these trials for now

% load  '/Users/mhagan/Dropbox/Les Figures/LesData/RFAnalysis/RFResults.mat'

% Monkeys = {'F1594','M1323'}; % add in Liz's controls
%
% for iMonkey = 1:length(Monkeys)
%     monkeyname = Monkeys{iMonkey};


% for istim = 1:length(stimsizes)
clear stimtrials
%     stimsize = stimsizes(istim);
stimname = ['stim' num2str(stimsize)];

if exist(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFlatency/LatencyResults_' monkeyname '_' stimname '.mat'],'file')
    load(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFlatency/LatencyResults_' monkeyname '_' stimname '.mat'],'LatResults')
    startchan = LatResults.(monkeyname).(stimname).completed + 1;
else
    startchan = 1;
end

if startchan <= nch
    
    trials = dbSelectTrials(tract);
    tasktrials = trials(strcmp({trials.Task}, Task));
    
    stimtrials = tasktrials([tasktrials.StimSize] == stimsize);
    
    whitetrials = stimtrials([stimtrials.Phase] == 90);
    blacktrials = stimtrials([stimtrials.Phase] ~= 90);
    
    load(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFanalysis/RFResults_' monkeyname '_' stimname '_p01.mat'],'RFResults')
    
    
    allRows = unique([stimtrials.Ydov]);
    allRows = sort(allRows,'descend'); % so that top of screen is in right place
    allCols = unique([stimtrials.Xdov]);
    
    nRows = numel(allRows);
    nCols = numel(allCols);
    
    % for CJ182
    bad_ch = [2 4 12 17 18 30 45 51 53 89];
    
    stimdur = 400; %calculate the latency within 400 ms
    
    for ich = startchan:nch
        clu = 1; % all multiunit
        
        if ~strcmp(monkeyname,'CJ182') || (strcmp(monkeyname,'CJ182') && ~ismember(ich,bad_ch))
            
            % is there a receptive field?
            
            RFMap = sq(RFResults.(monkeyname).(stimname).rfdata.RFbinary(ich,:,:));
            
            if sum(RFMap(:))
                
                % white trials
                if ~isempty(whitetrials)
                    
                    % get the trial indicies for the RF
                    rf_trial_inds = [];
                    for iCol = 1:nCols
                        for iRow = 1:nRows
                            if RFMap(iRow,iCol)
                                Ycoord = allRows(iRow);
                                Xcoord = allCols(iCol);
                                rf_inds = find([whitetrials.Xcoord] == Xcoord & [whitetrials.Ycoord] == Ycoord);
                                rf_trial_inds = [rf_trial_inds, rf_inds];
                            end
                        end
                    end
                    
                    RFtrials = whitetrials(rf_trial_inds);
                    
                    %calculate the spike times for the RF trials
                    RFspikes = trialSpike(RFtrials,array,ich, 1,clu,'StimOn',[0,stimdur]);
                    
                    rate = psth(RFspikes,[0,stimdur],10);
                    %         spikecount = round(rate);
                    spikes = [];
                    for itr = 1:length(RFspikes)
                        x = RFspikes{itr}';
                        spikes = [spikes x];
                    end
                    
                    e = 0:1:stimdur; % 1 ms bins
                    spikecount = histcounts(spikes,e);
                    
                    % spikeCount =  Number of spikes in each millisecond. The first entry in
                    %               this vector is considered time zero.
                    % minTheta   = The first index in spikeCount at which a response may have
                    %               started.
                    % maxTheta   = The last entry in spikeCount at which a response may have
                    %               started.
                    % maxKappa   =  Last index at which the peak of the response could be
                    %
                    % delta       =  The minimum time between response onset and sustained
                    %                   onset.
                    % responseSign = Set to 1 to detect only increases in firing
                    %                Set to -1 to detect only decreases in firing
                    %                Set to 0 to detect both increases and decreases in firing.
                    % graphics     = Set to true to show a graph of the procedure.
                    
                    [latencyEstimate,peakEstimate,pPoisson,baseRate,responseRate] = friedmanpriebe(spikecount','minTheta',10,'maxTheta',stimdur-1);
                    
                    % latency results
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.RFspikes{ich} = RFspikes;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.Spikecount(ich,:) = spikecount;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.PSTH(ich,:) = rate;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.Latency(ich) = latencyEstimate;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.Peak(ich) = peakEstimate;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.pval(ich) = pPoisson;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.baseRate(ich) = baseRate;
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.responseRate(ich) = responseRate;
                else
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.RFspikes{ich} = cell(1);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.Spikecount(ich,:) = nan(1,stimdur);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.PSTH(ich,:) = nan(1,stimdur+1);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.Latency(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.Peak(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.pval(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.baseRate(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.whitetrials.responseRate(ich) = nan(1);
                end
                
                % black trials
                if ~isempty(blacktrials)
                    
                    % get the trial indicies for the RF
                    rf_trial_inds = [];
                    for iCol = 1:nCols
                        for iRow = 1:nRows
                            if RFMap(iRow,iCol)
                                Ycoord = allRows(iRow);
                                Xcoord = allCols(iCol);
                                rf_inds = find([blacktrials.Xcoord] == Xcoord & [blacktrials.Ycoord] == Ycoord);
                                rf_trial_inds = [rf_trial_inds, rf_inds];
                            end
                        end
                    end
                    
                    RFtrials = blacktrials(rf_trial_inds);
                    
                    %calculate the spike times for the RF trials
                    RFspikes = trialSpike(RFtrials,array,ich, 1,clu,'StimOn',[0,stimdur]);
                    
                    rate = psth(RFspikes,[0,stimdur],10);
                    %         spikecount = round(rate);
                    spikes = [];
                    for itr = 1:length(RFspikes)
                        x = RFspikes{itr}';
                        spikes = [spikes x];
                    end
                    
                    e = 0:1:stimdur; % 1 ms bins
                    spikecount = histcounts(spikes,e);
                    
                    [latencyEstimate,peakEstimate,pPoisson,baseRate,responseRate] = friedmanpriebe(spikecount','minTheta',10,'maxTheta',stimdur-1);
                    
                    % latency results
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.RFspikes{ich} = RFspikes;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.Spikecount(ich,:) = spikecount;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.PSTH(ich,:) = rate;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.Latency(ich) = latencyEstimate;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.Peak(ich) = peakEstimate;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.pval(ich) = pPoisson;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.baseRate(ich) = baseRate;
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.responseRate(ich) = responseRate;
                else
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.RFspikes{ich} = cell(1);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.Spikecount(ich,:) = nan(1,stimdur);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.PSTH(ich,:) = nan(1,stimdur+1);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.Latency(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.Peak(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.pval(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.baseRate(ich) = nan(1);
                    LatResults.(monkeyname).(stimname).latdata.blacktrials.responseRate(ich) = nan(1);
                end
                
            else
                
                LatResults.(monkeyname).(stimname).latdata.whitetrials.RFspikes{ich} = cell(1);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.Spikecount(ich,:) = nan(1,stimdur);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.PSTH(ich,:) = nan(1,stimdur+1);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.Latency(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.Peak(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.pval(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.baseRate(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.whitetrials.responseRate(ich) = nan(1);
                
                LatResults.(monkeyname).(stimname).latdata.blacktrials.RFspikes{ich} = cell(1);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.Spikecount(ich,:) = nan(1,stimdur);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.PSTH(ich,:) = nan(1,stimdur+1);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.Latency(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.Peak(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.pval(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.baseRate(ich) = nan(1);
                LatResults.(monkeyname).(stimname).latdata.blacktrials.responseRate(ich) = nan(1);
            end
        else
            LatResults.(monkeyname).(stimname).latdata.whitetrials.RFspikes{ich} = cell(1);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.Spikecount(ich,:) = nan(1,stimdur);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.PSTH(ich,:) = nan(1,stimdur+1);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.Latency(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.Peak(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.pval(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.baseRate(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.whitetrials.responseRate(ich) = nan(1);
            
            LatResults.(monkeyname).(stimname).latdata.blacktrials.RFspikes{ich} = cell(1);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.Spikecount(ich,:) = nan(1,stimdur);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.PSTH(ich,:) = nan(1,stimdur+1);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.Latency(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.Peak(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.pval(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.baseRate(ich) = nan(1);
            LatResults.(monkeyname).(stimname).latdata.blacktrials.responseRate(ich) = nan(1);
        end
        
        LatResults.(monkeyname).(stimname).completed = ich;
        LatResults.(monkeyname).(stimname).xLabels = allCols;
        LatResults.(monkeyname).(stimname).yLabels = allRows;
        
        save(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFlatency/LatecyResults_' monkeyname '_' stimname '.mat'],'LatResults','-v7.3')
        
        disp(['completed! ' monkeyname ' ' stimname ' ch ' num2str(ich)])
    end
end
