function RFResults = RFAnalysis_utah(monkeyname, stimsize)

% parallize slurm to run all monkeys

% RF Analysis
% for all monkeys
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

stimsizes = [2 4 6]; % focus on these trials for now

% load  '/Users/mhagan/Dropbox/Les Figures/LesData/RFAnalysis/RFResults.mat'

% Monkeys = {'F1594','M1323'}; % add in Liz's controls
%
% for iMonkey = 1:length(Monkeys)
%     monkeyname = Monkeys{iMonkey};

% eval(['lesrf_' monkeyname '_massive'])
trials = dbSelectTrials(tract);
tasktrials = trials(strcmp({trials.Task}, Task));


% for istim = 1:length(stimsizes)
clear stimtrials
%     stimsize = stimsizes(istim);
stimname = ['stim' num2str(stimsize)];

if exist(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFanalysis/RFResults_' monkeyname '_' stimname '_p01.mat'])
    load(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFanalysis/RFResults_' monkeyname '_' stimname '_p01.mat'])
    startchan = RFResults.(monkeyname).(stimname).completed + 1;
else
    startchan = 1;
end

if startchan <= nch
    stimtrials = tasktrials([tasktrials.StimSize] == stimsize);
    stimdur = stimtrials(1).stimDur;
    
    if stimdur < 200, stimdur = 400; end  % for fef use larger time window because of black and white flashes!
    
    allRecs = unique({stimtrials.Rec});
    
    % need to linearize number of rows and columns since there
    % are two screen heights
    
    allRows = unique([stimtrials.Ydov]);
    allRows = sort(allRows,'descend'); % so that top of screen is in right place
    allCols = unique([stimtrials.Xdov]);
    
%     newRows = allRows(1):-stimsize:allRows(end);
%     newCols = allCols(1):stimsize:allCols(end);

% fef studies - stim size is bigger than stepsize!! dont need new rows
% anyway, only one screen height

    newRows = allRows;
    newCols = allCols;
    
    nRows = length(newRows);
    nCols = length(newCols);
    
    % for CJ182
    bad_ch = [2 4 12 17 18 30 45 51 53 89];
    
    for ich = startchan:nch
        clu = 1; % all multiunit
        
        if ~strcmp(monkeyname,'CJ182') || (strcmp(monkeyname,'CJ182') && ~ismember(ich,bad_ch))
        
        taskRates = trialRate(stimtrials,array,ich, 1,clu,'StimOn',[0,stimdur]);
        spontRates = trialRate(stimtrials,array,ich, 1,clu,'StimOn',[-stimdur,0]);
        %
        allSpikes = cell(nRows,nCols);
        nTrials = nan(nRows,nCols);
        meanSpikes = nan(nRows,nCols);
        stdSpikes = nan(nRows,nCols);
        
        for iCol = 1:nCols
            colind = [stimtrials.Xdov] > allCols(iCol)-stimsize & [stimtrials.Xdov] <= allCols(iCol);
            for iRow = 1:nRows
                rowind = [stimtrials.Ydov] > allRows(iRow)-stimsize & [stimtrials.Ydov] <= allRows(iRow);
                allSpikes{iRow,iCol} = taskRates(rowind & colind)';
                nTrials(iRow,iCol) = numel(taskRates(rowind & colind)');
                meanSpikes(iRow,iCol) = nanmean(taskRates(rowind & colind));
                stdSpikes(iRow,iCol) = nanstd(taskRates(rowind & colind));
            end
        end
        
        spontSpikes = spontRates;
        
        spontMeanSpikeRate = mean(spontSpikes);
        spontStdSpikeRate  = std(spontSpikes);
        
        %             %%
        %             xCentroidMin = round(min(newCols));
        %             xCentroidMax = round(max(newCols));
        %             yCentroidMin = round(min(newRows));
        %             yCentroidMax = round(max(newRows));
        %
        %             xLabels = round(linspace(xCentroidMin, xCentroidMax, nCols));
        %             yLabels = round(linspace(yCentroidMax, yCentroidMin, nRows));
        
        %% 3. cluster analysis
        % determine the min cluster size for the data using random permutations
        % switched to only using stimulus times - spont firing rate
        % too variable. everything is significant.
        
        nPerm = 1e3;
        cluster = zeros(1,nPerm);
        alldata = taskRates;
        ntr2 = size(alldata);
        for iPerm = 1:nPerm
            %disp(['Perm #' num2str(iPerm)])
            tmpPvals = zeros(nRows,nCols);
            for iCol = 1:nCols
                for iRow = 1:nRows
                    ntr1 = nTrials(iRow,iCol);
                    tmpdata = shuffle(alldata);
                    tmp1 = tmpdata(1:ntr1); tmp1 = tmp1(~isnan(tmp1));
                    tmp2 = tmpdata(ntr1+1:ntr1+1+ntr2); tmp2 = tmp2(~isnan(tmp2));
                    if ~isempty(tmp1) && ~isempty(tmp2)
                        tmpPvals(iRow,iCol) = ranksum(tmp1, tmp2,'tail','right');
                    else tmpPvals(iRow,iCol) = 1;
                    end
                end
            end
            tmp2 = zeros(size(tmpPvals));
            tmp2(tmpPvals<0.01) = 1;
            [L,n] = bwlabel(tmp2);
            bins = hist(reshape(L,size(L,1)*size(L,2),1),n+1);
            if length(bins)>1
                cluster(iPerm) = max(bins(2:end));
            else
                cluster(iPerm) = 0;
            end
        end
        cluster = sort(cluster,'descend');
        mincluster = cluster(10); % 1st percentile for cluster, p<0.01
        
        
        % find the actual clusters in the real data
        Pvals = nan(nRows,nCols);
        for iRow = 1:nRows
            for iCol = 1:nCols
                if sum(~isnan(allSpikes{iRow,iCol})) && sum(~isnan(alldata))
                    Pvals(iRow,iCol) = ranksum(allSpikes{iRow,iCol},alldata,'tail','right');
                else Pvals(iRow,iCol) = 1;
                end
            end
        end
        tmp = zeros(size(Pvals));
        tmp(Pvals<0.01) = 1;
        [L,n] = bwlabel(tmp);
        b = 0:1:100;
        [bins,bind] = hist(reshape(L,size(L,1)*size(L,2),1),n+1,b);
        bind = round(bind);
        %bins(1) = 0;
        sigclusters = find(bins>mincluster);
        % first bin is non sig clusters, so this is always the biggest, so make it go away
        if length(sigclusters>1), sigclusters = sigclusters(2:end); else sigclusters = []; end
        
        RFmap = zeros(size(Pvals)); % main cluster (RF)
        Sigmap = zeros(size(Pvals)); % all sig clusters
        if ~isempty(sigclusters)
            RFcluster = bind(find(bins == max(bins(sigclusters)))); % call the RF the biggest sig cluster, we'll keep the others just in case;
            for iCol = 1:nCols
                for iRow = 1:nRows
                    if L(iRow,iCol) == RFcluster, RFmap(iRow,iCol) = 1; end
                    if ismember(L(iRow,iCol),bind(sigclusters)), Sigmap(iRow,iCol) = 1; end
                end
            end
        end
        
        %% RF properties
        if sum(RFmap(:)) > 0
            
            %center of mass & peak firing rate
            RFMapSpikes = nan(nRows,nCols);
            RFMapSpikes(RFmap == 1) = meanSpikes(RFmap == 1);
            %peak
            peakFR = max(RFMapSpikes(:));
            [px,py] = find(RFMapSpikes == peakFR);
            if numel(px) > 1
                peak = [mean(newCols(py)), mean(newRows(px))];
            else
                peak = [newCols(py), newRows(px)];
            end
            
            %center of mass
            [rows, cols] = ndgrid(newRows, newCols);
            rowcenter = sum(rows(RFmap==1) .* meanSpikes(RFmap==1)) / sum(meanSpikes(RFmap==1));
            colcenter = sum(cols(RFmap==1) .* meanSpikes(RFmap==1)) / sum(meanSpikes(RFmap==1));
            
            com = [colcenter, rowcenter];
            
            % center of mass with normalized firing rates
            normSpikes = meanSpikes./peakFR;
            rowcenter = sum(rows(RFmap==1) .* normSpikes(RFmap==1)) / sum(normSpikes(RFmap==1));
            colcenter = sum(cols(RFmap==1) .* normSpikes(RFmap==1)) / sum(normSpikes(RFmap==1));
            
            norm_com = [colcenter, rowcenter];
            
            
            % distortion
            dist = sqrt((peak(2)-com(2))^2 + (peak(1)-com(1))^2);
            norm_dist = sqrt((peak(2)-norm_com(2))^2 + (peak(1)-norm_com(1))^2);
            
            RFarea = sum(RFmap(:))*stimsize; % area of flashes
            RFdiam = 2*sqrt(RFarea/pi); % rough diameter
        else
            peakFR = nan(1);
            peak = nan(1,2);
            com = nan(1,2);
            norm_com = nan(1,2);
            dist = nan(1);
            norm_dist = nan(1);
            RFarea = nan(1);
            RFdiam = nan(1);
        end
        
        
        %%   put it all together
        RFResults.(monkeyname).(stimname).completed = ich;
        RFResults.(monkeyname).(stimname).oldX = allCols;
        RFResults.(monkeyname).(stimname).oldY = allRows;
        RFResults.(monkeyname).(stimname).xLabels = newCols;
        RFResults.(monkeyname).(stimname).yLabels = newRows;
        % firing rate data
        RFResults.(monkeyname).(stimname).chdata(ich).allSpikes = allSpikes;
        RFResults.(monkeyname).(stimname).chdata(ich).meanSpikes = meanSpikes;
        RFResults.(monkeyname).(stimname).chdata(ich).stdSpikes = stdSpikes;
        RFResults.(monkeyname).(stimname).chdata(ich).spontSpikes = spontSpikes;
        % cluster results
        RFResults.(monkeyname).(stimname).rfdata.RFbinary(ich,:,:) = RFmap;
        RFResults.(monkeyname).(stimname).rfdata.Sigbinary(ich,:,:) = Sigmap;
        % RF properites
        RFResults.(monkeyname).(stimname).rfdata.peakFR(ich) = peakFR;
        RFResults.(monkeyname).(stimname).rfdata.area(ich) = RFarea;
        RFResults.(monkeyname).(stimname).rfdata.diam(ich) = RFdiam;
        RFResults.(monkeyname).(stimname).rfdata.dist(ich) = dist;
        RFResults.(monkeyname).(stimname).rfdata.ndist(ich) = norm_dist;
        RFResults.(monkeyname).(stimname).rfdata.peak(ich,:,:) = peak;
        RFResults.(monkeyname).(stimname).rfdata.com(ich,:,:) = com;
        RFResults.(monkeyname).(stimname).rfdata.ncom(ich,:,:) = norm_com;
        disp(['completed! ' monkeyname ' ' stimname ' ch ' num2str(ich)])
        save(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFanalysis/RFResults_' monkeyname '_' stimname '_p01.mat'],'RFResults','-v7.3')
        else
                 %%   put it all together
        RFResults.(monkeyname).(stimname).completed = ich;
        RFResults.(monkeyname).(stimname).oldX = allCols;
        RFResults.(monkeyname).(stimname).oldY = allRows;
        RFResults.(monkeyname).(stimname).xLabels = newCols;
        RFResults.(monkeyname).(stimname).yLabels = newRows;
        % firing rate data
        RFResults.(monkeyname).(stimname).chdata(ich).allSpikes = nan(nRows,nCols);
        RFResults.(monkeyname).(stimname).chdata(ich).meanSpikes = nan(nRows,nCols);
        RFResults.(monkeyname).(stimname).chdata(ich).stdSpikes = nan(nRows,nCols);
        RFResults.(monkeyname).(stimname).chdata(ich).spontSpikes = nan(nRows,nCols);
        % cluster results
        RFResults.(monkeyname).(stimname).rfdata.RFbinary(ich,:,:) = nan(nRows,nCols);
        RFResults.(monkeyname).(stimname).rfdata.Sigbinary(ich,:,:) = nan(nRows,nCols);
        % RF properites
        RFResults.(monkeyname).(stimname).rfdata.peakFR(ich) = nan(1);
        RFResults.(monkeyname).(stimname).rfdata.area(ich) = nan(1);
        RFResults.(monkeyname).(stimname).rfdata.diam(ich) = nan(1);
        RFResults.(monkeyname).(stimname).rfdata.dist(ich) = nan(1);
        RFResults.(monkeyname).(stimname).rfdata.ndist(ich) = nan(1);
        RFResults.(monkeyname).(stimname).rfdata.peak(ich,:,:) = nan(1,2);
        RFResults.(monkeyname).(stimname).rfdata.com(ich,:,:) = nan(1,2);
        RFResults.(monkeyname).(stimname).rfdata.ncom(ich,:,:) = nan(1,2);
        disp(['completed! ' monkeyname ' ' stimname ' ch ' num2str(ich)])
        save(['/home/mhagan/zb33_scratch/Data/fefrf/analysis/RFanalysis/RFResults_' monkeyname '_' stimname '_p01.mat'],'RFResults','-v7.3')
    end
end
%     end
end