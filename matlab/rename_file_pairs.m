% moving files to new directories for spikesorting
%shared drive path
% gvfs-mount smb://storage.erc.monash.edu.au/shares/R-MNHS-Lui-lab/
% ~/.gvfs/smb-share:server=storage.erc.monash.edu.au,share=shares/R-ENG-Wong-lab/Shared/CJ188/og
%% for MT files anyway...
monkeyname = 'CJ182';

cd('/home/mhagan/zb33_scratch/Data/fefrf/matlab')
eval(['fefrf_' monkeyname '_massive'])

global MONKEYDIR
super = 1;
if(super)
    filepath = '/scratch/zb33/Data/fefrf/';
else
    filepath = '/Users/ywon0005/Documents/Data/';
end


%%
if strcmp(monkeyname, 'CJ170')
    %     monkeyname = 'CJ177';
    pairs =    {'1','0041','028';
        '1','0007','009'};
    
elseif strcmp(monkeyname, 'CJ174')
    %     monkeyname = 'CJ179';
    pairs =    {'1','0010','012';
        '1','0016','017'};
    
    elseif strcmp(monkeyname, 'CJ176')
    pairs =    {'1','0028','029'};
    
    elseif strcmp(monkeyname, 'CJ182')
    pairs =    {'1','0027','011';
        '1','0038','020'};
end
% MONKEYDIR = [filepath '/' monkeyname];


%%
% copy files
cd([filepath monkeyname '/'])

areas = {'FEF'};
% loop thorugh each

recs = size(pairs,1);

for iRec = 1:recs
    day = pairs{iRec,1};
    if ~isfolder(day)
        mkdir(day)
    end
    for iArea = 1:size(areas,2)
        rec = pairs{iRec,3};
        if ~isfolder([day '/' rec])
            mkdir([day '/' rec])
        end
        
        disp(['now working on ' areas{iArea} ' rec ' pairs{iRec,3}])
        
%         % copy files
%         if (isfile([day '/' rec '/rec' rec '.' areas{iArea} '.ns6']))
%             disp(['file already exists: ' day '/' rec '/rec' rec '.' areas{iArea} '.ns6'])
%         else
%             copyfile(['og/datafile' pairs{iRec,3} '.ns6'], [day '/' rec '/rec' rec '.' areas{iArea} '.ns6']);
%         end
        
%         if (isfile([day '/' rec '/rec' rec '.' areas{iArea} '.ns3']))
%             disp(['file already exists: ' day '/' rec '/rec' rec '.' areas{iArea} '.ns3'])
%         else
%             copyfile(['og/datafile' pairs{iRec,3} '.ns3'], [day '/' rec '/rec' rec '.' areas{iArea} '.ns3']);
%         end
        
        if (isfile([day '/' rec '/rec' rec '.' areas{iArea} '.nev']))
            disp(['file already exists: ' day '/' rec '/rec' rec '.' areas{iArea} '.nev'])
        else
           if isfile([filepath 'og/' monkeyname '/' monkeyname '_datafile' pairs{iRec,3} '.nev'])
               copyfile([filepath 'og/' monkeyname '/' monkeyname '_datafile' pairs{iRec,3} '.nev'], [day '/' rec '/rec' rec '.' areas{iArea} '.nev']);
           end
        end
        
        
    end
    % need to check for the different task types
    if (isfile([day '/' rec '/rec' rec '.trials.mat']))
        disp(['file already exists: ' day '/' rec '/rec' rec '.trials.mat'])
    elseif (isfile([filepath 'og/' monkeyname '/' monkeyname '_RFMap' pairs{iRec,2} '.mat'])) 
        copyfile([filepath 'og/' monkeyname '/' monkeyname '_RFMap' pairs{iRec,2} '.mat'], [day '/' rec '/rec' rec '.trials.mat']);
    elseif (isfile([filepath 'og/' monkeyname '/' monkeyname '_RFmap_' pairs{iRec,2} '.mat'])) 
        copyfile([filepath 'og/' monkeyname '/' monkeyname '_RFmap_' pairs{iRec,2} '.mat'], [day '/' rec '/rec' rec '.trials.mat']);
    elseif (isfile([filepath 'og/' monkeyname '/' monkeyname '_PlaidRat_' pairs{iRec,2} '.mat'])) 
        copyfile([filepath 'og/' monkeyname '/' monkeyname '_PlaidRat_' pairs{iRec,2} '.mat'], [day '/' rec '/rec' rec '.trials.mat']);
    elseif (isfile([filepath 'og/' monkeyname '/' monkeyname '_PlaidRat' pairs{iRec,2} '.mat'])) 
        copyfile([filepath 'og/' monkeyname '/' monkeyname '_PlaidRat' pairs{iRec,2} '.mat'], [day '/' rec '/rec' rec '.trials.mat']);
   
    end
    
    if (isfile([day '/' rec '/rec' rec '.experiment.mat']))
        disp(['file already exists: ' day '/' rec '/rec' rec '.experiment.mat'])
    else
        if(isfile([day '/mat/experiment.mat']))
            copyfile([day '/mat/experiment.mat'], [day '/' rec '/rec' rec '.experiment.mat']);
        else
            
            % make the experiment.mat
            cd([MONKEYDIR '/' day]);
            if ~isfolder('mat')
                mkdir('mat')
            end
            addpath([filepath  '/experiment'])
            %             addpath('/Users/ywon0005/Documents/Matlab/experiment/')
            
            cd([MONKEYDIR '/' day '/mat/']);
            if ~isfile('experiment.mat')
                acquire_gui_definition_file
                hardware_definition_file
                histograms_definition_file
                recording_definition_file
                software_definition_file
                save('experiment.mat','experiment')
            end
            cd([filepath monkeyname '/'])
            copyfile([day '/mat/experiment.mat'], [day '/' rec '/rec' rec '.experiment.mat']);
        end
    end
end

%% fixing experiment files
% forgot to change drive name to FEF
% for some reason we only recorded from 95 channels...
monkeyname = 'CJ182';
day = '1';
cd([filepath monkeyname '/'])
clear experiment

cd('/scratch/zb33/Data/fefrf/experiment/');
% new experiment file
hardware_definition_file
acquire_gui_definition_file
histograms_definition_file
recording_definition_file
software_definition_file

cd([MONKEYDIR '/' day '/mat/']);
save('experiment.mat','experiment')

areas = {'MT'};
% loop thorugh each

recs = size(pairs,1);

for iRec = 1:recs
    day = pairs{iRec,1};
    
    for iArea = 1:size(areas,2)
        rec = pairs{iRec,3};       
        disp(['now working on ' areas{iArea} ' rec ' pairs{iRec,3}])
            cd([filepath monkeyname '/'])
            copyfile([day '/mat/experiment.mat'], [day '/' rec '/rec' rec '.experiment.mat']);
    end
end
