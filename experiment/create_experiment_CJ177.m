%% make experiment file
clear all
close all
monkeyname = 'CJ156';
tract = '1';
%cd(['/Users/mhagan/Documents/MTutahData/mt_utah'])
cd('C:\Users\wongy4\Mo_Google_Drive\mt_utah\');
% mt_utah_CJ156

acquire_gui_definition_file
hardware_definition_file
histograms_definition_file
recording_definition_file
software_definition_file

%save(['/Users/mhagan/Documents/MTutahData/' monkeyname '/' tract '/mat/experiment.mat'], 'experiment')
save(['C:\Users\wongy4\Mo_Google_Drive\mt_utah\' monkeyname '\' tract '\mat\experiment.mat'], 'experiment')