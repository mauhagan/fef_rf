
%
% The software definition file defines support for all the experimental
% software relevant to the data acquisition GUI.
%
%

%
% This definition files defines the paths for Reggie, when  
% 4 electrodes are used.
%

% global experiment

experiment.software.base_path = '/mnt/raid/Hobbes/';
experiment.software.mex_path_suffix = '/C/DoubleMT/2x4/mex/';
experiment.software.mex_path = [experiment.software.base_path experiment.software.mex_path_suffix];
experiment.software.release_base_path = '/mnt/raid/Hobbes/';
experiment.software.labview.path =  '/Labview/';
experiment.software.display_client.path = '/Display Client/';
experiment.software.matlab.path = '/matlab/';
experiment.software.GUI.path = '/GUI';
experiment.software.mat.path =  '/mat/';
experiment.software.m.path = '/m/';
experiment.software.C.path = '/C/';


%TDT
experiment.software.engine.terminal_title = '2x4 broker';
experiment.software.engine.path = '/mnt/raid/Hobbes/C/DoubleMT/2x4/broker';
experiment.software.module(1).terminal_title = 'module analogstate';
experiment.software.module(1).path = '/mnt/raid/Hobbes/C/DoubleMT/2x4/module_analogstate';


experiment.software.packages = {'labview','matlab','display_client','GUI','m','C'};

experiment.software.analyze = {'procDatMonash','procSpMuLfp','procDayEventsMonash'};

%%%
% Variables below are not used as yet
% For when the broker doesn't provide the functionality
%%%
experiment.software.electrodes(1).raw.available = 1;
experiment.software.electrodes(1).mu.available =  1;
experiment.software.electrodes(1).lfp.available =  1;
experiment.software.electrodes(1).spike.available =  1;
experiment.software.electrodes(2).raw.available = 1;
experiment.software.electrodes(2).mu.available =  1;
experiment.software.electrodes(2).lfp.available =  1;
experiment.software.electrodes(2).spike.available =  1;
experiment.software.electrodes(2).raw.available = 1;
experiment.software.electrodes(2).mu.available =  1;
experiment.software.electrodes(2).lfp.available =  1;
experiment.software.electrodes(2).spike.available =  1;
experiment.software.electrodes(2).raw.available = 1;
experiment.software.electrodes(2).mu.available =  1;
experiment.software.electrodes(2).lfp.available =  1;
experiment.software.electrodes(2).spike.available =  1;
experiment.software.behavior.state.available =  1;
experiment.software.behavior.eye.available =  1;
experiment.software.behavior.hand.available =  1;
experiment.software.behavior.stimulation.available = 1; 
experiment.software.set_filters.available =  1;
experiment.software.set_output_volume.available = 1;
experiment.software.set_output_channel.available =  1;
experiment.software.send_parallel.available =  1;
experiment.software.receive_parallel.available =  1;
experiment.software.control_recording.available = 1;
experiment.software.microdrive(1).available =  1;
experiment.software.microdrive(1).type =  1;
experiment.software.microdrive(1).control.available = 1;
experiment.software.microdrive(2).available =  1;
experiment.software.microdrive(2).type =  1;
experiment.software.microdrive(2).control.available = 1;
