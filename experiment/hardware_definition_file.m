%
% The hardware definition file defines support for all the experimental
% hardware relevant to the acquire data acquisition GUI
%
% This definition file is for the Double MT system.


% global experiment

%
% 
% Definitions to configure the microdrives
%
%

% MT 
%experiment.hardware.type(2) = 'Cerebus';
for i = 1:96
    experiment.hardware.microdrive(1).electrodes(i).motorid = 1;
    experiment.hardware.microdrive(1).electrodes(i).acquisitionid = i;
    experiment.hardware.microdrive(1).electrodes(i).channelid = i;
    experiment.hardware.microdrive(1).electrodes(i).type = 'utah';
    experiment.hardware.microdrive(1).electrodes(i).label = ['mt' num2str(i)];
    experiment.hardware.microdrive(1).electrodes(i).gain = 10000;
    experiment.hardware.microdrive(1).electrodes(i).numContacts = 1;
end
experiment.hardware.microdrive(1).name = 'FEF';
experiment.hardware.microdrive(1).type = 'Cerebus';
experiment.hardware.microdrive(1).coordinate.x = -13;
experiment.hardware.microdrive(1).coordinate.y = 15;
experiment.hardware.microdrive(1).coordinate.pitch = -10;
experiment.hardware.microdrive(1).coordinate.yaw = 0;


%
%
% Definitions to configure the behaviour data acquired from the A/D
%
%

experiment.hardware.behavior(1).acquisitionid = 1;
experiment.hardware.behavior(1).channelid = 1;
experiment.hardware.behavior(1).label = 'State';
experiment.hardware.behavior(1).type = 'LabviewAnalogOut';

experiment.hardware.behavior(1).acquisitionid = 1;
experiment.hardware.behavior(1).channelid = 2;
experiment.hardware.behavior(1).label = 'EyeX';
experiment.hardware.behavior(1).type = 'LabviewAnalogOut';

experiment.hardware.behavior(1).acquisitionid = 1;
experiment.hardware.behavior(1).channelid = 3;
experiment.hardware.behavior(1).label = 'EyeY';
experiment.hardware.behavior(1).type = 'LabviewAnalogOut';

experiment.hardware.behavior(1).acquisitionid = 1;
experiment.hardware.behavior(1).channelid = 4;
experiment.hardware.behavior(1).label = 'HandX';
experiment.hardware.behavior(1).type = 'LabviewAnalogOut';

experiment.hardware.behavior(1).acquisitionid = 1;
experiment.hardware.behavior(1).channelid = 5;
experiment.hardware.behavior(1).label = 'HandY';
experiment.hardware.behavior(1).type = 'LabviewAnalogOut';

experiment.hardware.behavior(1).acquisitionid = 1;
experiment.hardware.behavior(1).channelid = 6;
experiment.hardware.behavior(1).label = 'Stim';
experiment.hardware.behavior(1).type = 'LabviewAnalogOut';


%
%
% Definition file to configure the acquisiton hardware
%
%

experiment.hardware.acquisition(1).samplingrate = 30000;
experiment.hardware.acquisition(1).data_type = 'Analog';
experiment.hardware.acquisition(1).num_channels = 8;
experiment.hardware.acquisition(1).data_format = 'short'; 
experiment.hardware.acquisition(1).type = 'Blackrock'; 
experiment.hardware.acquisition(1).voltage_range = [-5,5];
experiment.hardware.acquisition(1).AD_range = [0,4095];
experiment.hardware.acquisition(1).AD_neural_gain = 10/4096;
experiment.hardware.acquisition(1).channel(1).label = 'State';
experiment.hardware.acquisition(1).channel(1).filters = [0,10000];
experiment.hardware.acquisition(1).channel(2).label = 'EyeX';
experiment.hardware.acquisition(1).channel(2).filters = [0,10000];
experiment.hardware.acquisition(1).channel(3).label = 'EyeY';
experiment.hardware.acquisition(1).channel(3).filters = [0,10000];
experiment.hardware.acquisition(1).channel(4).label = 'HandX';
experiment.hardware.acquisition(1).channel(4).filters = [0,10000];
experiment.hardware.acquisition(1).channel(5).label = 'HandY';
experiment.hardware.acquisition(1).channel(5).filters = [0,10000];
experiment.hardware.acquisition(1).channel(6).label = 'Stim';
experiment.hardware.acquisition(1).channel(6).filters = [0,10000];
experiment.hardware.acquisition(1).channel(7).label = 'Raw1';
experiment.hardware.acquisition(1).channel(7).filters = [0,10000];
experiment.hardware.acquisition(1).channel(8).label = 'Raw2';
experiment.hardware.acquisition(1).channel(8).filters = [0,10000];
experiment.hardware.acquisition(1).channel(9).label = 'Raw3';
experiment.hardware.acquisition(1).channel(9).filters = [0,10000];
experiment.hardware.acquisition(1).channel(10).label = 'Raw4';
experiment.hardware.acquisition(1).channel(10).filters = [0,10000];
experiment.hardware.acquisition(1).channel(11).label = 'Raw5';
experiment.hardware.acquisition(1).channel(11).filters = [0,10000];
experiment.hardware.acquisition(1).channel(12).label = 'Raw6';
experiment.hardware.acquisition(1).channel(12).filters = [0,10000];
experiment.hardware.acquisition(1).channel(13).label = 'Raw7';
experiment.hardware.acquisition(1).channel(13).filters = [0,10000];
experiment.hardware.acquisition(1).channel(14).label = 'Raw8';
experiment.hardware.acquisition(1).channel(14).filters = [0,10000];
experiment.hardware.acquisition(1).channel(15).label = 'Mu1';
experiment.hardware.acquisition(1).channel(15).filters = [300,3000];

