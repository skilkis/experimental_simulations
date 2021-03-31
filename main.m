clc
clear all
close all

%% Setup Experiment Geometry
geom.b = 1.4*cosd(4); % span [m]
geom.cR = 0.222; % root chord [m]
geom.cT = 0.089; % tip chord [m]
geom.S = geom.b/2*(geom.cT+geom.cR); % reference area [m^2]
geom.taper = geom.cT/geom.cR; % taper ratio
geom.c = 2*geom.cR/3*(1+geom.taper+geom.taper^2)/(1+geom.taper); % mean aerodynamic chord [m]
geom.D = 0.2032; % propeller diameter [m]
geom.R = geom.D/2; % propeller radius [m]
geom.XmRefB = [0,0,0.0465/geom.c]; % moment reference points (x,y,z coordinates) in balance reference system [1/c] 
geom.testSec = 5; % test-section number
geom.nBlades = 6;  % Number of propeller blades

MicData = processing.mic.processdata(...
    fullfile('data', 'MICS'),...
    'remove_zero', true...
);

raw_balance_files = {...
    'raw_rudder0.txt',...
    'raw_rudder10.txt',...
    'raw_j_sweep.txt',...
};
BalData = processing.bal.processdata(...
    fullfile('.', 'data', 'BAL'),...
    raw_balance_files,...
    {'zer_ 20210303_edited.txt'},...
    geom...
);
PresData = processing.pres.processdata(fullfile('data', 'PRES'));

% Normalizing Sound Pressure
MicData = processing.mic.normalize(MicData, BalData, geom.D, geom.nBlades);

% Applying Corrections
BalDataCorr = processing.correction.applycorrection(BalData);

% Saving results to data/processed.mat
cache_file = sprintf('%s.mat', fullfile('data', 'processed'));
[cache_path, ~, ~] = fileparts(cache_file);
if ~exist(cache_path, 'dir')
    mkdir(cache_path)
end
save(cache_file,'BalData', 'BalDataCorr', 'MicData', 'PresData');