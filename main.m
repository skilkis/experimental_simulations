clc
clear all
close all

% MicData = processing.mic.processdata(fullfile('data', 'MICS'));

raw_balance_files = {...
    'raw_rudder0.txt',...
    'raw_rudder10.txt',...
    'raw_j_sweep.txt',...
};
BalData = processing.bal.processdata(...
    fullfile('.', 'data', 'BAL'),...
    raw_balance_files,...
    {'zer_ 20210303_edited.txt'}...
);
% PresData = processing.pres.processdata(fullfile('data', 'PRES'));

% Applying Corrections
[BalDataCorr,BalDataCorr_Indiv] = processing.correction.applycorrection(BalData);
% % Saving results to data/processed.mat
% cache_file = sprintf('%s.mat', fullfile('data', 'processed'));
% [cache_path, ~, ~] = fileparts(cache_file);
% if ~exist(cache_path, 'dir')
%     mkdir(cache_path)
% end
% save(cache_file,'BalData', 'BalDataCorr', 'MicData', 'PresData');