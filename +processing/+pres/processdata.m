function [PresData] = processdata(data_dir)

% indices of relevant channels in raw  data files
idxPfus1 = 1; % index of 1st fuselage pressure tap in pressure data
idxPfus2 = 2; % index of 2nd fuselage pressure tap in pressure data
idxPfus3 = 3; % index of 3rd fuselage pressure tap in pressure data
idxPfus4 = 4; % index of 4th fuselage pressure tap in pressure data
idxPfus5 = 5; % index of 5th fuselage pressure tap in pressure data
idxPfus6 = 6; % index of 6th fuselage pressure tap in pressure data
idxPSus  = 19; % index of upstream static pressure measurement in pressure data
idxPTus  = 20; % index of upstream total pressure measurement in pressure data

PresFiles = dir(fullfile(data_dir, '**', '*.cal'));
relativePaths = arrayfun(@util.relativepath,{PresFiles.folder});
[PresFiles.folder] = relativePaths{:};
presFolders = unique({PresFiles.folder});

%% Load raw data
for folder = presFolders % loop over all measurement files
    
    [~,fname,~] = fileparts(folder);
    inFolder = ismember({PresFiles.folder}, folder);
    filename = PresFiles(inFolder).name;
    
    % load raw data
    fid = fopen(fullfile(folder{:}, filename));
    data = textscan(fid,['%f-%f-%f_%f-%f-%f ',repmat('%f ',1,32)]);
    fclose(fid);
   
    % index data and separate into relevant variables
    p = cell2mat(data(7:end));                   % measured pressures [Pa]
    PresData.(fname).Pfus1 = p(:,idxPfus1);   % 1st fus. pressure tap [Pa]
    PresData.(fname).Pfus2 = p(:,idxPfus2);   % 2nd fus. pressure tap [Pa]
    PresData.(fname).Pfus3 = p(:,idxPfus3);   % 3rd fus. pressure tap [Pa]
    PresData.(fname).Pfus4 = p(:,idxPfus4);   % 4th fus. pressure tap [Pa]
    PresData.(fname).Pfus5 = p(:,idxPfus5);   % 5th fus. pressure tap [Pa]
    PresData.(fname).Pfus6 = p(:,idxPfus6);   % 6th fus. pressure tap [Pa]
    PresData.(fname).PSus  = p(:,idxPSus);    % upstream measurement [Pa]
    PresData.(fname).PTus  = p(:,idxPTus);    % upstream measurement [Pa]

    PresData.(fname).year  = data{1};   % year at time of measurement
    PresData.(fname).month = data{2};   % month at time of measurement
    PresData.(fname).day   = data{3};   % day at time of measurement
    PresData.(fname).hour  = data{4};   % hour at time of measurement
    PresData.(fname).min   = data{5};   % min at time of measurement
    PresData.(fname).sec   = data{6};   % sec at time of measurement
     
 end % end for loop over all measurement files