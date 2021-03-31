function [MicData] = processdata(data_dir, varargin)
%% Process Optional Arguments
arg = inputParser; % Analyzes passed arguments
addOptional(arg, 'remove_zero', false, @islogical);
parse(arg, varargin{:});
remove_zero = arg.Results.remove_zero;

%% Define Constants
% freestream conditions ONLY FOR GROUPS 1-13 ! 
% add 1 entry per measurement point ! 
operManual.vInf = [40,40,40];

% settings for spectral analysis (pwelch function)
nB       = 94;      % number of bins
nOverlap = 0;       % number of samples for overlap
fs       = 51200;   % sampling frequency [Hz]

% propeller geometry
D = 0.2032; % propeller diameter [m]

% Disable warnings due to rounding
warning('off', 'signal:check_order:InvalidOrderRounding')

%% Load microphone calibration data
[parent_dir, ~, ~] = fileparts(mfilename("fullpath"));
micCal = load(fullfile(parent_dir, "calibrationData.mat"));

%% Parsing inputs and checking directory given by data_dir for .tdms files
% Parse input to obtain all TDMS files in the data_dir using glob pattern
if exist(data_dir, 'dir')
    TDMSFiles = dir(fullfile(data_dir, '**', '*.tdms'));
    relativePaths = arrayfun(@util.relativepath,{TDMSFiles.folder});
    [TDMSFiles.folder] = relativePaths{:};
    tdmsFolders = unique({TDMSFiles.folder});
else
    if ~isempty(getenv('CI'))
        % Handle GitHub Actions not having access to MIC files
        TDMSFiles = NaN;
        tdmsFolders = {...
            fullfile('.', 'data', 'MICS', 'rudder0'),...
            fullfile('.', 'data', 'MICS', 'rudder10'),...
            fullfile('.', 'data', 'MICS', 'j_sweep'),...
        };
    else
        error('%s is not a valid directory', data_dir)
    end
end

%% Iterating through all unique folders with .tdms files
for folder = tdmsFolders
    [~, fname, ~] = fileparts(folder);
    cache_file = fullfile('.cache','MICS',sprintf("%s.mat", fname));
    if ~exist(cache_file, 'file')
        % Indices of files that match the current folder
        inFolder = ismember({TDMSFiles.folder}, folder);
        % Run processing on cache-miss and save processed results to disk
        [PXX,SPL,f,opp] = processtdmsfiles(folder{:}, TDMSFiles(inFolder));
        [cache_path, ~, ~] = fileparts(cache_file);
        if ~exist(cache_path, 'dir')
            mkdir(cache_path)
        end
        save(cache_file,'PXX','SPL','f','opp')
        fprintf('Saved processed results to %s\n', cache_file)
    else
        % Retrieve persisted values if cache exists
        cache_results = load(cache_file);
        PXX = cache_results.PXX;
        SPL = cache_results.SPL;
        f = cache_results.f;
        opp = cache_results.opp;
    end
    % Filter zero measurements
    if remove_zero
        % For some reason tunnel measurements at V=0 show up as ~1.3 m/s
        non_zero_idx = opp.vInf > 2;
        % Hardcoding removal of duplicate data point
        if strcmp(fname, 'rudder0')
            non_zero_idx(9) = false;
        end
        PXX = PXX(:, non_zero_idx);
        SPL = SPL(:, non_zero_idx);
        f = f(:, non_zero_idx);
        for field = fieldnames(opp)'
            opp.(field{:}) = opp.(field{:})(non_zero_idx);
        end
    end
    % Assign output to MicData struct by foldername
    MicData.(fname).PXX = PXX;
    MicData.(fname).SPL = SPL;
    MicData.(fname).f = f;
    MicData.(fname).opp = opp;
end

%% Loop over all TDMS files of name "Measurement_i.tdms" in data_dir
    function [PXX,SPL,f,opp] = processtdmsfiles(folder, files)
        fprintf('Processing TDMS Files in %s\n', folder)
        nFiles = length(files);
        for i = 1:nFiles
            % Construct filenames (mic data and operating conditions)
            fn = ['Measurement_' num2str(i) '.tdms'];
            TDMSpath = fullfile(folder, fn);
            OPPpath = [TDMSpath(1:end-4),'txt'];

            % Check if this file exists
            if exist(TDMSpath,'file')

                % Load data (dataOut should have just one channel group)
                % keep aneye on the amount of RAM used, it may fill up
                dataOut{i} = processing.mic.readtdms(TDMSpath);
                dataOpp(:,i) = load(OPPpath);
                
                fprintf('\t...Processing file (%d/%d)\n', i, nFiles)
               % Apply calibration
                for j=1:6 % loop over the microphones     
                    [pMic{i}(j,:),~,~] = processing.mic.applycalcurve(...
                        fs,...
                        dataOut{i}{1}(:,j)-mean(dataOut{i}{1}(:,j)),...
                        micCal.f_oct,...
                        micCal.F_mics(j,:)...
                    );
                end     

                % perform spectral analysis
                w = hann(length(pMic{i})/nB); % window
                wOvrlp = length(w)*nOverlap; % overlap window
                [PXX{i},SPL{i},f{i}] = processing.mic.spectralanalysis(...
                    pMic{i}.',...
                    w,...
                    wOvrlp,...
                    fs,...
                    []...
                );
            end
        end
        % extract some data from the operating-condition file
        if size(dataOpp,1)==33
            opp.vInf = dataOpp(6,:);    % freestream velocity [m/s]
            opp.AoA = dataOpp(12,:);    % angle of attack [deg]
            opp.AoS = dataOpp(13,:);    % angle of sideslip [deg]
            opp.RPS_M1 = dataOpp(14,:); % RPS motor 1 [Hz]
            opp.RPS_M2 = dataOpp(21,:); % RPS motor 2 [Hz]
        else
            % no tunnel data stored --> assume input freestream velocity
            % values (provide as input at top of file)
            opp.vInf = operManual.vInf; % freestream velocity [m/s]
            opp.AoA = dataOpp(1,:);     % angle of attack [deg]
            opp.AoS = dataOpp(2,:);     % angle of sideslip [deg]
            opp.RPS_M1 = dataOpp(3,:);  % RPS motor 1 [Hz]
            opp.RPS_M2 = dataOpp(10,:); % RPS motor 2 [Hz]
        end
        opp.J_M1 = opp.vInf./(opp.RPS_M1*D); % Advance ratio M1
        opp.J_M2 = opp.vInf./(opp.RPS_M2*D); % Advance ratio M2

        if isempty(getenv('CI'))
            figure('Name','Spectra')
            for i=1:6
                subplot(2,3,i), box on, hold on
                plot(f{1},SPL{1}(:,i),'b')
                plot(f{2},SPL{2}(:,i),'r')
                xlim([0 2000])
                xlabel('Frequency f [Hz]')
                ylabel('SPL [dB]')
                title(['Mic ',num2str(i)])
                ylim([50 120])
            end
        end
    end
end