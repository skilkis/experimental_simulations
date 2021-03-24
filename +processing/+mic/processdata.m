function [PXX,SPL,f,opp] = processdata(data_dir)

cache_file = fullfile(...
    '.cache',...
    'processdata',...
    sprintf("%s.mat",util.hashstring(data_dir))...
);
if exist(cache_file, 'file')
    % Retrieve persisted values if cache exists
    cache_results = load(cache_file);
    PXX = cache_results.PXX;
    SPL = cache_results.SPL;
    f = cache_results.f;
    opp = cache_results.opp;
else

% Disable warnings due to rounding
warning('off', 'signal:check_order:InvalidOrderRounding')
% freestream conditions ONLY FOR GROUPS 1-13 ! 
% add 1 entry per measurement point ! 
operManual.vInf = [40,40,40];

% settings for spectral analysis (pwelch function)
nB       = 94;      % number of bins
nOverlap = 0;       % number of samples for overlap
fs       = 51200;   % sampling frequency [Hz]

% propeller geometry
D = 0.2032; % propeller diameter [m]

%% Load microphone calibration data
[parent_dir, ~, ~] = fileparts(mfilename("fullpath"));
micCal = load(fullfile(parent_dir, "calibrationData.mat"));

%% Loop over all TDMS files of name "Measurement_i.tdms" in data_dir
done = 0;
idx =0;
while done == 0

    % Update counter 
    idx = idx +1;
        
    % Construct filenames (mic data and operating conditions)
    fn = ['Measurement_' num2str(idx) '.tdms'];
    TDMSpath = fullfile(data_dir, fn);
    OPPpath = [TDMSpath(1:end-4),'txt'];
    
    % Check if this file exists
    if exist(TDMSpath,'file')
        
        % Load data (dataOut should have just one channel group) keep an
        % eye on the amount of RAM used, it may fill up
        dataOut{idx} = processing.mic.readtdms(TDMSpath);
        dataOpp(:,idx) = load(OPPpath);        
        disp(['Loaded file ' num2str(idx)])
        
       % Apply calibration
        for i=1:6 % loop over the microphones     
            [pMic{idx}(i,:),~,~] = processing.mic.applycalcurve(...
                fs,...
                dataOut{idx}{1}(:,i)-mean(dataOut{idx}{1}(:,i)),...
                micCal.f_oct,...
                micCal.F_mics(i,:)...
            );
        end     
        
        % perform spectral analysis
        w        = hann(length(pMic{idx})/nB); % window
        wOvrlp   = length(w)*nOverlap; % overlap window
        [PXX{idx},SPL{idx},f{idx}] = processing.mic.spectralanalysis(...
            pMic{idx}.',...
            w,...
            wOvrlp,...
            fs,...
            []...
        );
            
    else % If no file has been found for this index, then break
        done = 1;
        
        % extract some data from the operating-condition file
        if size(dataOpp,1)==33
            opp.vInf   = dataOpp(6,:);   % freestream velocity [m/s]
            opp.AoA    = dataOpp(12,:);  % angle of attack [deg]
            opp.AoS    = dataOpp(13,:);  % angle of sideslip [deg]
            opp.RPS_M1 = dataOpp(14,:);  % RPS motor 1 [Hz] 
            opp.RPS_M2 = dataOpp(21,:);  % RPS motor 2 [Hz] 
        else
            % no tunnel data stored --> assume input freestream velocity
            % values (provide as input at top of file)
            opp.vInf   = operManual.vInf;   % freestream velocity [m/s]
            opp.AoA    = dataOpp(1,:);  % angle of attack [deg]
            opp.AoS    = dataOpp(2,:);  % angle of sideslip [deg]
            opp.RPS_M1 = dataOpp(3,:); % RPS motor 1 [Hz]  
            opp.RPS_M2 = dataOpp(10,:); % RPS motor 2 [Hz] 
        end
        opp.J_M1   = opp.vInf./(opp.RPS_M1*D); % advance ratio motor 1
        opp.J_M2   = opp.vInf./(opp.RPS_M2*D); % advance ratio motor 2
        
    end % end if statement check whether file exists
end % end while loop over files

% Save processed results to disk
[cache_path, ~, ~] = fileparts(cache_file);
if ~exist(cache_path, 'dir')
    mkdir(cache_path)
end
save(cache_file,'PXX','SPL','f','opp')
end
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