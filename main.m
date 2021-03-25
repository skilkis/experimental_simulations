[PXX,SPL,f,opp] = processing.mic.processdata('.\data\MICS\polar1\');

raw_balance_files = {...
    'raw_rudder0.txt',...
    'raw_rudder10.txt',...
    'raw_j_sweep.txt',...
};
BAL = processing.bal.processdata(...
    fullfile(".", 'data', 'BAL'),...
    raw_balance_files,...
    {'zer_ 20210303_edited.txt'}...
);
BAL_corr = processing.correction.applycorrection(BAL);