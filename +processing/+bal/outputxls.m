%% Functionxls_file_output.m
% Write Balance data in an excel file
% =========================================================================
% Tomas Sinnige - T.Sinnige@tudelft.nl 
% TU Delft - LR - AWEP - Flight Performance and Propulsion
%
% Version: 0
% Last updated:  26 January 2021
% First version: 26 January 2021
% =========================================================================
% | Version |    Date   |   Author  |              Changelog              |
% |---------|-----------|-----------|-------------------------------------|
% |   0.0   | 26/01/'21 |Pedro Lopez| First version                       |
% |         |           |    (TA)   |                                     |
% |---------|-----------|-----------|-------------------------------------|
% =========================================================================
% Inputs:  BAL      - structure containing measurement data balance
%          fnBAL    - filename(s) of the raw BAL data file(s)
%          idxB     - structure containing indices balance data (W3D)
% -------------------------------------------------------------------------
function outputxls(BAL,fnBAL)
for i=1:length(fnBAL)
BAL.windOn.(BAL.config{i}) = rmfield( BAL.windOn.(BAL.config{i}) , 'B' );
BAL.windOn.(BAL.config{i}) = rmfield( BAL.windOn.(BAL.config{i}) , 'B16zeroed' );
BAL.windOn.(BAL.config{i}) = rmfield( BAL.windOn.(BAL.config{i}) , 'b' );
BAL.windOn.(BAL.config{i}) = rmfield( BAL.windOn.(BAL.config{i}) , 'c' );
BAL.windOn.(BAL.config{i}) = rmfield( BAL.windOn.(BAL.config{i}) , 'S' );
% writetable(struct2table(BAL.windOn.(BAL.config{i}),'AsArray',0),'OUTPUT.xls','AutoFitWidth',true,'Sheet', BAL.config{i});
writetable(struct2table(BAL.windOn.(BAL.config{i}),'AsArray',0),'OUTPUT.xls','Sheet', BAL.config{i});

end