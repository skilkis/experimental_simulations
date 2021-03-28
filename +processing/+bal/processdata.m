%% Function process.m
% Processes LTT balance and pressure data
% =========================================================================
% Tomas Sinnige - T.Sinnige@tudelft.nl 
% TU Delft - LR - AWEP - Flight Performance and Propulsion
%
% Version: 2.0
% Last updated:  17 Jan 2021
% First version: 17 October 2017
% =========================================================================
% | Version |    Date   |   Author  |              Changelog              |
% |---------|-----------|-----------|-------------------------------------|
% |   2.1   | 26/02/'21 |  P.Lopez  |  Removed PRS code                   |
% |---------|-----------|-----------|-------------------------------------|
% |   2.0   | 17/01/'21 |  P.Lopez  |  Change BAL and PRS to structure    |
% |         |           |    (TA)   |  type                               |
% |---------|-----------|-----------|-------------------------------------|
% |   1.0   | 28/05/'19 | T.Sinnige | Made code more generally applicable |
% |         |           |           | by removing part of previous code   |
% |         |           |           | related to isolated propeller setup.|
% |---------|-----------|-----------|-------------------------------------|
% |   0.0   | 17/10/'17 | T.Sinnige | First version                       |
% |---------|-----------|-----------|-------------------------------------|
% =========================================================================
% Inputs:  diskPath - root path on disk where data is stored
%          fnBAL    - filename(s) of the raw BAL data file(s)
%          fn0      - filename(s) of the zero measurement data file(s)
%          idxB     - structure containing indices balance data (W3D)
%          D        - propeller diameter [m]
%          S        - reference area [m^2]
%          b        - wing span [m]
%          c        - mean aerodynamic chord [m]
%          XmRefB   - moment reference points (x,y,z coordinates) in
%                     balance reference system [1/c] 
%          XmRefM   - moment reference points (x,y,z coordinates) in model 
%                     reference system [1/c] 
%          dAoA     - angle of attack offset (subtracted from measured 
%                     values)   
%          dAoS     - angle of sideslip offset (subtracted from measured 
%                     values)
%          modelType- type of wind-tunnel model
%          modelPos - orientation of wind-tunnel model
%          testSec  - LTT test-section number   
% -------------------------------------------------------------------------
% Outputs: BAL      - structure containing all the configurations.
%                     Each configuration contains an structure with the
%                     data of the different parameters.
% =========================================================================
function [BAL] = processdata(diskPath,fnBAL,fn0)

% check whether inputs for raw data and zero data have the same length, and
% if only one zero file is given, apply it to all data sets, otherwise,
% throw an error
if length(fnBAL)~=length(fn0)
    if length(fn0) == 1
        fn0 = repelem(fn0,length(fnBAL));
    else
        error([...
            'Invalid inputs for fnBAL and fn0. Enter 1 zero offset file'...
            ' (fn0) for each raw data file (fnBAL).'...
        ])
    end
end

%% Setup Operating Condition
% get indices balance data files
[idxB, ~] = processing.bal.mapping;

% wing geometry
b     = 1.4*cosd(4); % span [m]
cR    = 0.222; % root chord [m]
cT    = 0.089; % tip chord [m]
S     = b/2*(cT+cR);   % reference area [m^2]
taper = cT/cR; % taper ratio
c     = 2*cR/3*(1+taper+taper^2)/(1+taper); % mean aerodynamic chord [m]

% prop geometry
D        = 0.2032; % propeller diameter [m]
R        = D/2;   % propeller radius [m]

% moment reference points
% Use XMRefB and not XMRefM
% To get the CG to be at 0.55 we can subtract 0.3 on the x coordinate
XmRefB    = [0,0,0.0465/c]; % moment reference points (x,y,z coordinates) in balance reference system [1/c] 
XmRefM    = [0.25,0,0];     % moment reference points (x,y,z coordinates) in model reference system [1/c] 

% incidence angle settings
dAoA      = 0.0; % angle of attack offset (subtracted from measured values)   
dAoS      = 0.0;  % angle of sideslip offset (subtracted from measured values)
modelType = 'aircraft'; % options: aircraft, 3dwing, halfwing
modelPos  = 'inverted'; % options: normal, inverted
testSec   = 5;    % test-section number

%% Process data
for i=1:length(fnBAL) % loop over the files that are to be loaded
    
    % extract configuration name
    BAL.config(i) = extractBetween(char(fnBAL{i}),'raw_','.txt');
    
    if isletter(BAL.config{i}(1))
        % first character is letter so can be used as field name -> continue
    else
        if strcmpi(BAL.config{i}(1),'+')
            BAL.config{i}(1) = 'p';
        elseif strcmpi(BAL.config{i}(1),'-')
            BAL.config{i}(1) = 'm';
        else
           error('Unexpected character used as first character of filename. Please add if statement here that catches this exception and changes the character into a letter.') 
        end
    end
    
    % give status update
    disp([...
        'Processing balance data; configuration ''',...
        BAL.config{i},...
        '''; filename ''',...
        fnBAL{i},...
        '''.'...
    ]);
    % load zero-measurement data (BAL0)
    BAL.windOff.(BAL.config{i}) = processing.bal.readzerodata(...
        diskPath,...
        fn0{i},...
        idxB...
    );
    % load measurement data (BAL)
    BAL.windOn.(BAL.config{i}) = processing.bal.readdata(...
        diskPath,...
        fnBAL{i},...
        idxB...
    );
    
    % compute calibrated balance forces and moments
    BAL.windOn.(BAL.config{i}) = processing.bal.calcforces(BAL.windOn.(BAL.config{i}),BAL.windOff.(BAL.config{i}),D,S,b,c,XmRefB,XmRefM,dAoA,dAoS,modelType,modelPos,testSec);    
    
end % end for loop over filenames

end % end of function process