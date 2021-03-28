%% Function calcforces
% Computes balance forces (LTT)
% =========================================================================
% Tomas Sinnige - T.Sinnige@tudelft.nl 
% TU Delft - LR - AWEP - Flight Performance and Propulsion
%
% Version: 1.2
% Last updated:  26 February 2021
% First version: 12 October 2017
% =========================================================================
% | Version |    Date   |   Author  |              Changelog              |
% |---------|-----------|-----------|-------------------------------------|
% |   1.2   | 26/02/'21 | T.Sinnige | -) Removed PRS references           |
% |---------|-----------|-----------|-------------------------------------|
% |   1.1   | 17/10/'17 | T.Sinnige | -) Added computation advance ratio  |
% |         |           |           | -) Commented code 'intp' zeroMode   |
% |         |           |           |     -> the sweep option works fine  |
% |---------|-----------|-----------|-------------------------------------|
% |   1.0   | 16/10/'17 | T.Sinnige | -) Added interpolation (in time)    |
% |         |           |           |    between multiple zero-m'ments    |
% |---------|-----------|-----------|-------------------------------------|
% |   0.0   | 12/10/'17 | T.Sinnige | First version                       |
% |---------|-----------|-----------|-------------------------------------|
% =========================================================================
% Inputs:  BAL     - structure containing measurement data balance
%          BAL0    - structure containing zero-measurement data balance
%          idxB    - structure containing indices balance data (W3D)
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
%          testSec - LTT test section number
% -------------------------------------------------------------------------
% Outputs: BAL     - structure containing measurement data balance
%                    (with processed forces and moments added)
%                      BAL.raw       -> raw data (from raw_.. file)
%                      BAL.BAL0      -> zero-measurement data
%                      BAL.B16zeroed -> raw data of B1-6 (steps), with zero
%                                       offset removed
% =========================================================================
function [BAL] = calcforces(BAL,BAL0,D,S,b,c,XmRefB,XmRefM,dAoA,dAoS,modelType,modelPos,testSec)

%% Compute operating conditions
% get operating conditions 
oper = processing.bal.getoper(BAL,testSec);

% write updated data to BAL variable
BAL.q     = oper.qInf; % freestream dynamic pressure [Pa] 
BAL.rho   = oper.rho;  % freestream air density [kg/m^3]
BAL.V     = oper.vInf; % freestream velocity [m/s]
BAL.temp  = oper.tInf; % freestream temperature [K]
BAL.pInf  = oper.pInf; % freestream static pressure [Pa]
BAL.pBar  = oper.pBar; % barometric pressure [Pa]
BAL.nu    = oper.nu;   % freestream kinematic viscosity [m^2/s]

%% Compute advance ratio
BAL.J_M1 = BAL.V./(BAL.rpsM1*D);
BAL.J_M2 = BAL.V./(BAL.rpsM2*D);

%% Compute Reynolds number 
BAL.Re = BAL.V*c./BAL.nu;

%% Get balance calibration factors
[p,pnl,arm,FX_cor,x_bend,y_bend,e] = processing.bal.getcalfactors;

%% Remove zero offset from balance data
BAL = processing.bal.subtractzero(BAL,BAL0);

%% Perform balance calibration
for i=1:size(BAL.Bzeroed,1)
    [F(i,:),M(i,:)] = processing.bal.calibrate(BAL.Bzeroed(i,:),p,pnl,arm,FX_cor,x_bend,y_bend,e);
end

%% Compute nondimensional forces and moments
CF(:,1) = F(:,1) ./ (oper.qInf*S);
CF(:,2) = F(:,2) ./ (oper.qInf*S);
CF(:,3) = F(:,3) ./ (oper.qInf*S);
CM(:,1) = M(:,1) ./ (oper.qInf*S*b);
CM(:,2) = M(:,2) ./ (oper.qInf*S*c);
CM(:,3) = M(:,3) ./ (oper.qInf*S*b);

%% Redefinition of the aerodynamic coefficients in model reference system
% get angles of incidence
AoA = BAL.AoA - dAoA; % correct AoA for constant AoA offset [deg]
AoS = BAL.AoS - dAoS; % correct AoS for constant AoS offset [deg]
BAL.AoA = AoA;
BAL.AoS = AoS;

% Estimate thrust from graph of tc vs ct

% V = 20
idx_V20 = find(abs(BAL.V - 20) <= 3);
BAL.TC1 = zeros(length(BAL.V),1);
BAL.TC2 = zeros(length(BAL.V),1);

% Motor 1
C_T1_V20 = processing.correction.estimatethrust(20,BAL.J_M1(idx_V20));
BAL.TC1(idx_V20) = C_T1_V20.*(BAL.rho(idx_V20).*BAL.rpsM1(idx_V20).^2*D.^4)./(oper.qInf(idx_V20)*S);

% Motor 2
C_T2_V20 = processing.correction.estimatethrust(20,BAL.J_M2(idx_V20));
BAL.TC2(idx_V20) = C_T2_V20.*(BAL.rho(idx_V20).*BAL.rpsM2(idx_V20).^2*D.^4)./(oper.qInf(idx_V20)*S);

% V = 30
idx_V30 = find(abs(BAL.V - 30) <= 3);

% Motor 1
C_T1_V30 = processing.correction.estimatethrust(30,BAL.J_M1(idx_V30));
BAL.TC1(idx_V30) = C_T1_V30.*(BAL.rho(idx_V30).*BAL.rpsM1(idx_V30).^2*D.^4)./(oper.qInf(idx_V30)*S);

% Motor 2
C_T2_V30 = processing.correction.estimatethrust(30,BAL.J_M2(idx_V30));
BAL.TC2(idx_V30) = C_T2_V30.*(BAL.rho(idx_V30).*BAL.rpsM2(idx_V30).^2*D.^4)./(oper.qInf(idx_V30)*S);

% V = 40
idx_V40 = find(abs(BAL.V - 40) <= 3);

% Motor 1
C_T1_V40 = processing.correction.estimatethrust(40,BAL.J_M1(idx_V40));
BAL.TC1(idx_V40) = C_T1_V40.*(BAL.rho(idx_V40).*BAL.rpsM1(idx_V40).^2*D.^4)./(oper.qInf(idx_V40)*S);

% Motor 2
C_T2_V40 = processing.correction.estimatethrust(40,BAL.J_M2(idx_V40));
BAL.TC2(idx_V40) = C_T2_V40.*(BAL.rho(idx_V40).*BAL.rpsM2(idx_V40).^2*D.^4)./(oper.qInf(idx_V40)*S);


% compute forces and moments
if strcmpi(modelType,'aircraft') || strcmpi(modelType,'3dwing')
    % forces in balance axis system
%     CFt = CF(:,1).*cosd(AoA) - CF(:,3).*sind(AoA); % tangential force [N]
    CFt = CF(:,1).*cosd(AoA) - CF(:,3).*sind(AoA) - BAL.TC1 - BAL.TC2; % tangential force without thrust [N]
    CFn = CF(:,3).*cosd(AoA) + CF(:,1).*sind(AoA); % normal force [N]
    CFs = -CF(:,2); % side force [N]
    
    % moments with respect to the moment point in the balance axis system (XmRef)
    CMr = CM(:,1) - XmRefB(2)*CF(:,3) + XmRefB(3)*CF(:,2); % rolling moment [Nm]
    CMp = CM(:,2) - XmRefB(3)*CF(:,1) + XmRefB(1)*CF(:,3); % pitching moment [Nm]
    CMy = CM(:,3) + XmRefB(2)*CF(:,1) - XmRefB(1)*CF(:,2); % yawing moment [Nm] 

    % recalculate moments in the airplane axis system
    CMr = CMr.*cosd(AoA) - CMy.*sind(AoA);
    CMp = CMp;
    CMy = CMy.*cosd(AoA) + CMr.*sind(AoA);
    
    % account for model orientation (upper surface pointing down or up in
    % the wind tunnel)
    if strcmpi(modelPos,'normal')
        CFt = +CF(:,1).*cosd(AoA) + CF(:,3).*sind(AoA) - BAL.TC1 - BAL.TC2; % tangential force without thrust [N]
        CFn = -CF(:,3).*cosd(AoA) + CF(:,1).*sind(AoA); % normal force [N]
        CMr = -CMr;
        CMp = -CMp;
        CMy = -CMy;
    end
    
elseif strcmpi(modelType,'halfwing') % for this case, only the case of angle-of-attack variations is programmed (sideslip always zero)

    % forces
    CFt = CF(:,1) - BAL.TC1 - BAL.TC2; % tangential force without thrust[N]
    CFn = CF(:,2); % normal force [N]
    CFs = CF(:,3); % side force [N]
    
    % moments
    CMr = CM(:,1)          + XmRefB(3)*CF(:,2)*(c/b) - XmRefB(2)*CF(:,3)*(c/b); % rolling moment [Nm]
    CMp = -CM(:,3)*(c/b) - XmRefB(2)*CF(:,1)         + XmRefB(1)*CF(:,2);         % pitching moment [Nm]
    CMy = CM(:,2)*(c/b)  + XmRefB(1)*CF(:,3)*(c/b) - XmRefB(3)*CF(:,1)*(c/b); % yawing moment [Nm] 

    % account for model orientation (upper surface pointing down or up in
    % the wind tunnel)
    if strcmpi(modelPos,'inverted')
        CFn = -CFn; % normal force [N]
        CMr = -CMr;
        CMp = -CMp;
        CMy = -CMy;
    end
   
end % end if statement model orientation

%% Compute lift/drag/pitching moment
CFl   =  (CFn.*cosd(AoA)-CFt.*sind(AoA)); % lift
CFd   =  (CFn.*sind(AoA)+CFt.*cosd(AoA)).*cosd(AoS)+CFs.*sind(AoS); % drag
CFyaw = -(CFn.*sind(AoA)+CFt.*cosd(AoA)).*sind(AoS)+CFs.*cosd(AoS); % sideforce
CMp25c = CMp + CFn*(0.25-XmRefM(1))-CFt*(0.0-XmRefM(3));
% The reference system is such that the model MAC is at (0, 0, 0)
% https://brightspace.tudelft.nl/d2l/le/285568/discussions/threads/79381/View
    
%% Write forces to BAL data structure
BAL.FX  = F(:,1);
BAL.FY  = F(:,2);
BAL.FZ  = F(:,3);
BAL.MX  = M(:,1);
BAL.MY  = M(:,2);
BAL.MZ  = M(:,3);

BAL.CFX = CF(:,1);
BAL.CFY = CF(:,2);
BAL.CFZ = CF(:,3);
BAL.CMX = CM(:,1);
BAL.CMY = CM(:,2);
BAL.CMZ = CM(:,3);


% BAL.N   = Fn;
% BAL.T   = Ft;
% BAL.Y   = Ft;
BAL.CN  = CFn;
BAL.CT  = CFt;
BAL.CY  = CFs;

% BAL.CTp  = CT;
% BAL.TCp  = TC;
% BAL.TCw  = TCw;

% BAL.L   = Fl;
% BAL.D   = Fd;
BAL.CL  = CFl;
BAL.CD  = CFd;
BAL.CYaw  = CFyaw;

% BAL.Mr   = Mr;
% BAL.Mp   = Mp;
% BAL.My   = My;
BAL.CMr  = CMr;
BAL.CMp  = CMp;
BAL.CMp25c  = CMp25c;
BAL.CMy  = CMy;

BAL.b  = b;
BAL.c  = c;
BAL.S  = S;

end % end of function calcforces