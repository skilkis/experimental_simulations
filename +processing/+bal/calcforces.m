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
function [BAL] = calcforces(BAL,BAL0,D,S,b,c,XmRefB,XmRefM,YmProp,dAoA,dAoS,modelType,modelPos,testSec)

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

%% Estimate thrust from graph of tc vs ct

% Iterating through velocities in nearest tenths (i.e. 20, 30, 40)
BAL.TC1 = zeros(length(BAL.V),1);
BAL.TC2 = zeros(length(BAL.V),1);
BAL.FT1 = zeros(length(BAL.V),1);
BAL.FT2 = zeros(length(BAL.V),1);
BAL.CTM1 = zeros(length(BAL.V),1);
BAL.CTM2 = zeros(length(BAL.V),1);
Vmeas = unique(round(BAL.V, -1))'; %Measured velocities (i.e. 20, 30, 40)

% Estimating thrust coefficient and normalized thrust force per motor
% TODO check if the correct diameter was used
for V = Vmeas
    idx = find(abs(BAL.V - V) <= 3);
    BAL.TC1(idx) = processing.correction.estimatethrust(V,BAL.J_M1(idx));
    BAL.TC2(idx) = processing.correction.estimatethrust(V,BAL.J_M2(idx));
    % Thrust Force due to motor 1 and 2:
    BAL.FT1(idx) = BAL.TC1(idx).*(BAL.rho(idx).*BAL.rpsM1(idx).^2*D.^4);
    BAL.FT2(idx) = BAL.TC2(idx).*(BAL.rho(idx).*BAL.rpsM2(idx).^2*D.^4);
    % Normalized tangential force (thrust) due to motor 1 and 2:
    BAL.CTM1(idx) = BAL.FT1(idx)./(oper.qInf(idx)*S);
    BAL.CTM2(idx) = BAL.FT2(idx)./(oper.qInf(idx)*S);
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

% compute forces and moments
% if strcmpi(modelType,'aircraft') || strcmpi(modelType,'3dwing')
    % forces in balance axis system
%     CFt = CF(:,1).*cosd(AoA) - CF(:,3).*sind(AoA); % tangential force [N]
CFt = CF(:,1).*cosd(AoA) - CF(:,3).*sind(AoA) - BAL.CTM1 - BAL.CTM2; % tangential force without thrust [N]
CFn = CF(:,3).*cosd(AoA) + CF(:,1).*sind(AoA); % normal force [N]
CFs = -CF(:,2); % side force [N]

% moments with respect to the moment point in the balance axis system (XmRef)
CMr = CM(:,1) - XmRefB(2)*CF(:,3) + XmRefB(3)*CF(:,2); % rolling moment [Nm]
CMp = CM(:,2) - XmRefB(3)*CF(:,1) + XmRefB(1)*CF(:,3); % pitching moment [Nm]
CMy = CM(:,3) + XmRefB(2)*CF(:,1) - XmRefB(1)*CF(:,2); % yawing moment [Nm]
% Isolating yawing moment due to aerodynamic effects and not thrust
% Motor 1 creates a positive yawing moment, while motor 2 creates
% a negative one. Therefore, the yawing moment due to thrust is
% added to compensate for M2, while it is subtracted for M1.
CMya = CMy + (BAL.CTM2 - BAL.CTM1) * YmProp;  % Aerodynamic yawing moment

% recalculate moments in the airplane axis system
CMr = CMr.*cosd(AoA) - CMy.*sind(AoA);
CMp = CMp;
CMy = CMy.*cosd(AoA) + CMr.*sind(AoA);
CMya = CMya.*cosd(AoA) + CMr.*sind(AoA);
    
    % account for model orientation (upper surface pointing down or up in
    % the wind tunnel)
%     if strcmpi(modelPos,'normal')
%         CFt = +CF(:,1).*cosd(AoA) + CF(:,3).*sind(AoA) - BAL.CFtM1 - BAL.CFtM2; % tangential force without thrust [N]
%         CFn = -CF(:,3).*cosd(AoA) + CF(:,1).*sind(AoA); % normal force [N]
%         CMr = -CMr;
%         CMp = -CMp;
%         CMy = -CMy;
%         CMya = -CMya;
%     end
    
% elseif strcmpi(modelType,'halfwing') % for this case, only the case of angle-of-attack variations is programmed (sideslip always zero)
% 
%     % forces
%     CFt = CF(:,1) - BAL.CFtM1 - BAL.CFtM2; % tangential force without thrust[N]
%     CFn = CF(:,2); % normal force [N]
%     CFs = CF(:,3); % side force [N]
%     
%     % moments
%     CMr = CM(:,1)          + XmRefB(3)*CF(:,2)*(c/b) - XmRefB(2)*CF(:,3)*(c/b); % rolling moment [Nm]
%     CMp = -CM(:,3)*(c/b) - XmRefB(2)*CF(:,1)         + XmRefB(1)*CF(:,2);         % pitching moment [Nm]
%     CMy = CM(:,2)*(c/b)  + XmRefB(1)*CF(:,3)*(c/b) - XmRefB(3)*CF(:,1)*(c/b); % yawing moment [Nm] 
% 
%     % account for model orientation (upper surface pointing down or up in
%     % the wind tunnel)
%     if strcmpi(modelPos,'inverted')
%         CFn = -CFn; % normal force [N]
%         CMr = -CMr;
%         CMp = -CMp;
%         CMy = -CMy;
%     end
%    
% end % end if statement model orientation

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
BAL.CMya = CMya;

BAL.b  = b;
BAL.c  = c;
BAL.S  = S;

end % end of function calcforces