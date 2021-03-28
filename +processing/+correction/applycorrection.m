function [BalDataCorr] = applycorrection(BAL)
    % Outputs: BAL_new - updated measurement file

    % Inputs: BAL - To be corrected measurement file, idxB - Master index file
    % from SUP_getIdx, tailoff - 1 if this is tailoff data, powered - 1 if
    % this is powered data

    % Constants
    s_wing = 0.6985;
    c_wing = 0.165;
    t_wing = 0.1583*c_wing;
    V_wing = 0.0030229;
    V_fuselage = 0.0160632 + 0.0004491 + 0.0035296;
    V_tail = 0.0024485 + 0.0003546;
    K1_wing = V_wing / (2*s_wing*c_wing*t_wing);
    K1_tail = 1.03;
    K3 = 0.91;
    Tau_fuselage = 0.8;
    Tau_wing = 0.88;
    Tau_tail = 0.86;
    A_tunnel = 2.044;
    A_VTail = 0.0415;
    V_total = V_wing + V_fuselage + V_tail;  % Volume not velocity

    for field = fieldnames(BAL.windOn)'
        % Making copies of measured data to mutate differently
        % Base = Base data without any corrections
        % SB = Solid Blockage
        % WB = Wake Blockage
        % SS = Slipstream Blockage
        % SL = Streamline Correction due to lift interference
        % DW = Downwash Correciton due to vertical tail
        % Total = All corrections applied
        Base = BAL.windOn.(field{:});
        SB = BAL.windOn.(field{:});
        WB = BAL.windOn.(field{:});
        SS = BAL.windOn.(field{:});
        SL = BAL.windOn.(field{:});
        DW = BAL.windOn.(field{:});
        Total = BAL.windOn.(field{:});

        % Blockage

        % Solid Blockage

        % Compute factors based on constants provided
        Epsilon_wing = ones(size(Base.CD)) * (K1_wing*Tau_wing*V_wing)/(A_tunnel)^(3/2);
        Epsilon_fuselage = ones(size(Base.CD)) * (K3*Tau_fuselage*V_fuselage)/(A_tunnel)^(3/2);
        Epsilon_tail = ones(size(Base.CD)) * (K1_tail*Tau_tail*V_tail)/(A_tunnel)^(3/2);

        % Wake Blockage
        S = 0.2172;
        Cd_u = Base.CD;
        Cd_0 = 0.0583;          % pre-recorded for a sideslip variation and hard-coded
        Cl_u_2 = Base.CL.^2;

        % Slope of Cd vs Cl^2 graph used to estimate induced drag - (for this data set, there is no variation of V or RPM in one polar)
        Slope = mean(diff(Cd_u)) ./ mean(diff(Cl_u_2));

        % Compute Cd_0 and induced drag, check if separated flow
        Cl_u = Base.CL;
        Cd_i = (Cl_u.^2 * Slope) - Cd_0;
        Cd_s = Cd_u - Cd_i - Cd_0;

        % Check if separated and assign non-negative value
        if Cd_s < 0
            Cd_s = 0;
        end

        Epsilon_Wake_Blockage = (S*Cd_0 + 5*S*Cd_s)/(4*A_tunnel);

        % Slipstream blockage
        Breadth_tunnel = 1.8;
        Height_tunnel = 1.24;
        Sp = 0.0314;
        
        % Perform slipstream blockage correction with thrust coefficient
        Thrust_coefficient = Base.TC1 + Base.TC2;
        Epsilon_Slipstream = -Sp*(sqrt(1+(8*Thrust_coefficient)/pi())-1)/(2*pi()*Breadth_tunnel*Height_tunnel);
        Epsilon_total = Epsilon_wing + Epsilon_fuselage + Epsilon_tail + Epsilon_Wake_Blockage + Epsilon_Slipstream;

        SB.V = Base.V.*(1 + Epsilon_wing + Epsilon_fuselage + Epsilon_tail);
        WB.V = Base.V.*(1+Epsilon_Wake_Blockage);
        SS.V = Base.V.*(1+Epsilon_Slipstream);
        Total.V = Base.V.*(1 + Epsilon_total);
        
        % Correct Solid Blockages
        SB.CD = SB.CD.*Base.V.^2 ./ SB.V.^2;
        SB.CL = SB.CL.*Base.V.^2 ./ SB.V.^2;
        SB.CY = SB.CY.*Base.V.^2 ./ SB.V.^2;
        SB.CMy = SB.CMy.*Base.V.^2 ./ SB.V.^2;
        
        % Correct Only Wake Blockage
        WB.CD = WB.CD.*Base.V.^2 ./ WB.V.^2;
        WB.CL = WB.CL.*Base.V.^2 ./ WB.V.^2;
        WB.CY = WB.CY.*Base.V.^2 ./ WB.V.^2;
        WB.CMy = WB.CMy.*Base.V.^2 ./ WB.V.^2;
        
        % Correct Only Slipstream Blockage
        SS.CD = SS.CD.*Base.V.^2 ./ SS.V.^2;
        SS.CL = SS.CL.*Base.V.^2 ./ SS.V.^2;
        SS.CY = SS.CY.*Base.V.^2 ./ SS.V.^2;
        SS.CMy = SS.CMy.*Base.V.^2 ./ SS.V.^2;
        
        % Correct coefficents due to all blockages
        Total.CD = Total.CD.*Base.V.^2 ./ Total.V.^2;
        Total.CL = Total.CL.*Base.V.^2 ./ Total.V.^2;
        Total.CY = Total.CY.*Base.V.^2 ./ Total.V.^2;
        Total.CMy = Total.CMy.*Base.V.^2 ./ Total.V.^2;

        CN_Tail = Total.CY;

        % Streamline Curvature correction (lift interference) for the
        % vertical tail due to its own circulation
        Tau_2 = 0.07;           % Taper = 0.68, MAC = 0.17, Lt = 0.085, Lt/2H = 0.11, Fig. 10.37 in the Barlow book
        Delta = 0.12;          % Elliptic loading for a closed elliptic jet, Fig. 10.32 in the Barlow book
        diff_AOS = mean(diff(Total.AoS));
        if diff_AOS == 0
            diff_AOS = 0.3778;
        end
        Slope_CN_vs_Beta = mean(diff(CN_Tail)) ./ diff_AOS;        Delta_Beta_Total_SC = (1+Tau_2)*Delta*(A_VTail/A_tunnel)*Total.CY;
        Delta_Cm = 0.125*Delta_Beta_Total_SC*Slope_CN_vs_Beta;               % moment = arm*db*dCn/db (Barlow page 377)
        Total.CD = Total.CD + Delta*(A_VTail/A_tunnel)*Total.CL.^2;           % Equation 2.10 pre-test report

        Total.CMy = Total.CMy + Delta_Cm;
        
        % Downwash correction
        Beta_u = Total.AoS;
        Beta_downwash = Beta_u+(Delta*(A_VTail/A_tunnel)*CN_Tail*57.3);

        % Correct beta for both
        Total.AoS = Delta_Beta_Total_SC + Beta_downwash;
        
        % Only slipstream correction
        CN_Tail = SL.CY;
        diff_AOS = mean(diff(SL.AoS));
        if diff_AOS == 0
            diff_AOS = 0.3778;
        end
        Slope_CN_vs_Beta = mean(diff(CN_Tail)) ./ diff_AOS;
        Delta_Beta_Total_SC = (1+Tau_2)*Delta*(A_VTail/A_tunnel)*SL.CY;
        Delta_Cm = 0.125*Delta_Beta_Total_SC*Slope_CN_vs_Beta;               % moment = arm*db*dCn/db (Barlow page 377)
        SL.CD = SL.CD + Delta*(A_VTail/A_tunnel)*SL.CL.^2;           % Equation 2.10 pre-test report

        SL.CMy = SL.CMy + Delta_Cm;
        SL.AoS = SL.AoS + Delta_Beta_Total_SC;
        
        % Downwash correction
        Beta_u = DW.AoS;
        Beta_downwash = Beta_u+(Delta*(A_VTail/A_tunnel)*CN_Tail*57.3);
        DW.AoS = Beta_downwash;
        
%         % Buoyancy correction, Equation 2.2 pre-test report
%         dp_by_dv = ;                        % from pressure taps along the fuselage
%         del_D = dp_by_dv * V_total;
%         BAL_new.CD = BAL_new.CD - del_D/(0.5*BAL_new.rho.*BAL_new.V.^2.*BAL_new.S);
        BalDataCorr.(field{:}).SB = SB;
        BalDataCorr.(field{:}).WB = WB;
        BalDataCorr.(field{:}).SS = SS;
        BalDataCorr.(field{:}).SL = SL;
        BalDataCorr.(field{:}).DW = DW;
        BalDataCorr.(field{:}).Total = Total;
    end
end
