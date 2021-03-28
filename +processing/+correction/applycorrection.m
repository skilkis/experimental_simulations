function [BAL_new,BAL_new_individual_return] = applycorrection(BAL)
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
    V_total = V_wing + V_fuselage + V_tail;

    for field = fieldnames(BAL.windOn)'
        % TODO find more semantic name for BAL_new
        BAL_new = BAL.windOn.(field{:});
        BAL_new_individual(1:5) = BAL_new;

        % Blockage

        % Solid Blockage

        % Compute factors based on constants provided
        Epsilon_wing = ones(size(BAL_new.CD)) * (K1_wing*Tau_wing*V_wing)/(A_tunnel)^(3/2);
        Epsilon_fuselage = ones(size(BAL_new.CD)) * (K3*Tau_fuselage*V_fuselage)/(A_tunnel)^(3/2);
        Epsilon_tail = ones(size(BAL_new.CD)) * (K1_tail*Tau_tail*V_tail)/(A_tunnel)^(3/2);

        % Wake Blockage
        S = 0.2172;
        Cd_u = BAL_new.CD;
        Cd_0 = 0.0583;          % pre-recorded for a sideslip variation and hard-coded
        Cl_u_2 = BAL_new.CL.^2;

        % Slope of Cd vs Cl^2 graph used to estimate induced drag - (for this data set, there is no variation of V or RPM in one polar)
        Slope = mean(diff(Cd_u)) ./ mean(diff(Cl_u_2));

        % Compute Cd_0 and induced drag, check if separated flow
        Cl_u = BAL_new.CL;
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

        V = BAL_new.V;

        % Perform slipstream blockage correction with thrust coefficient
        % value
        Thrust_coefficient = BAL_new.TC1 + BAL_new.TC2;
        Epsilon_Slipstream = -Sp*(sqrt(1+(8*Thrust_coefficient)/pi())-1)/(2*pi()*Breadth_tunnel*Height_tunnel);

        Epsilon_total = Epsilon_wing + Epsilon_fuselage + Epsilon_tail + Epsilon_Wake_Blockage + Epsilon_Slipstream;

        % Update velocity due to blokages
        BAL_new.V = V.*(1 + Epsilon_total);
        BAL_new_individual(1).V = V.*(1 + Epsilon_wing + Epsilon_fuselage + Epsilon_tail);
        BAL_new_individual(2).V = V.*(1+Epsilon_Wake_Blockage);
        BAL_new_individual(3).V = V.*(1+Epsilon_Slipstream);

        % Correct important coefficents due to blockages
        BAL_new.CD = BAL_new.CD.*V.^2 ./ BAL_new.V.^2;
        BAL_new.CL = BAL_new.CL.*V.^2 ./ BAL_new.V.^2;
        BAL_new.CY = BAL_new.CY.*V.^2 ./ BAL_new.V.^2;
        BAL_new.CMy = BAL_new.CMy.*V.^2 ./ BAL_new.V.^2;
        
        BAL_new_individual(1).CD = BAL_new_individual(1).CD.*V.^2 ./ BAL_new_individual(1).V.^2;
        BAL_new_individual(1).CL = BAL_new_individual(1).CL.*V.^2 ./ BAL_new_individual(1).V.^2;
        BAL_new_individual(1).CY = BAL_new_individual(1).CY.*V.^2 ./ BAL_new_individual(1).V.^2;
        BAL_new_individual(1).CMy = BAL_new_individual(1).CMy.*V.^2 ./ BAL_new_individual(1).V.^2;
        
        BAL_new_individual(2).CD = BAL_new_individual(2).CD.*V.^2 ./ BAL_new_individual(2).V.^2;
        BAL_new_individual(2).CL = BAL_new_individual(2).CL.*V.^2 ./ BAL_new_individual(2).V.^2;
        BAL_new_individual(2).CY = BAL_new_individual(2).CY.*V.^2 ./ BAL_new_individual(2).V.^2;
        BAL_new_individual(2).CMy = BAL_new_individual(2).CMy.*V.^2 ./ BAL_new_individual(2).V.^2;
        
        BAL_new_individual(3).CD = BAL_new_individual(3).CD.*V.^2 ./ BAL_new_individual(3).V.^2;
        BAL_new_individual(3).CL = BAL_new_individual(3).CL.*V.^2 ./ BAL_new_individual(3).V.^2;
        BAL_new_individual(3).CY = BAL_new_individual(3).CY.*V.^2 ./ BAL_new_individual(3).V.^2;
        BAL_new_individual(3).CMy = BAL_new_individual(3).CMy.*V.^2 ./ BAL_new_individual(3).V.^2;
        

        CN_Tail = BAL_new.CY;

        % CN_Tail = CN_Model;

        % Streamline Curvature correction (lift interference) for the vertical tail due to its
        % own circulation
        Tau_2 = 0.07;           % Taper = 0.68, MAC = 0.17, Lt = 0.085, Lt/2H = 0.11, Fig. 10.37 in the Barlow book
        Delta = 0.12;          % Elliptic loading for a closed elliptic jet, Fig. 10.32 in the Barlow book
        Slope_CN_vs_Beta = mean(diff(CN_Tail)) ./ mean(diff(BAL_new.AoS));
        Delta_Beta_Total_SC = (1+Tau_2)*Delta*(A_VTail/A_tunnel)*BAL_new.CY;
        Delta_Cm = 0.125*Delta_Beta_Total_SC*Slope_CN_vs_Beta;               % moment = arm*db*dCn/db (Barlow page 377)
        BAL_new.CD = BAL_new.CD + Delta*(A_VTail/A_tunnel)*BAL_new.CL.^2;           % Equation 2.10 pre-test report

        BAL_new.CMy = BAL_new.CMy + Delta_Cm;
        
        % Downwash correction
        Beta_u = BAL_new.AoS;
        Beta_downwash = Beta_u+(Delta*(A_VTail/A_tunnel)*CN_Tail*57.3);

        % Correct beta for both
        BAL_new.AoS = Delta_Beta_Total_SC + Beta_downwash;
        
        
        % Only slipstream correction
        CN_Tail = BAL_new_individual(4).CY;
        Slope_CN_vs_Beta = mean(diff(CN_Tail)) ./ mean(diff(BAL_new_individual(4).AoS));
        Delta_Beta_Total_SC = (1+Tau_2)*Delta*(A_VTail/A_tunnel)*BAL_new_individual(4).CY;
        Delta_Cm = 0.125*Delta_Beta_Total_SC*Slope_CN_vs_Beta;               % moment = arm*db*dCn/db (Barlow page 377)
        BAL_new_individual(4).CD = BAL_new_individual(4).CD + Delta*(A_VTail/A_tunnel)*BAL_new_individual(4).CL.^2;           % Equation 2.10 pre-test report

        BAL_new_individual(4).CMy = BAL_new_individual(4).CMy + Delta_Cm;
        BAL_new_individual(4).AoS = BAL_new_individual(4).AoS + Delta_Beta_Total_SC;
        
        % Downwash correction
        Beta_u = BAL_new_individual(5).AoS;
        Beta_downwash = Beta_u+(Delta*(A_VTail/A_tunnel)*CN_Tail*57.3);
        BAL_new_individual(5).AoS = Beta_downwash;
        
%         % Buoyancy correction, Equation 2.2 pre-test report
%         dp_by_dv = ;                        % from pressure taps along the fuselage
%         del_D = dp_by_dv * V_total;
%         BAL_new.CD = BAL_new.CD - del_D/(0.5*BAL_new.rho.*BAL_new.V.^2.*BAL_new.S);

        % Update original file
        BAL.windOn.(field{:}) = BAL_new;
        BAL_new_individual_return.(field{:}) = BAL_new_individual;
    end
end
