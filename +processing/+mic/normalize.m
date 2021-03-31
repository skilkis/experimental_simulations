function [MicData] = normalize(MicData, BalData, D)
    %% Adds sound pressures normalized by the propeller thrust into MicData
    %
    % Args:
    %   MicData: MicData structure
    %   BalData: Uncorected Balance Data
    %   D: Prop Diameter

    totalPropArea = pi * D / 2;
    for field = fieldnames(MicData)'
        thrust = BalData.windOn.(field{:}).FT1 + BalData.windOn.(field{:}).FT2;
        soundPressureData = MicData.(field{:}).PXX;
        for i = 1:length(soundPressureData)
            soundPressure = sqrt(soundPressureData{i});
            % TODO consider using RMS pressure
            thrustPressure = thrust(i) / totalPropArea;
            if thrustPressure <= 0
                thrustPressure = 1;
            end
            thrustPressure(thrustPressure <= 0) = 1;  % Correcting for negative thrust
            MicData.(field{:}).Pn{i} = soundPressure ./ thrustPressure;
        end
    end
end