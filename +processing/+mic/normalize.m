function [MicData] = normalize(MicData, BalData, D, nBlades)
    %% Adds sound pressures normalized by the propeller thrust into MicData
    %
    % Args:
    %   MicData: MicData structure
    %   BalData: Uncorected Balance Data
    %   D: Prop Diameter

    propArea = pi * D / 4;
    for field = fieldnames(MicData)'
        thrust = BalData.windOn.(field{:}).FT1;
        soundPressureData = MicData.(field{:}).PXX;
        for i = 1:length(soundPressureData)
            soundPressure = sqrt(soundPressureData{i});
            % TODO consider using RMS pressure
            thrustPressure = thrust(i) / propArea;
            if thrustPressure <= 0
                thrustPressure = 1;
            end
            thrustPressure(thrustPressure <= 0) = 1;  % Correcting for negative thrust
            MicData.(field{:}).Pn{i} = soundPressure ./ thrustPressure;
        end
    end
end