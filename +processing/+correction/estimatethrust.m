function [Tc] = estimatethrust(velocity,advanceRatio)
%ESTIMATETHRUST Estimates thrust coefficient based on emperical data

[parentDir, ~, ~] = fileparts(mfilename("fullpath"));
ThrustData = load(fullfile(parentDir, 'thrustdata.mat'));

try
    field = sprintf('V%2d', round(velocity));
    thrustCoefficients = ThrustData.(field).Tc;
    advanceRatios = ThrustData.(field).J;
catch
    error('No data exists for V=%f', velocity)
end

% Tc = interp1(advanceRatios, thrustCoefficients, advanceRatio, 'makima');
p = polyfit(advanceRatios,thrustCoefficients,2);
Tc = polyval(p, advanceRatio);


figure('Name', 'Thrust Interpolation')
plot(advanceRatios, thrustCoefficients, '-o', 'DisplayName', 'Data')
hold on
plot(advanceRatio, Tc, '-or', 'DisplayName', 'Interpolated');
legend('Location', 'Best')
xlabel('Advance Ratio')
ylabel('Thrust Coefficient')
title(...
    sprintf('T_c vs. J, V = %d [m/s]', round(velocity))...
);
hold off

end

