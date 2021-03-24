function BAL = subtractzero(BAL,BAL0)
    nRuns = length(BAL.run);

    % Rounding measured Angle of Attack and Angle of Sideslip
    AoAr = round(BAL.AoA);
    AoSr = round(BAL.AoS);
    AoA0r = round(BAL0.AoA);
    AoS0r = round(BAL0.AoS);

        function [B] = find_zero_measurements(AoA, AoS)
            % Find correct zero measurement using AoA and AoS scalar input
            match = (AoA == AoA0r) & (AoS == AoS0r);
            if all(match == 0)
                error([...
                    'No calibration data found'...
                    'for AoA=%0.2f, AoS=%0.2f'...
                ], AoA, AoS)
            end
            B = [...
                BAL0.B1(match),...
                BAL0.B2(match),...
                BAL0.B3(match),...
                BAL0.B4(match),...
                BAL0.B5(match),...
                BAL0.B6(match),...
            ];
        end

    BAL0.intp = zeros(nRuns, 6);
    for i = 1:nRuns
        try
            BAL0.intp(i,:) = find_zero_measurements(AoAr(i), AoSr(i));
        catch
            fprintf('Run %d has to applicable zero measurement\n',i)
        end
    end

    % Subtract zero measurement data from the measured data
    BAL.Bzeroed = BAL.B-BAL0.intp;

end