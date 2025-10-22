function [ppm,peak_locs] = calcClockDrift(t_nco_values)
% CALCCLOCKDRIFT Calculate the clock drift in parts per million (ppm) from t_nco values.
%
% Parameters:
%   t_nco_values : matrix or vector
%       An array containing the relative time delay values provided to the NCO.
%
% Returns:
%   ppm : scalar or vector
%       The clock deviation in parts per million (ppm).

    % Ensure t_nco_values is a column vector if input is 1D
    if isvector(t_nco_values)
        t_nco_values = t_nco_values(:);
    end

    % Calculate timing error by removing the mean
    timingError = t_nco_values - mean(t_nco_values, 1);

    % Create time index vector
    t = (0:size(timingError, 1)-1)';

    % Get number of modes (columns)
    nModes = size(t_nco_values, 2);
    ppm = zeros(1, nModes);

    for indMode = 1:nModes
        % Calculate absolute difference of timing error
        diff_timingError = abs(diff(timingError(:, indMode)));

        % Find peaks with minimum height of 0.5
        [peaks, ~] = findpeaks(diff_timingError, 'MinPeakHeight', 0.4);

        if length(peaks) < 2
            warning('Not enough peaks found for mode %d. Cannot calculate period.', indMode);
            ppm(indMode) = NaN;
            continue;
        end

        % Find peak locations
        peak_locs = find(ismember(diff_timingError, peaks));

        % Calculate mean period between peaks
        mean_period = mean(diff(t(peak_locs)));

        % Calculate frequency
        fo = 1 / mean_period;

        % Calculate ppm with sign based on mean t_nco_values
        ppm(indMode) = sign(mean(t_nco_values(:, indMode))) * fo * 1e6;
    end

    % If single mode, return scalar instead of vector
    if nModes == 1
        ppm = ppm(1);
    end
end