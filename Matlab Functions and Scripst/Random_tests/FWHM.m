function [FWHM_samples, FWHM_hz] = FWHM(data,freq_v)


halfMaxValue = (min(data) + max(data)) / 2; % Find the half max value.

% Find indexes of power where the power first and last is above the half max value.
leftIndex = find(data >= halfMaxValue, 1, 'first');
rightIndex = find(data >= halfMaxValue, 1, 'last');

FWHM_samples = rightIndex-leftIndex + 1; % FWHM in indexes.

% OR, if you have an x vector
FWHM_hz = freq_v(rightIndex) - freq_v(leftIndex);
end