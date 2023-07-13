function [acw_0, acw_50, acf, lags] = acw_3d(x, fs)
%% Calculate ACW-0 and ACW-50 of ACF in a 3D array of channels x time x trials
% Get one ACF per trial, average ACFs, do the calculation

nroi = size(x, 1);
ntime = size(x, 2);
ntrial = size(x, 3);

decayrate = zeros(nroi, 1);
ACF = zeros(nroi, ntime, ntrial);

for i = 1:ntrial
    for j = 1:nroi
        [~, ~, ACF(j, :, i), lags] = acw(x(j, : ,i), fs);
    end
end

acf = mean(ACF, 3);

[~, ACW_50_i] = max(acf<=0.5, [], 2);
acw_50 = ACW_50_i / fs; % Divide by fs to convert to seconds
[~, ACW_0_i] = max(acf<=0, [], 2);
acw_0 = ACW_0_i / fs;
lags = lags / fs;
end