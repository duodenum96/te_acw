function [acw_0, acw_50, acf, lags] = acw(x, fs, varargin)
    %% Calculate ACW-0 and ACW-50
    % This function doesn't follow Honey et al.'s version (averaging ACW
    % values in 50% overlapping 20s segments). Instead, it's a direct version that
    % is more suitable for fMRI data due to it's low temporal resolution.
    % x is input time series. If 'plot' is an input, then plots ACW_0 and
    % ACW_50. fs is sampling rate in seconds. if not specified, fs defaults to
    % 1.
    %
    % Authored by Yasir Çatal a.k.a. Duodenum
    
    % Check if econometrics toolbox installed
    a = ver('econ');
    
    if ~isempty(a) % Method using autocorr
        [acf, lags] = autocorr(x, length(x)-1);
        [~, ACW_50_i] = max(acf<=0.5);
        acw_50 = ACW_50_i / fs; % Divide by fs to convert to seconds
        [~, ACW_0_i] = max(acf<=0);
        acw_0 = ACW_0_i / fs;
        lags = lags / fs;
    else % Method using xcorr
        [acf, lags] = xcorr(x, 'coeff');
        index = find(acf == 1); % Get rid of the left-side
        acf = acf(index:end); lags = lags(index:end);
    
        % Rest of the calculations are same as above
    
        [~, ACW_50_i] = max(acf<=0.5);
        acw_50 = ACW_50_i / fs; % Divide by fs to convert to seconds
        [~, ACW_0_i] = max(acf<=0);
        acw_0 = ACW_0_i / fs;
        lags = lags / fs;
    end
    
    if CheckInput(varargin,'plot')
        plot(lags,acf,'k')
        xlim([0 max(lags)])
        hold on
        area(lags(1:ACW_50_i), acf(1:ACW_50_i),'FaceColor','r','FaceAlpha',0.3);
        area(lags(1:ACW_0_i), acf(1:ACW_0_i),'FaceColor','m','FaceAlpha',0.3);
        title(['ACW-0 = ', num2str(acw_0, '%.1f'), ' ACW-50 = ', num2str(acw_50, '%.1f')])
        xlabel('Lags (s)')
        ylabel('Autocorrelation')
    end
    end