function [acw_0, acw_50, acf, lags] = acw_windowed(x, fs, window, overlap, isplot)
    %% Calculate ACF in windows, average, calculate acw-0 and 50 in averaged ACF
    % Default settings are as in Honey et al 2012 (20 second windows with 50%
    % overlap)
    % Similar to Welch method really
    % Input:
    %   x: Signal (1D vector)
    %   fs: Sampling rate of signal
    %   window: Window size in seconds
    %   overlap: Overlap in percentage (e.g. 50 means 50%)
    %   isplot: Logical argument to plot or not
    % Output:
    %   acw_0, acw_50, acf, lags: Pretty intuitive
    % If you want to use with default settings, use as
    % acw_windowed(x, fs)
    % Authored by Yasir Çatal aka Duodenum
    
    
    window=window*fs;
    ACF = [];
    % Do the sliding window, adapted from ACW_estimation from dynameas toolbox
    ii=1;   % Windows counter
    while true
        % Begining and ending of the current window
        SWindow=[1+(ii-1)*window*overlap/100, window+(ii-1)*window*overlap/100];
        % Chek if index exceeds vector dimensions. If so, break!
        if SWindow(2)>=length(x), break; end
        % ACF computation into the window (normalized between -1 and 1)
        [~, ~, ACF(ii,:), lags] = acw(x(SWindow(1):SWindow(2)), fs);
        % Next window
        ii=ii+1;
    end
    
    % Average ACF of windows to get acf
    acf = mean(ACF,1);
    % Rest of the calculations are the same as acw function
    [~, ACW_50_i] = max(acf<=0.5);
    acw_50 = ACW_50_i / fs; % Divide by fs to convert to seconds
    [~, ACW_0_i] = max(acf<=0);
    acw_0 = ACW_0_i / fs;
    
    if isplot
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