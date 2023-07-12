function [fr3d, time2d] = trialize(v_E, time, offsets, dur, dt)
    %% Turn the 2d firing rates to 3d roi x time x trial array for erps
    % discard the first 50 seconds of offsets (consistent with chaudhuri.m)
    offsets(offsets < 50*1000) = [];
    offsets = offsets - 50*1000;
    offsets(offsets == 0) = [];
    
    % do the rest of the stuff
    fs = (1/dt)*1000;
    fr3d = zeros([size(v_E, 1), dur*fs+1, length(offsets)]);
    time2d = zeros([dur*fs+1, length(offsets)]);
    
    for i = 1:length(offsets)
        o = offsets(i);
        fr3d(:, :, i) = v_E(:, o:round(o+dur*fs));
        time2d(:, i) = time(o:round(o+dur*fs));
    end
    end