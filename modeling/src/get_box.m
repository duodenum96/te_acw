function [boxcar, offsets] = get_box(starttime, iti, dur, amplitude, endtime, dt, ntime)
    fs = (1 / dt) * 1000; % ms to s
    boxcar = zeros(1, ntime);
    offsets = (starttime*fs):round(iti*fs):(endtime*fs);
    
    for i = 1:length(offsets)
        o = offsets(i);
        boxcar(o:round(o+dur*fs)) = amplitude;
    end
end