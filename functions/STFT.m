function [time, freq, time_freq_mat] = STFT(audio, Fs, window_time)
    window_length = Fs*window_time;
    window_num = floor(length(audio)/(window_length/2));
    time_freq_mat = zeros(1+floor(window_length/2),window_num-1);
    % calculating fft using an overlapping sliding window
    temp = zeros(1+floor(window_length/2),1);
    for i = 1:window_num-1
        %%% time_freq_mat
        first = 1+(i-1)*floor(window_length/2);
        last = (i-1)*floor(window_length/2)+window_length;
        temp = FFT(audio(first:last));
        time_freq_mat(:,i) = temp + rand(1+floor(window_length/2),1) ./ 1000;
    end
    % time and freq vectors
    time = (window_time/2)*(1:(window_num-1));
    freq = Fs * (0:floor(window_length/2)) / window_length;
end