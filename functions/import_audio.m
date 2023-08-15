function [downsampled_Fs, resampled_audio] = import_audio(path, song_num, format)
    % import the audio 
    [audio, Fs] = audioread([path, 'music', num2str(song_num), format]);
    % getting mean over right and left channels
    audio_size = size(audio);
    if audio_size(2) == 2
        shifted_audio = (audio(:,1) + audio(:,2))/2;
    else
        shifted_audio = audio;
    end
    %%% audioMono
    % downsample the audio to 8 KHz
    downsampled_Fs = 8000;
    Len = length(shifted_audio);
    time1 = 0 : 1/Fs : 1/Fs * (Len-1);
    main_audio = timeseries(shifted_audio,time1);
    desired_time = 0 : 1/downsampled_Fs : 1/Fs * (Len-1);
    downsampled_audio = resample(main_audio,desired_time);
    resampled_audio = downsampled_audio.Data;
    %%% resampled_audio
end