%% Part 7
clear; close all; clc;
addpath('functions/');
addpath('database/');
addpath('musics/');
addpath('test_musics/');

database = load('database/database.mat').database;

% importing an audio
path = 'test_musics/'; % test musics path
format = '.wav';
check = zeros(50,1);
for song_num = 1:10
[downsampled_Fs, audioMono] = import_audio(path, song_num, format);

% creating the time-freq matrix of the audio using fft and an overlapping sliding window with the length of "window_time"
window_time = 0.1;
[time, freq, time_freq_mat] = STFT(audioMono, downsampled_Fs, window_time);

% a full screen figure for plots
figure('Units','normalized','Position',[0 0 1 1])

% plotting the stft
subplot(1,2,1);
pcolor(time, freq, time_freq_mat);
shading interp
colorbar;
xlabel('time(s)');
ylabel('frequency(Hz)');
title('STFT(dB)');

% finding the anchor points from time_freq_mat using a sliding window with the size of 2dt*2df
df = floor(0.1*size(time_freq_mat, 1)/4);
dt = 2/window_time;
% finding anchor points
anchor_points = find_anchor_points(time_freq_mat, dt, df);
% plotting the anchor points
subplot(1,2,2);
scatter(time(anchor_points(:, 2)), freq(anchor_points(:, 1)),'x');
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
title("anchor points",'interpreter','latex');
xlim([time(1) time(end)]);
ylim([freq(1) freq(end)]);
grid on; grid minor;

% creating the hash tags using a window with the size of dt*2df for each anchor point
df_hash = floor(0.1*size(time_freq_mat,1));
dt_hash = 20/window_time;
% creating hash-keys and hash-values for each pair of anchor points
% Key format: (f1*f2*(t2-t1)) 
% Value format: (song_name*time_from_start)
[hash_key, hash_value] = create_hash_tags(anchor_points, df_hash, dt_hash, 0);

list = []; 

% searching for found hash-keys in the database
for i = 1:length(hash_key)
    key_tag = [num2str(hash_key(i, 1)), '*', num2str(hash_key(i, 2)), '*', num2str(hash_key(i, 3))];
    if (isKey(database, key_tag))
        temp1 = split(database(key_tag),'+');
        for j = 1:length(temp1)
            temp2 = split(temp1{j},'*');
            list = [list; [str2num(temp2{1}),str2num(temp2{2}),hash_value(i,2)]];
        end
    end
end
score = scoring(list);
detected_song = score(1,1);

start_time = [];
for i=1:length(list)
    if list(i,1) == detected_song
        start_time = [start_time;list(i,2)];
    end
end
fprintf("\ndetected song is : %d\n",detected_song);
fprintf("begin time of test music %d (s) : %f \n----------------------\n",song_num,median(start_time)/20);
end


%% Part 8
clear; close all; clc;
addpath('functions/');
addpath('database/');
addpath('musics/');
addpath('test_musics/');
database = load('database/database.mat').database;
% importing an audio
path = 'musics/'; % test musics path
format = '.mp3';
test_songs = 19;
for song_num = test_songs
[downsampled_Fs, audioMono] = import_audio(path, song_num, format);
detected_songs = zeros(100,1);
similarity = zeros(100,1);
begin = randi([1,length(audioMono)-20*downsampled_Fs]);
main_sample = audioMono(begin:begin+20*downsampled_Fs);
sample = main_sample;
for k=1:100
    disp(k);
    snr = 0.1*k;
    sample = awgn(main_sample,snr);
    window_time = 0.1;
    [time, freq, time_freq_mat] = STFT(sample, downsampled_Fs, window_time);

% finding the anchor points from time_freq_mat using a sliding window with the size of 2dt*2df
df = floor(0.1*size(time_freq_mat, 1)/4);
dt = 2/window_time;
% finding anchor points
anchor_points = find_anchor_points(time_freq_mat, dt, df);
% creating the hash tags using a window with the size of dt*2df for each anchor point
df_hash = floor(0.1*size(time_freq_mat,1));
dt_hash = 20/window_time;
% creating hash-keys and hash-values for each pair of anchor points
% Key format: (f1*f2*(t2-t1)) 
% Value format: (song_name*time_from_start)
[hash_key, hash_value] = create_hash_tags(anchor_points, df_hash, dt_hash, 0);

list = []; 
% searching for found hash-keys in the database
for i = 1:length(hash_key)
    key_tag = [num2str(hash_key(i, 1)), '*', num2str(hash_key(i, 2)), '*', num2str(hash_key(i, 3))];
    if (isKey(database, key_tag))
        temp1 = split(database(key_tag),'+');
        for j = 1:length(temp1)
            temp2 = split(temp1{j},'*');
            list = [list; [str2num(temp2{1}),str2num(temp2{2}),hash_value(i,2)]];
        end
    end
end
score = scoring(list);
detected_songs(k) = score(1,1);
similarity(k) = score(1,2);
end
SNR = 0.1:0.1:10;
figure;
subplot(2,1,1);
plot(SNR,similarity,LineWidth=1);
xlabel("$snr$",Interpreter="latex",Fontsize=18);
ylabel("$P$",Interpreter="latex",Fontsize=18);
title("Probability of similarity",Interpreter="latex",Fontsize=18);
grid minor;

subplot(2,1,2);
scatter(SNR,detected_songs,'filled');
xlabel("$snr$",Interpreter="latex",Fontsize=18);
ylabel("$n$",Interpreter="latex",Fontsize=18);
title("detected song number",Interpreter="latex",Fontsize=18);
grid minor;
sgtitle("Song number "+song_num,fontSize=20);
end


%% Part 9
clear; close all; clc;
addpath('functions/');
addpath('database/');
addpath('musics/');
addpath('test_musics/');
database = load('database/database.mat').database;
% importing an audio
path = 'musics/'; % test musics path
format = '.mp3';
test_songs = randperm(50,5);
help = zeros(1,100);
prob = zeros(1,10);
for song_num = test_songs
[downsampled_Fs, audioMono] = import_audio(path, song_num, format);
for k=0:9
    snr = 10*(1 + k);
for l=1:100
    begin = randi([1,length(audioMono)-20*downsampled_Fs]);
    sample = audioMono(begin:begin+20*downsampled_Fs);
    sample = awgn(sample,snr);
    window_time = 0.1;
[time, freq, time_freq_mat] = STFT(sample, downsampled_Fs, window_time);

% finding the anchor points from time_freq_mat using a sliding window with the size of 2dt*2df
df = floor(0.1*size(time_freq_mat, 1)/4);
dt = 2/window_time;
% finding anchor points
anchor_points = find_anchor_points(time_freq_mat, dt, df);
% creating the hash tags using a window with the size of dt*2df for each anchor point
df_hash = floor(0.1*size(time_freq_mat,1));
dt_hash = 20/window_time;
% creating hash-keys and hash-values for each pair of anchor points
% Key format: (f1*f2*(t2-t1)) 
% Value format: (song_name*time_from_start)
[hash_key, hash_value] = create_hash_tags(anchor_points, df_hash, dt_hash, 0);

list = []; 
% searching for found hash-keys in the database
for i = 1:length(hash_key)
    key_tag = [num2str(hash_key(i, 1)), '*', num2str(hash_key(i, 2)), '*', num2str(hash_key(i, 3))];
    if (isKey(database, key_tag))
        temp1 = split(database(key_tag),'+');
        for j = 1:length(temp1)
            temp2 = split(temp1{j},'*');
            list = [list; [str2num(temp2{1}),str2num(temp2{2}),hash_value(i,2)]];
        end
    end
end
score = scoring(list);
help(l) = score(1,2);
end
prob(k+1) = mean(help);
end
SNR = 10:10:100;
figure;
plot(SNR,prob,LineWidth=1);
xlabel("$snr$",Interpreter="latex",Fontsize=18);
ylabel("$P$",Interpreter="latex",Fontsize=18);
title("Probability of similarity for Song "+song_num,Interpreter="latex",Fontsize=18);
grid minor;
end
%% Part 10
clear; close all; clc;

addpath('functions/');
addpath('database/');
addpath('musics/');
addpath('test_musics/');

database = load('database/database.mat').database;
% importing an audio
path = 'test_musics/'; % test musics path
format = '.wav';
for song_num = [11,12]
[downsampled_Fs, audioMono] = import_audio(path, song_num, format);
% creating the time-freq matrix of the audio using fft and an overlapping sliding window with the length of "window_time"
window_time = 0.1;
[time, freq, time_freq_mat] = STFT(audioMono, downsampled_Fs, window_time);

% a full screen figure for plots
figure('Units','normalized','Position',[0 0 1 1])

% plotting the stft
subplot(1,2,1);
pcolor(time, freq, time_freq_mat);
shading interp
colorbar;
xlabel('time(s)');
ylabel('frequency(Hz)');
title('STFT(dB)');

% finding the anchor points from time_freq_mat using a sliding window with the size of 2dt*2df
df = floor(0.1*size(time_freq_mat, 1)/4);
dt = 2/window_time;
% finding anchor points
anchor_points = find_anchor_points(time_freq_mat, dt, df);
% plotting the anchor points
subplot(1,2,2);
scatter(time(anchor_points(:, 2)), freq(anchor_points(:, 1)),'x');
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
title("anchor points",'interpreter','latex');
xlim([time(1) time(end)]);
ylim([freq(1) freq(end)]);
grid on; grid minor;

% creating the hash tags using a window with the size of dt*2df for each anchor point
df_hash = floor(0.1*size(time_freq_mat,1));
dt_hash = 20/window_time;
% creating hash-keys and hash-values for each pair of anchor points
% Key format: (f1*f2*(t2-t1)) 
% Value format: (song_name*time_from_start)
[hash_key, hash_value] = create_hash_tags(anchor_points, df_hash, dt_hash, 0);
list = []; 

% searching for found hash-keys in the database
for i = 1:length(hash_key)
    key_tag = [num2str(hash_key(i, 1)), '*', num2str(hash_key(i, 2)), '*', num2str(hash_key(i, 3))];
    if (isKey(database, key_tag))
        temp1 = split(database(key_tag),'+');
        for j = 1:length(temp1)
            temp2 = split(temp1{j},'*');
            list = [list; [str2num(temp2{1}),str2num(temp2{2}),hash_value(i,2)]];
        end
    end
end

score = scoring(list);
detected_song = score(1,1);
start_time = [];
for i=1:length(list)
    if list(i,1) == detected_song
        start_time = [start_time;list(i,2)];
    end
end
fprintf("\ndetected song is : %d\n",detected_song);
fprintf("begin time of test music %d (s) : %f \n----------------------\n",song_num,median(start_time)/20);
end


%% Part 11
clear; close all; clc;
addpath('functions/');
addpath('database/');
addpath('musics/');
addpath('test_musics/');
database = load('database/database.mat').database;
% importing an audio
path = 'musics/'; % test musics path
format = '.mp3';

[downsampled_Fs, audioMono1] = import_audio(path, 11, format);
[downsampled_Fs2, audioMono2] = import_audio(path, 22, format);

begin1 = randi([1,length(audioMono1)-20*downsampled_Fs]);
main_sample1 = audioMono1(begin1:begin1+20*downsampled_Fs);

begin2 = randi([1,length(audioMono2)-20*downsampled_Fs]);
main_sample2 = audioMono2(begin2:begin2+20*downsampled_Fs);

k = 1;
prob11 = zeros(11,1);
prob22 = zeros(11,1);
detected_songs = zeros(11,1);
for alpha = 0:0.1:1
   sample = alpha * main_sample1 + (1-alpha) * main_sample2;
   window_time = 0.1;
   [time, freq, time_freq_mat] = STFT(sample, downsampled_Fs, window_time);

% a full screen figure for plots
figure('Units','normalized','Position',[0 0 1 1])

% plotting the stft
subplot(1,2,1);
pcolor(time, freq, time_freq_mat);
shading interp
colorbar;
xlabel('time(s)');
ylabel('frequency(Hz)');
title('STFT(dB)');

% finding the anchor points from time_freq_mat using a sliding window with the size of 2dt*2df
df = floor(0.1*size(time_freq_mat, 1)/4);
dt = 2/window_time;
% finding anchor points
anchor_points = find_anchor_points(time_freq_mat, dt, df);
% plotting the anchor points
subplot(1,2,2);
scatter(time(anchor_points(:, 2)), freq(anchor_points(:, 1)),'x');
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
title("anchor points",'interpreter','latex');
xlim([time(1) time(end)]);
ylim([freq(1) freq(end)]);
grid on; grid minor;

% creating the hash tags using a window with the size of dt*2df for each anchor point
df_hash = floor(0.1*size(time_freq_mat,1));
dt_hash = 20/window_time;
% creating hash-keys and hash-values for each pair of anchor points
% Key format: (f1*f2*(t2-t1)) 
% Value format: (song_name*time_from_start)
[hash_key, hash_value] = create_hash_tags(anchor_points, df_hash, dt_hash, 0);
list = []; 
% searching for found hash-keys in the database
for i = 1:length(hash_key)
    key_tag = [num2str(hash_key(i, 1)), '*', num2str(hash_key(i, 2)), '*', num2str(hash_key(i, 3))];
    if (isKey(database, key_tag))
        temp1 = split(database(key_tag),'+');
        for j = 1:length(temp1)
            temp2 = split(temp1{j},'*');
            list = [list; [str2num(temp2{1}),str2num(temp2{2}),hash_value(i,2)]];
        end
    end
end
score = scoring(list);
check11 = false;
check22 = false;
index11 = 0;
index22 = 0;
counter = 1;
temp = 0;
while check11 ~= true || check22 ~= true
    if score(counter,1) == 11
        check11 = true;
        index11 = counter;
    end
    if score(counter,1) == 22
        check22 = true;
        index22 = counter;
    end
    counter = counter + 1;
end
prob11(k) = score(index11,2);
prob22(k) = score(index22,2);
detected_songs(k) = score(1,1);
k = k+1;
end
ALPHA = 0:0.1:1;
figure;
subplot(2,1,1);
plot(ALPHA,prob11,ALPHA,prob22,LineWidth=1);
xlabel("$\alpha$",Interpreter="latex",Fontsize=18);
ylabel("$P$",Interpreter="latex",Fontsize=18);
title("Probability of similarity",Interpreter="latex",Fontsize=18);
grid minor;
legend(["music 11","music 22"]);

subplot(2,1,2);
scatter(ALPHA,detected_songs,'filled');
xlabel("$\alpha$",Interpreter="latex",Fontsize=18);
ylabel("$n$",Interpreter="latex",Fontsize=18);
title("detected song number",Interpreter="latex",Fontsize=18);
grid minor;
sgtitle("Song 11 and Song 22",fontSize=20);