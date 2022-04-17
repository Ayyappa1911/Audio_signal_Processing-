
% Mixing two recorded signals in such a way that they are recorded at a
% time.

% This file contains the code for :-
% 1) Converting the recorded signals from mp3 format into .wav
% 2) Seperating the channels of the recorded signals.
% 3) Taking one channel from each signal and add them as channels to
%    generate a signal which looks like they are recorded at a time.
% 4) Similar to third point all channels of two signals were added
%    together.

% extracting data from first recorded mp3 file.
signal = audioread("baby elephant 60.mp3");
info = audioinfo("baby elephant 60.mp3");
Fs = info.SampleRate;
display(info);

% Converting the extracted data into .wav extension.
audiowrite("wav_inputs/input_Elephant.wav",signal,Fs);


signal = audioread("wav_inputs/input_Elephant.wav");

% Signals recorded in voice recorder will have two channels which vary in
% Power and Amplitude.
% 
% Here, We are seperating two channels and saving them again as .wav files.
column1 = signal(:, 1);
column2 = signal(:, 2);

audiowrite("input1_cnls/test1.wav",column1,Fs);
audiowrite("input1_cnls/test2.wav",column2,Fs);

%Similar process is done for second recorded signal.
signal1 = audioread("pink_panther.mp3");
info1 = audioinfo("pink_panther.mp3");
Fs1 = info1.SampleRate;
display(info1);

audiowrite("wav_inputs/input_pink_panther.wav",signal1,Fs);

signal1 = audioread("wav_inputs/input_pink_panther.wav");

column1_1 = signal1(:, 1);
column2_1 = signal1(:, 2);

audiowrite("input2_cnls/test1_1.wav",column1_1,Fs);
audiowrite("input2_cnls/test2_1.wav",column2_1,Fs);

%display(length(column2));
%display(length(column2_1));

a = length(column2);
b = length(column2_1);

% we will chose one of the channels from each signal to mix.
%
% Before mixing we ensure that they are of same length(same time duration)
% by padding zeroes at the end.
%
% Here, I have used the signals having higher power(just for convenience.)
if(a-b>0)
    for a = 1 : 1 : a-b
        column2_1(end+1) = 0;
        column1_1(end+1) = 0;

    end
else
    for a = 1 : 1 : b-a
        column2(end+1) = 0;
        column1(end+1) = 0;
    end
end


%display(length(column2));
%display(length(column2_1)); 
%display(length(column1));
%display(length(column1_1));

% We add the Signals in the form of channels taking one as left channel and
% other as right channel.
Output_signal = [column2,column2_1];

audiowrite("Output/Output_signal.wav",Output_signal,Fs);

info = audioinfo("Output/Output_signal.wav");

display(info);


% For the second case, we are adding add the four channels of recorded
% signal.
Output_signal1 = [column2,column2_1,column1,column1_1];

audiowrite("Output/Output_signal1.wav",Output_signal1,Fs);

info1 = audioinfo("Output/Output_signal1.wav");

display(info1);
