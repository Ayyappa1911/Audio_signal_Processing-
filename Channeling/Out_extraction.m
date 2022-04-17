
% Extracting the signals from the mixed signals.

% This file contains the code for the seperating the channels of 4 channel
% mixed signal to seperate all the recorded signals.

% Data extraction from the four channel mixed output signal.
[Output_signal1,Out_Fs] = audioread("Output/Output_signal1.wav");
info = audioinfo("Output/Output_signal1.wav");
display(info);

% Seperating the data into channels.
output1_1 = Output_signal1(:,1);
output1_2 = Output_signal1(:,2);
output1_3 = Output_signal1(:,3);
output1_4 = Output_signal1(:,4);

% Writing the extracted channels in .wav extension files.
audiowrite("splt_audio/Output1_1.wav",output1_1,Out_Fs);
audiowrite("splt_audio/Output1_2.wav",output1_2,Out_Fs);
audiowrite("splt_audio/Output1_3.wav",output1_3,Out_Fs);
audiowrite("splt_audio/Output1_4.wav",output1_4,Out_Fs);