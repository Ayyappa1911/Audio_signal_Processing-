
% Finding the spectrum of input signals and output signals.

% This file contains the code for :-
% 1) Finding the Spectrum of recorded signals.
% 2) Finding the Output signal formed by merging one channel from each of the recorded
%    signals.
% 3) Finding the Output signal formed by merging all the channels of recorded signals. 


% Data from first recording is exracted.
signal = audioread("wav_inputs/input_Elephant.wav");
info = audioinfo("wav_inputs/input_Elephant.wav");
Fs = info.SampleRate;
ts = info.Duration * Fs;

f = (-ts/2:ts/2-1)*(Fs/ts);

subplot(4,1,1);
% Plotting the spectrum of first recorded signal.
plot(f,abs(fftshift(fft(signal))/length(signal)));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|X(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|X(f)|');
ylim([0,0.004]);
title("Spectrum of input1 signal");


% Data from the second recording is extracted.
signal = audioread("wav_inputs/input_pink_panther.wav");
info = audioinfo("wav_inputs/input_pink_panther.wav");
Fs = info.SampleRate;
ts = info.Duration * Fs;

f = (-ts/2:ts/2-1)*(Fs/ts);

subplot(4,1,2);
% Plotting the spectrum of second recorded signal.
plot(f,abs(fftshift(fft(signal))/length(signal)));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|X(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|X(f)|');
ylim([0,0.004]);
title("Spectrum of input2 signal");

% Data extracted from first output file.
signal = audioread("Output/Output_signal.wav");
info = audioinfo("Output/Output_signal.wav");
Fs = info.SampleRate;
ts = info.Duration * Fs;

f = (-ts/2:ts/2-1)*(Fs/ts);

subplot(4,1,3);
% plotting the spectrum of the first Output signal.
plot(f,abs(fftshift(fft(signal))/length(signal)));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|X(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|X(f)|');
ylim([0,0.004]);
title("Spectrum of 2 channel output signal");


% Data extracted from second output file.
signal = audioread("Output/Output_signal1.wav");
info = audioinfo("Output/Output_signal1.wav");
Fs = info.SampleRate;
ts = info.Duration * Fs;

f = (-ts/2:ts/2-1)*(Fs/ts);

subplot(4,1,4);
% plotting the spectrum of the second Output signal.
plot(f,abs(fftshift(fft(signal))/length(signal)));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|X(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|X(f)|');
ylim([0,0.004]);
title("Spectrum of 4 channel output signal");