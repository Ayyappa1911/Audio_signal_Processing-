
% Method 2 : Mixing the recorded signals and extracting the signals using
%  ICA - Independent Component Analysis.(Blind source Seperation).

% This file contains the code for mixing the recorded signals in the form
% of a1+*S1(t)+ a2*S2(t). And seperating using ICA.

signal = audioread("baby elephant 60.mp3");
info = audioinfo("baby elephant 60.mp3");
Fs = info.SampleRate;

display(info);

audiowrite("wav_inputs/input_Elephant.wav",signal,Fs);

signal = audioread("wav_inputs/input_Elephant.wav");
info = audioinfo("wav_inputs/input_Elephant.wav");
Fs = info.SampleRate;
ts = info.Duration * Fs;
ts1 = ts;
f = (-ts/2:ts/2-1)*(Fs/ts);

column1 = signal(:, 1);
column2 = signal(:, 2);

audiowrite("input1_cnls/test1.wav",column1,Fs);
audiowrite("input1_cnls/test2.wav",column2,Fs);


subplot(6,1,1);
plot(f,abs(fftshift(fft(signal))/length(signal)));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|X1(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|X1(f)|');
ylim([0,0.012]);
title("Spectrum of input1 signal");




signal1 = audioread("pink_panther.mp3");
info1 = audioinfo("pink_panther.mp3");
Fs1 = info1.SampleRate;
display(info1);

%display(signal1);

audiowrite("wav_inputs/input_pink_panther.wav",signal1,Fs);

signal1 = audioread("wav_inputs/input_pink_panther.wav");
info = audioinfo("wav_inputs/input_pink_panther.wav");
Fs = info.SampleRate;
ts = info.Duration * Fs;

f = (-ts/2:ts/2-1)*(Fs/ts);

subplot(6,1,2);
plot(f,abs(fftshift(fft(signal1))/length(signal1)));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|X2(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|X2(f)|');
ylim([0,0.012]);
title("Spectrum of input2 signal");

disp("plot")

column1_1 = signal1(:, 1);
column2_1 = signal1(:, 2);

audiowrite("input2_cnls/test1_1.wav",column1_1,Fs);
audiowrite("input2_cnls/test2_1.wav",column2_1,Fs);

a = length(column2);
b = length(column2_1);

%display(length(column2));
%display(length(column2_1));

if(a-b>0)
    max_t = ts1;
    for a = 1 : 1 : a-b
        column2_1(end+1) = 0;
        column1_1(end+1) = 0;

    end
else
    max_t = ts;
    for a = 1 : 1 : b-a
        column2(end+1) = 0;
        column1(end+1) = 0;
    end
end

%display(length(column2));
%display(length(column2_1)); 
%display(length(column1));
%display(length(column1_1));

S = zeros(length(column2),2);
S(:,1)  = column2;
S(:,2)  = column2_1;


rng default % For reproducibility
%mixdata = S*randn(2) + randn(1,2);
mixdata = zeros(length(column2),2);

% Here we are considering the values of a1,b1,a2, and b2(Equations mentioned at the starting of the code) as {2,1,1,2}
% respectively.
mixdata(:,1) = 2*S(:,1) + S(:,2);
mixdata(:,2) = S(:,1) + 2*S(:,2);

f = (-max_t/2:max_t/2-1)*(Fs/max_t);

subplot(6,1,3);
plot(f,abs(fftshift(fft(mixdata(:,1)))/length(mixdata(:,1))));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|M1(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|M1(f)|');
ylim([0,0.012]);
title("Spectrum of mix1 signal");

subplot(6,1,4);
plot(f,abs(fftshift(fft(mixdata(:,2)))/length(mixdata(:,2))));
grid on;
grid minor;
xlim([-10000,10000]);
legend('|M2(f)|');
legend boxoff;
xlabel('f(Hz)');
ylabel('|M2(f)|');
ylim([0,0.012]);
title("Spectrum of mix2 signal");




for i = 1:2
    disp(i);
    %sound(S(:,i));
    %pause;
end

for i = 1:2
    disp(i);
    %sound(mixdata(:,i));
    %pause;
end

audiowrite("ica_mix/mix1.wav",mixdata(:,1),Fs);
audiowrite("ica_mix/mix2.wav",mixdata(:,2),Fs);

% Prewhitening the mixed signals will make the mean of signals zero and
% unity variance. This is very heloful for applying ICA to decrease its
% computation complexity.
mixdata = prewhiten(mixdata);


% For applying ICA we have used a pre-existing function and sperated the
% signals.
Mdl = rica(mixdata,2,'NonGaussianityIndicator',ones(2,1));

% transform fucntion will be helpful for scaling the extracted signals.
unmix = transform(Mdl,mixdata);

audiowrite("unmix/unmix1.wav",unmix(:,1),Fs);
audiowrite("unmix/unmix2.wav",unmix(:,2),Fs);


% This is the prewhiten function used for making mean of the signal as zero and unity variance.
function [X, W] = prewhiten(X)
    X = X - repmat(mean(X, 1), [size(X, 1) 1]);
    W = inv(sqrtm(cov(X)));
    X = X * W; 
end

