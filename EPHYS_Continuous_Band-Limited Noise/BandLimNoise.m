function [BLnoiseGAM, N, f1, f2] = BandLimNoise(fs,f,BW_Oct,Pulse_Dur,ramp,Mod_Freq)
% Generate band-limited noise. Then gating and amplitude modulation in the
% end.

% get frequency range
f1 = f/(2^(BW_Oct/2));
f2 = f*(2^(BW_Oct/2));
% generate broadband noise
N1 = round(2*fs*Pulse_Dur/1000); % generate twice the pulse duration length
noise = randn(1,N1);
% filter design
filtDur = 2*10^-3; % 2msec for filter duration
filtN = round(fs*filtDur);
Hd = BPF_TDT(fs, filtN, f1, f2);
H = Hd.Numerator;
% filter the noise
BLnoise = filtfilt(H,1,noise);
Nmid = round(N1/2);
N = round(fs*Pulse_Dur/1000);
Nhalf = round(N/2);
BLnoise = BLnoise(Nmid-Nhalf:Nmid+Nhalf-1);
N = length(BLnoise);
% Signal Gating
omega = (1000/ramp)*(acos(sqrt(0))-acos(sqrt(1)));
tt = 0:1/fs:pi/2/omega;
tt=tt(1:end-1);
edge = cos(omega*tt).^2;
Nm = N - 2*length(edge);
gate = [fliplr(edge), ones(1,Nm),edge];
BLnoiseG = BLnoise.*gate;
% amplitude modulation
t = 0:1/fs:(N-1)/fs;
BLnoiseGAM = (1+sin(2*pi*Mod_Freq*t-pi/2)).*BLnoiseG;