clc
clear

fs = 195313;
T = 0.1; % 100msec
N = floor(fs*T);
t = linspace(0,T,N);

% prepare ramp
Tr = 0.005; % 5msec
Nr = floor(fs * Tr);
r = sin(linspace(0, pi/2, Nr));
r = [r, ones(1, N - Nr * 2), fliplr(r)];

freqs = 2:11; %kHz
PathName = sprintf('%s\\Sine Waves',pwd);
for f = freqs
    x = sin(2*pi*f*1000*t);
    xr = x.*r;
    filePath = sprintf('%s\\sineWave_%dkHz.wav',PathName,f);
    audiowrite(filePath,xr,fs)
end


 