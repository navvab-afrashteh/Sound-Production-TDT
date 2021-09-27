function [noiseG, N] = WideBandNoise(fs,Pulse_Dur,ramp)
% Generate band-limited noise. Then gating and amplitude modulation in the
% end.
% Pulse_Dur and ramp in msec
Pulse_Dur = Pulse_Dur/1000;
ramp = ramp/1000;
% generate broadband noise
N = round(fs*Pulse_Dur); % generate the pulse duration length
noise = randn(1,N);
% Signal Gating
omega = (1/ramp)*(acos(sqrt(0))-acos(sqrt(1))); % 0% to 100% in ramp duration
tt = 0:1/fs:pi/2/omega;
tt = tt(1:end-1);
edge = cos(omega*tt).^2;
Nm = N - 2*length(edge);
gate = [fliplr(edge), ones(1,Nm),edge];
noiseG = noise.*gate;