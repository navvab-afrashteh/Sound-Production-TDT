function [FMSweepG, N] = FMchirp(fs,min_Freq,max_Freq,Sdir,Pulse_Dur,ramp,type)
% Generate band-limited noise. Then gating and amplitude modulation in the
% end.
if ~exist('type','var')
    type = 'logarithmic';
end
Pulse_Dur = Pulse_Dur/1000; % convert to sec
N = round(Pulse_Dur*fs);
t = linspace(0,Pulse_Dur,N);

switch Sdir
    case 'up'
        f0 = min_Freq;
        fT = max_Freq;
    case 'down'
        f0 = max_Freq;
        fT = min_Freq;
end
f0 = f0*1000;
fT = fT*1000;
FMSweep = chirp(t,f0,Pulse_Dur,fT,type);

% Signal Gating
omega = (1000/ramp)*(acos(sqrt(0))-acos(sqrt(1)));
tt = 0:1/fs:pi/2/omega;
tt=tt(1:end-1);
edge = cos(omega*tt).^2;
Nm = N - 2*length(edge);
gate = [fliplr(edge), ones(1,Nm),edge];
FMSweepG = FMSweep.*gate;