function [Gate] = Cos2Gate(Fs,ramp,unit)
% Fs: sampling frequency
if ~exist('unit','var')
    unit = 'ms';
end
if strcmpi(unit(1),'m')
    ramp = ramp/1000;
end
% ramp (gating) duration in sec
if ramp < 1/Fs
    Gate = 1;
    return;
end
omega = (1/ramp)*(acos(sqrt(0))-acos(sqrt(1))); % we go from 0 to 1 in ramp msec
N = round(Fs*ramp);
t = linspace(0,pi/2/omega,N);
Gate = cos(omega*t).^2;
Gate = fliplr(Gate); % gate of onset
