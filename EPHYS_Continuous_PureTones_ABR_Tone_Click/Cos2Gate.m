function tone = Cos2Gate(fs,f,phi,dur,ramp)

% generates a sine wave with frequency f, duration dur, phase phi, and 
% sampling frequency fs. The sine wave is then gated on both ends with a 
% squared cosine gate with duration ramp. 

% sine wave without gating
dur = dur/1000;ramp = ramp/1000;
N = round(dur*fs);
t = linspace(0,dur,N);
tone = sin(2*pi*f*t+phi);
% squared cosine gating. 
omega = (1/ramp)*(acos(sqrt(0))-acos(sqrt(1)));
tt = 0:1/fs:pi/2/omega;
tt=tt(1:end-1);
edge = cos(omega*tt).^2;
tone(1:length(edge)) = tone(1:length(edge)) .* fliplr(edge);
tone(end-length(edge)+1 : end) = tone(end-length(edge)+1 : end) .* edge;
