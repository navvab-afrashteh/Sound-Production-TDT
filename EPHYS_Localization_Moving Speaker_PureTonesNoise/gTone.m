function tone = gTone(fs,f,dur,ramp)
% dur and ramp in msec
dur = dur/1000;
ramp = ramp/1000;

t = linspace(0, dur, fs*dur); % generate time stamps
tone = sin(2*pi*f*t);
% generate ramp
omega = (1/ramp)*(acos(sqrt(0))-acos(sqrt(1))); % 0% to 100% in ramp duration
tt = 0:1/fs:pi/2/omega;
tt = tt(1:end-1);
edge = cos(omega*tt).^2;
tone(1:length(edge)) = tone(1:length(edge)) .* fliplr(edge);
tone(end-length(edge)+1 : end) = tone(end-length(edge)+1 : end) .* edge;