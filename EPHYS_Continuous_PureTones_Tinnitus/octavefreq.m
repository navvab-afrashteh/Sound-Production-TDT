function F = octavefreq(Fspec)
% Fspec = [f1:f2,fpoints] or Fspec = [f1,f2,fpoints]
% Fspec = [1:82,3] or Fspec = [1,82,3]

f1 = Fspec(1); % the minimum frequency (linear scale, kHz)
f2 = Fspec(end-1); % the maxmimum frequency (linear scale, kHz)
fpoints = Fspec(end); % frequency per octave value 

if f1 < 1
    f1 = 1;
end
if f2 > 100
    f2 = 100;
end

F = [];
f = f1;
while f <= f2+0.001
	F = [F, f];
	f = f*2^(1/fpoints);
end
F = F.';

% n = nextpow2(f2/f1)+1;
% Foctave = 0:1/fpoints:n;
% F = f1*(2.^Foctave);
% F(F<f1) = [];
% F(F>f2) = [];
