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

% n1 = nextpow2(f1)-1;
% n2 = nextpow2(f2);
% Foctave = [];
% 
% for k = n1:n2
%     f = linspace(k,k+1,fpoints+1);
%     Foctave = [Foctave, f(1:end-1)];
% end
% 
% p = 1;
% F = [];
% for k=1:length(Foctave)
%     if 2^Foctave(k) >= f1 && 2^Foctave(k) <= f2 
%         F(p) = 2.^Foctave(k);
%         p = p+1;
%     end
% end

F = [];
f = f1;
while f <= f2+0.001
	F = [F, f];
	f = f*2^(1/fpoints);
end
F = F.';