function Amp = genAmps(Amps,ToneRep,seed)

if ~exist('seed','var')
    seed = round(rand*2^20);
end
rng(seed)
NTones = length(ToneRep);
nAmps = length(Amps);
Amp = cell(NTones,1);
for tone = 1:NTones
    idx = randi(nAmps,ToneRep(tone),1);
    Amp{tone} = Amps(idx);
end