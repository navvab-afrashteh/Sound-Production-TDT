function [TimeOrder, FreqOrder] = blockC(F, recDur, stimDur, ISI, jitter,Nstim,seed)
% F: freq list
% recDur, stimDur, ISI and jitter: total recording and stim durations,
% inter-stim-interval and maximum jitter time
% Nstim: number of all stims (for all freq list)
rng(2000*(seed));

nF = length(F);
NperStim = Nstim/nF;
FreqOrder = repmat(F(:)',NperStim,1);
pF = randperm(NperStim*nF);
FreqOrder = FreqOrder(pF);
FreqOrder = FreqOrder(:);
jitter = ISI*jitter;
TimeOrder = rand(Nstim,1)*jitter + ISI + stimDur;
d = recDur - sum(TimeOrder) - ISI - stimDur;
d = d/Nstim;
if d<0
    TimeOrder = TimeOrder + d;
end