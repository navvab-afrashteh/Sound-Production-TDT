function [TimeOrder, FreqOrder, idxMMN, FmainAB] = blockAB(FI,FII,minDist, MMPercent, ...
    recDur, PulseDur, ISI, jitter,NperCam,seed)

% fmain: dominant frequency
% fMM : mismatch frequency
% minSpace: minimum number of fmain in between two fMM
% MMFraction: fraction of all stims alocated to mismatch tone
% recDur, stimDur, ISI and jitter: total recording and stim durations,
% inter-stim-interval and maximum jitter time
% Nstim: number of all stims (fmain+fMM)

rng(1000*(seed));

f = [FI,FII];
idx = randi(2);
F1 = f(idx); f(idx) = [];
F2= f;
FmainAB = [F1,F2];
FMM = [F2,F1];

jitter = ISI*jitter;
NMM = round(NperCam*MMPercent/100);

idxMMN = nan(NMM,2);
FreqOrder = nan(NperCam,2);
TimeOrder = nan(NperCam,2);
for k = 1:2
    idxMMN(:,k) = randiMinDist(minDist+1,NperCam,minDist,NMM);
    idxMMN(:,k) = sort(idxMMN(:,k));
    FreqOrder(:,k) = FmainAB(k)*ones(NperCam,1);
    FreqOrder(idxMMN(:,k),k) = FMM(k);
    
    TimeOrder(:,k) = rand(NperCam,1)*jitter + ISI + PulseDur;
    d = recDur - sum(TimeOrder(:,k)) - ISI - PulseDur;
    d = d/NperCam;
    if d<0
        TimeOrder(:,k) = TimeOrder(:,k) + d;
    end
end

