function [F, FreqOrder, AmpOrder] = FATorder(Fspec, Amp_Array, ISI, jitter, Ntrials, randomized,varargin)

% Fspec : vector of frequencies specifications: [first freq, last freq, frequency per octave],in kHz
% Amp_Array : vector of amplitudes in dB
% ramp : ramp time for the start and the end of sine wave, in milli second
% dur : duration of tones, in milli second
% ISI : time between two consecutive tones, in milli second
% jitter :  jitter for silence period, in mili second
% Ntrials : number of all complex sounds, integer number
% randomized : indicating if the signal be generated in a randomized manner
% or not. 1 or 0
% Edited by Navvab Afrashteh, May 22, 2017

FsCam = nan;
if nargin > 6
    FsCam = varargin{1};
    TCam = 1/FsCam;
end
nFrames = nan;
if nargin > 7
    nFrames = varargin{2};
end
tot_recDur = nFrames/FsCam;
NCam = 5;
if nargin > 8
    NCam = varargin{3};
end
F = octavefreq(Fspec);
nF = length(F);
nA = length(Amp_Array);

ISI = ISI/1000;
jitter = jitter/1000;

FreqOrder = zeros(Ntrials*nF*nA,1);
AmpOrder = zeros(Ntrials*nF*nA,1);
TimeOrder  = zeros(Ntrials*nF*nA,1);

NtotCam = Ntrials*nF*nA/NCam;
NperCam = Ntrials/NCam;
Forder = zeros(NtotCam,1);
Aorder = zeros(NtotCam,1);
kk =  0;
for n = 1:NperCam
    for f = 1:nF
        for a = 1:nA
            kk = kk+1;
            Forder(kk) = F(f);
            Aorder(kk) = Amp_Array(a);
        end
    end
end
if randomized
    for ncam = 1:NCam
        randOrder = randperm(NtotCam);
        Ftemp = Forder(randOrder);
        Atemp = Aorder(randOrder);
        FreqOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = Ftemp;
        AmpOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = Atemp;
    end
else
    for ncam = 1:NCam
        FreqOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = Forder;
        AmpOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = Aorder;
    end
end

for p = 1:Ntrials*nF*nA
    Tsilence = ISI + (2*rand-1)*jitter;
    if ~isnan(FsCam)
        Tsilence = TCam*round(Tsilence/TCam);
    end
    TimeOrder(p) = Tsilence;
end

if ~isnan(tot_recDur)
    for trial = 1:NCam
        idx = 1+(trial-1)*NtotCam : trial*NtotCam;
        timeOrderTemp = TimeOrder(idx);
        s = sum(timeOrderTemp);
        if s > tot_recDur-ISI % sum of silence time be at least ISI sec less than total recording time
            d = (s+ISI-tot_recDur)/length(timeOrderTemp);
            timeOrderTemp = timeOrderTemp-d;
            TimeOrder(idx) = timeOrderTemp;
        end
    end
end
