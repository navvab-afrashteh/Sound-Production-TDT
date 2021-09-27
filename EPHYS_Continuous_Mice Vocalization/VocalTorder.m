function [VocalOrder, TimeOrder] = VocalTorder(Nvocal, NCam, ISI, jitter, Ntrials, randomized,varargin)

% Nvocal : number of vocalization files
% ISI : time between two consecutive tones, in milli second
% jitter :  jitter for silence period, in mili second
% Ntrials : number of all complex sounds, integer number
% randomized : indicating if the signal be generated in a randomized manner
% or not. 1 or 0
% Edited by Navvab Afrashteh, Sept 26, 2017

FsCam = nan;
if nargin > 5
    FsCam = varargin{1};
    TCam = 1/FsCam;
end
nFrames = nan;
if nargin > 6
    nFrames = varargin{2};
end
tot_recDur = nFrames/FsCam;

ISI = ISI/1000;
jitter = jitter/1000;

TimeOrder  = zeros(Ntrials*Nvocal,1);

p = 1;
for trial = 1:Ntrials*Nvocal
    Tsilence = ISI + (2*rand-1)*jitter;
    if ~isnan(FsCam)
        Tsilence = TCam*round(Tsilence/TCam);
    end
    TimeOrder(p) = Tsilence;
    p = p+1;
end

NtotCam = Nvocal*Ntrials/NCam;
VocalOrder = zeros(Ntrials*Nvocal,1);
if randomized
    for n = 1:NCam
        idx = 1+(n-1)*NtotCam : n*NtotCam;
        VocalOrder(idx) = mod(randperm(Ntrials*Nvocal/NCam),Nvocal)+1;
    end
end

if ~randomized
    for n = 1:Ntrials
        idx = 1+(n-1)*Nvocal : n*Nvocal;
        VocalOrder(idx) = 1:Nvocal;
    end
end

ISIend = max(ISI,4);
if ~isnan(tot_recDur)
    for trial = 1:NCam
        idx = 1+(trial-1)*NtotCam : trial*NtotCam;
        timeOrderTemp = TimeOrder(idx);
        s = sum(timeOrderTemp);
        if s > tot_recDur-ISIend % sum of silence time be at least ISI sec less than total recording time
            d = (s+ISIend-tot_recDur)/length(timeOrderTemp);
            timeOrderTemp = timeOrderTemp-d;
            TimeOrder(idx) = timeOrderTemp;
        end
    end
end

n=0;