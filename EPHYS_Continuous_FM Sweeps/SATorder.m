function [SweepTimes, SweepTimeOrder, SweepDirOrder, rampTimeOrder, AmpOrder, TimeOrder] = ...
    SATorder(Pulse_Dur_spec, ramp_spec, Amp_Array,...
	ISI, jitter, Ntrials, randomized, UpDown, Log_Levels, varargin)

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
if nargin > 8
    FsCam = varargin{1};
    TCam = 1/FsCam;
end
nFrames = nan;
if nargin > 9
    nFrames = varargin{2};
end
tot_recDur = nFrames/FsCam;
NCam = 5;
if nargin > 9
    NCam = varargin{3};
end
SweepTimes = octaveSweepTimes(Pulse_Dur_spec, Log_Levels);
rampTimes = octaveSweepTimes(ramp_spec, Log_Levels);

nST = length(SweepTimes);
nA = length(Amp_Array);
nSD = UpDown+1; % if UpDown=1 then there should be up-sweeps and down-sweeps
if UpDown
    SweepDir = {'up'; 'down'};
else
    SweepDir = {'up'};
end
ISI = ISI/1000;
jitter = jitter/1000;

SweepTimeOrder = zeros(Ntrials*nST*nA*nSD,1);
SweepDirOrder = cell(Ntrials*nST*nA*nSD,1);
rampTimeOrder = zeros(Ntrials*nST*nA*nSD,1);
AmpOrder = zeros(Ntrials*nST*nA*nSD,1);
TimeOrder  = zeros(Ntrials*nST*nA*nSD,1);

NtotCam = Ntrials*nST*nA*nSD/NCam;
NperCam = Ntrials/NCam;
STorder = zeros(NtotCam,1);
rTorder = zeros(NtotCam,1);
Aorder = zeros(NtotCam,1);
SDorder = cell(NtotCam,1);
kk =  0;
for n = 1:NperCam
    for t = 1:nST
        for a = 1:nA
            for sd = 1:nSD
                kk = kk+1;
                STorder(kk) = SweepTimes(t);
                rTorder(kk) = rampTimes(t);
                Aorder(kk) = Amp_Array(a);
                SDorder{kk} = SweepDir{sd};
            end
        end
    end
end
if randomized
    for ncam = 1:NCam
        randOrder = randperm(NtotCam);
        STtemp = STorder(randOrder);
        rTtemp = rTorder(randOrder);
        Atemp = Aorder(randOrder);
        SDtemp = SDorder(randOrder);
        SweepTimeOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = STtemp;
        rampTimeOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = rTtemp;
        AmpOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = Atemp;
        SweepDirOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = SDtemp;
    end
else
    for ncam = 1:NCam
        SweepTimeOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = STorder;
        rampTimeOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = rTorder;
        AmpOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = Aorder;
        SweepDirOrder((ncam-1)*NtotCam+1:ncam*NtotCam) = SDorder;
    end
end

for p = 1:Ntrials*nST*nA*nSD
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
