function [times, ToneRep]= genTimes(alldurs,allramps,Fs,FsCam,recDur,minOnset,minOffset,muISI,seed)
% generate exponentially distributed random times for onsets of stimuli

% Fs = 200000;
% nFrames = 7000;
% FsCam = 30;
% minOnset = 1;
% minOffset = 2;
% muISI = 0.5;
% recDur = nFrames/FsCam;
if ~exist('seed','var')
    seed = round(rand*2^20);
end
rng(seed);

minN = round(minOnset*Fs);
maxN = round((recDur-minOffset)*Fs);
muN = round(muISI*Fs);
avgNT = ceil(recDur/muISI);
NTones = size(alldurs,1);
times = cell(NTones,1);
ToneRep = zeros(NTones,1);
for tone = 1:NTones
    maxDur = alldurs(tone,end);
    maxRamp = allramps(tone,end);
    minISI = round(maxDur*Fs);
    dn = exprnd(muN,2*avgNT,1);
    dn = round(dn);
    dn(dn<minISI) = [];
    d = cumsum(dn);
    d(d<minN) = [];
    d(d>maxN) = [];
    t = d/Fs;
%     t = floor(t*FsCam)/FsCam; % make sure that times are multipliers of 1/FsCam
%     d = round(t*Fs); % update indices
    times{tone}.inds = d;
    times{tone}.times = t;
    if size(alldurs,2)>1
        minDur = alldurs(tone,1);
        times{tone}.dur = rand(length(d),1)*(maxDur-minDur) + minDur;
        if size(allramps,2)>1
            minRamp = allramps(tone,1);
            m = (maxRamp-minRamp)/(maxDur-minDur);
            times{tone}.ramp = (times{tone}.dur-minDur)*m+minRamp;
        else
            times{tone}.ramp = maxRamp*ones(length(d),1);
        end
    else
        times{tone}.dur = maxDur*ones(length(d),1);
        times{tone}.ramp = maxRamp*ones(length(d),1);
    end
    times{tone}.N = round(times{tone}.dur*Fs);
    ToneRep(tone) = length(d);
end

