close all
clear all
clc
% Hardware parameters
Fs   = 97656.25*2;     % TDT sampling frequency in Hz
FsCam        = 30;     % Camera sampling frequency in Hz
SpeakerName  = 'ES1-A'; % speaker name
% Stimulus parameters
nFrames = 7000; % number of frames
minOnset = 1; % minimum silence at begining in s
minOffset = 2; % minimum silence in the end in s
recDur = nFrames/FsCam; % recording duration in s
saveFlag     = 1;      % save?
%-------------------------------------------------------------------------------------------------------------------------------------
%% Vocalization parameters
amp = 70;
muISI = 0.5;
seed = 1;
VocalStimForSTRF(Fs,FsCam,amp,muISI,minOnset,minOffset,recDur,saveFlag,seed,SpeakerName)
close all
%-------------------------------------------------------------------------------------------------------------------------------------
%% DRC parameters
fmin = 2000; fmax = 96000; FreqPerOct = 8; Amps = [25:10:75];%Amps=70;
dur = [0.033, 0.066];%dur=0.05;
ramp = [0.002, 0.005];%ramp=0.000005;
muISI = 0.5;
seed = 1;
ToneStimForSTRF(Fs,FsCam,fmin,fmax,FreqPerOct,Amps,dur,ramp,muISI,minOnset,minOffset,recDur,saveFlag,seed,SpeakerName)
close all
%-------------------------------------------------------------------------------------------------------------------------------------
%% DMR parameters
MaxFreq      = 96000;  % max frequency in Hz
MinFreq      = 2000;    % min frequency in Hz
ramp         = 5;      % ramp (gating) duration in msec
FsFm         = 1.5;    % max time-varying temporal modulation change-rate
FsRd         = 3;      % max instantaneous ripple density change-rate
FreqPerOct   = 50;     % number of carrier freqq. per octave, Escabi and Schreiner used ~43
MaxFm        = 40;     % max. time-varying temporal modulation rate
MaxRd        = 4;      % max instantaneous ripple density rate
M            = 40;     % Modulation depth of envelope in dB
Amp          = 70;     % Mean sound pressure level
seed         = 1;      % seed value for reproducability
DMRStimForSTRF(Fs,FsCam,recDur,minOnset,minOffset,MaxFreq,MinFreq,ramp,FsFm,FsRd,FreqPerOct,MaxFm,MaxRd,M,saveFlag,seed,SpeakerName,Amp);
close all
%-------------------------------------------------------------------------------------------------------------------------------------
%% FM Sweeps parameters
fmin = 2000; fmax = 96000; type = 'logarithmic'; NB = 20; Amps = [25:10:75]; %Amps=70;
dur = [0.01, 0.1];%dur=0.033;
ramp = [0.002, 0.005];%ramp=0.002;
muISI = 0.5;
seed = 1;
FMStimForSTRF(Fs,FsCam,fmin,fmax,type,NB,Amps,dur,ramp,muISI,minOnset,minOffset,recDur,saveFlag,seed,SpeakerName)
close all
