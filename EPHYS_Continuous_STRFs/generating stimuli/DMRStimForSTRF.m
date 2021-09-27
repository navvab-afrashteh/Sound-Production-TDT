%{
 **************   Script_Generate_DMRs   ******************
Modified from bellow by Navvab Afrashteh.
Date  : 22-Aug-2018 

This script generates Dynamic Moving Ripple (DMR) stimuli with randomly
generated parameters. DMRs are an extension of moving Ripple by

  Monty A. Escabi and Christoph E. Schreiner: Nonlinear Spectrotemporal
  Sound Analysis by Neurons in the Auditory Midbrain. The Journal of
  Neuroscience, May 15, 2002, 22(10):4114-4131.

  The basic formula is

     S(x,t) = M/2 * sin( 2 * pi(omega of t + Omega of x) + Phase phi of t)

  where M is the modulation depth of the envelope in decibels (30dB or 45dB)
  and x = log2 (F/F0) where F0 is the basefreq (500 Hz) and F is the
  frequency of the tone.

  The phase phi is a function of time and controls the time-varying
  temporal modulation rate ( Phase Phi of t = Integral from 0 till t of
  Fm of t as a function of t where Fm is the time-varying temporal
  modulation rate).

  Fm and Omega are slowly time variying signals (Fm <= 1.5 Hz and Omega <=
  3 Hz)

See also Script_Generate_FMBanks

Author: arne.f.meyer@uni-oldenburg.de
Date  : 22-Nov-2010 09:41:44
%}

function DMRStimForSTRF(Fs,FsCam,recDur,minOnset,minOffset,MaxFreq,MinFreq,ramp,FsFm,FsRd,FreqPerOct,MaxFm,MaxRd,M,saveFlag,seed,SpeakerName,Amp)
%% -----------------------------------------------
% Script settings
% -----------------------------------------------
if ~exist('seed','var')
    seed = round(rand*2^20);
end
seedOrig = seed;
sDataDir   = fullfile(fileparts(mfilename('fullpath')), 'DMRs');
if saveFlag && ~exist(sDataDir, 'dir')
    mkdir(sDataDir);
end
sDataDir   = sprintf('%s\\seed#%d',sDataDir,seed);
if saveFlag && ~exist(sDataDir, 'dir')
    mkdir(sDataDir);
end
StimLen  = recDur-minOnset-minOffset;    % stimulus length is s
nOct     = log2(MaxFreq/MinFreq);        % number of octaves covered
% Define some variables from basic parameters
Time    = 1/Fs   : 1/Fs : StimLen;       % time vector
TimeCam = 1/FsCam   : 1/FsCam : StimLen; % camera time vector
TimeFm  = 0 : 1/(2*FsFm) : StimLen;      % vector for sampling of temporal modulation rate
TimeRD  = 0 : 1/(2*FsRd) : StimLen;      % vector for sampling of inst. ripple density
% Vector containing carrier frequencies
Freqs = MinFreq * 2.^linspace(0, nOct, FreqPerOct*nOct)';
% Use an uniform distribution for frequency modulation
rng(seed); seed = seed + seed*1024 + 1;
Fm        = rand(size(TimeFm));
FmUpSamp  = interp1(TimeFm, Fm, Time, 'pchip');
FmUpSampCam  = interp1(TimeFm, Fm, TimeCam, 'pchip');
Fm        = (FmUpSamp.* MaxFm*2) - MaxFm;     % take Fm from [-MaxFm, MaxFm] Hz interval
FmCam     = (FmUpSampCam.* MaxFm*2) - MaxFm;  % take Fm from [-MaxFm, MaxFm] Hz interval
% Computation of Phase of t (integral of temporal modulation rate)
PhaseFm = 2*pi*cumsum(Fm)/Fs;
PhaseFmCam = 2*pi*cumsum(FmCam)/FsCam;
% Use an uniform distribution for ripple density
rng(seed); seed = seed + seed*1024 + 1;
RD        = rand(size(TimeRD));
RDUpSamp  = interp1(TimeRD,RD,Time,'pchip');
RDUpSampCam  = interp1(TimeRD,RD,TimeCam,'pchip');
RD        = (RDUpSamp.* MaxRd);     % take RD from [0, MaxRd] cyc/octaves interval
RDCam     = (RDUpSampCam.* MaxRd);  % take RD from [0, MaxRd] cyc/octaves interval
% generate a random phase for each frequency
rng(seed); seed = seed + seed*1024 + 1;
nFreqs = size(Freqs,1);
Phases    = rand([1,nFreqs])*2*pi;
LogFreq   = log2(Freqs/MinFreq);
%% -----------------------------------------------
% Generate waveform
% -----------------------------------------------
% Compute the stimulus by multiplying carriers with the (exponentialized)
% envelope.
if exist('SpeakerName','var') && exist('Amp','var')
    AmpOrder = Amp*ones(nFreqs,1);
    AmpTDT = getTDTAmp(Freqs/1000, Freqs/1000, AmpOrder, SpeakerName);
else
    AmpTDT = ones(nFreqs,1);
end
StimTDT = zeros(size(Time)); % output sound sequence to be played by TDT
Stim = zeros(size(Time));    % output sound sequence
for f = 1:nFreqs
    StimTDT = StimTDT + AmpTDT(f) * sin( 2*pi*Freqs(f) .* Time + Phases(f) ) .* ...
        10.^( ((M/2) * sin( 2*pi*LogFreq(f)*RD + PhaseFm) - M/2) /20 );
    Stim = Stim + sin( 2*pi*Freqs(f) .* Time + Phases(f) ) .* ...
        10.^( ((M/2) * sin( 2*pi*LogFreq(f)*RD + PhaseFm) - M/2) /20 );
end
Gate = Cos2Gate(Fs,ramp); % onset and offset gate sequences
StimTDT = applyGate(StimTDT,Gate); % apply gate sequence to onset and offset of the stimulus
% adjust amplitudes for required dB SPL
if exist('SpeakerName','var') && exist('Amp','var')
    [calib_V, calib_level] = getCalib(SpeakerName);
    [StimTDT, dTDT] = adjAmpSPL(StimTDT,Amp,calib_level,calib_V);
    [Stim, d] = adjAmpSPL(Stim,Amp,calib_level,calib_V);
end
% spectrogram data (in dB)
Stx = -inf(nFreqs,length(TimeCam)); % actual spectrogram data (downsampled)
for f = 1:nFreqs
    Stx(f,:) = (M/2) *sin( 2*pi*LogFreq(f)*RDCam + PhaseFmCam)- M/2 +d + Amp; % d is added to be comparable to othe stim
end
% add silence to start and end
N1 = round(minOnset*Fs);
N2 = round(minOffset*Fs);
StimTDT = [zeros(1,N1), StimTDT, zeros(1,N2)];
Stim = [zeros(1,N1), Stim, zeros(1,N2)];
N1Cam = round(minOnset*FsCam);
N2Cam = round(minOffset*FsCam);
Stx = [-inf(nFreqs,N1Cam), Stx, -inf(nFreqs,N2Cam)];
% update time stamps
TimeCam   = (1:size(Stx,2))/FsCam;    % camera time vector
%% -----------------------------------------------
% Plotting
% -----------------------------------------------
% plot actual stimulus spectrogram
figure;
imagesc(TimeCam, Freqs, Stx,[d-M, d]+Amp)
axis xy; colorbar;
xlabel('Time / s', 'FontSize', 16), xlim([0 recDur])
ylabel('Frequency / Hz', 'FontSize', 16), ylim([MinFreq MaxFreq])
set(gca(), 'FontSize', 14)
if saveFlag
%     saveas(gcf(), fullfile(sDataDir, 'Spectrogram_Actual'), 'pdf');
    figName = 'Spectrogram_Actual';
    dpi = 600;
    pdfFilePath = sprintf('%s\\%s.pdf',sDataDir,figName);
    save2pdf(pdfFilePath,gcf,dpi);
end
% plot estimated stimulus spectrogram
if ~isempty(SpeakerName)
    [calib_V, calib_level, calib_freq, calib_amp,FIRcoeff] = getCalib(SpeakerName,Fs);
else
    calib_level = 14.75;
    calib_V = nan; calib_freq = nan; calib_amp = nan; FIRcoeff = nan;
end
nWinLen               = 0.0025;
nOverlap              = nWinLen / 2;
[spS, spFreq, spTime,spP] = spectrogram(StimTDT, kaiser(ceil(nWinLen * Fs),5), ...
    round(nOverlap * Fs), 2^11, Fs);
spPdB = 10*log10(abs(spP));
cPdB = interp1(calib_freq,calib_amp,spFreq,'pchip');
cPdB(spFreq>calib_freq(end)) = calib_amp(end);
cPdB(spFreq<calib_freq(1)) = calib_amp(1);
cPdB = repmat(cPdB,1,length(spTime));
spPdB = spPdB-cPdB+calib_level-10*log10(nWinLen)-14.75;

nFsSpec               = 1 / (spTime(2) - spTime(1));
figure('Name','Stimulus spectrogram');
imagesc(spTime, spFreq, spPdB,[d-M, d]+Amp)
axis xy; colorbar;
xlabel('Time / s', 'FontSize', 16), xlim([0 recDur])
ylabel('Frequency / Hz', 'FontSize', 16), ylim([MinFreq MaxFreq])
set(gca(), 'FontSize', 14)
if saveFlag
%     saveas(gcf(), fullfile(sDataDir, 'Spectrogram_Estimated'), 'pdf');
    figName = 'Spectrogram_Estimated';
    dpi = 600;
    pdfFilePath = sprintf('%s\\%s.pdf',sDataDir,figName);
%     save2pdf(pdfFilePath,gcf,dpi);
end
% Actual Frequency-wise autocorrelation function
StxLin = 10.^(Stx/20);
MaxLagSec = 0.5;
nMaxLag    = ceil(MaxLagSec * FsCam);
Corr      = zeros(size(StxLin,1), 2*nMaxLag+1);
tLags = (-nMaxLag:nMaxLag) / FsCam;
for iFreq=1:size(StxLin,1)
    Corr(iFreq,:) = xcorr(abs(StxLin(iFreq,:)), abs(StxLin(iFreq,:)), nMaxLag);
    Corr(iFreq,:) = Corr(iFreq,:) / max(Corr(iFreq,:));
end
figure('Name','Frequency-wise stimulus autocorrelation');
imagesc(tLags, spFreq, Corr,[0,1]);
axis xy; colorbar
xlabel('Time Lag (s)', 'FontSize', 16), xlim(MaxLagSec * [-1 1])
ylabel('Frequency (Hz)', 'FontSize', 16), ylim([MinFreq MaxFreq])
set(gca(), 'FontSize', 14)
if saveFlag
%     saveas(gcf(), fullfile(sDataDir, 'AutoCorrelation_Actual'), 'pdf');
    figName = 'AutoCorrelation_Actual';
    dpi = 600;
    pdfFilePath = sprintf('%s\\%s.pdf',sDataDir,figName);
    save2pdf(pdfFilePath,gcf,dpi);
end
ActualCorr.Corr = Corr;
ActualCorr.MaxLagSec = MaxLagSec;
ActualCorr.tLags = tLags;
ActualCorr.Freqs = spFreq;
% Estimated Frequency-wise autocorrelation function
nMaxLag    = ceil(MaxLagSec * nFsSpec);
Corr      = zeros(size(spP,1), 2*nMaxLag+1);
tLags = (-nMaxLag:nMaxLag) / nFsSpec;
for iFreq=1:size(spP,1)
    Corr(iFreq,:) = xcorr(abs(spP(iFreq,:)), abs(spP(iFreq,:)), nMaxLag);
    Corr(iFreq,:) = Corr(iFreq,:) / max(Corr(iFreq,:));
end
figure('Name','Frequency-wise stimulus autocorrelation');
imagesc(tLags, spFreq, Corr,[0,1]);
axis xy; colorbar
xlabel('Time Lag (s)', 'FontSize', 16), xlim(MaxLagSec * [-1 1])
ylabel('Frequency (Hz)', 'FontSize', 16), ylim([MinFreq MaxFreq])
set(gca(), 'FontSize', 14)
if saveFlag
%     saveas(gcf(), fullfile(sDataDir, 'AutoCorrelation_Estimated'), 'pdf');
    figName = 'AutoCorrelation_Estimated';
    dpi = 600;
    pdfFilePath = sprintf('%s\\%s.pdf',sDataDir,figName);
    save2pdf(pdfFilePath,gcf,dpi);
end
EstCorr.Corr = Corr;
EstCorr.MaxLagSec = MaxLagSec;
EstCorr.tLags = tLags;
EstCorr.Freqs = spFreq;
%% -----------------------------------------------
% Save data
% -----------------------------------------------
if ~exist('SpeakerName','var')
    SpeakerName = NaN;
end
if ~exist('Amp','var')
    Amp = NaN;
end
spec.recDur       = recDur;      % recording length is s
spec.minOnset     = minOnset;    % begining silence length is s
spec.minOffset    = minOffset;   % end silence length is s
spec.StimLen      = StimLen;     % stimulus length is s
spec.Fs           = Fs;          % TDT sampling frequency in Hz
spec.MaxFreq      = MaxFreq;     % max frequency in Hz
spec.MinFreq      = MinFreq;     % min frequency in Hz
spec.ramp         = ramp;        % ramp (gating) duration in msec
spec.FsFm         = FsFm;        % max time-varying temporal modulation change-rate
spec.FsRd         = FsRd;        % max instantaneous ripple density change-rate
spec.FreqPerOct   = FreqPerOct;  % number of carrier freqq. per octave, Escabi and Schreiner used 230 carriers over 5.32 octaves (0.5 - 20kHz)
spec.MaxFm        = MaxFm;       % max. time-varying temporal modulation rate
spec.MaxRd        = MaxRd;       % max instantaneous ripple density rate
spec.M            = M;           % Modulation depth of envelope in dB
spec.FsCam        = FsCam;       % Camera sampling frequency in Hz
spec.SpeakerName  = SpeakerName; % speaker name
spec.Amp          = Amp;         % Mean sound pressure level
spec.AmpTDT       = AmpTDT;
spec.seed         = seedOrig;    % seed value for reproducability
spec.Freqs        = Freqs;
spec.TimeCam      = TimeCam;
% speaker info
Speaker.SpeakerName = SpeakerName;
Speaker.calib_V     = calib_V;
Speaker.calib_level = calib_level;
Speaker.calib_freq  = calib_freq;
Speaker.calib_amp   = calib_amp;
Speaker.FIRcoeff    = FIRcoeff;
% spectrogram outputs
Spectrogram.spS    = spS;
Spectrogram.spFreq = spFreq;
Spectrogram.spTime = spTime;
Spectrogram.spP    = spP;
Spectrogram.spPdB  = spPdB;

matFile = sprintf('%s\\Stim.mat',sDataDir);
if saveFlag
    save(matFile,'Stim','StimTDT','Stx','spec','Speaker','Spectrogram','ActualCorr','EstCorr','-v7.3');
%     audiowrite(sprintf('%s\\Stim.wav', sDataDir), StimTDT, round(Fs));
end

