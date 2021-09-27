function FMStimForSTRF(Fs,FsCam,MinFreq,MaxFreq,SweepType,NBlock,Amps,dur,ramp,muISI,minOnset,minOffset,recDur,saveFlag,seed,SpeakerName)
if ~exist('seed','var')
    seed = round(rand*2^20);
end
seedOrig = seed;
sDataDir   = fullfile(fileparts(mfilename('fullpath')), 'FM Sweeps');
if saveFlag && ~exist(sDataDir, 'dir')
    mkdir(sDataDir);
end
sDataDir   = sprintf('%s\\seed#%d',sDataDir,seed);
if saveFlag && ~exist(sDataDir, 'dir')
    mkdir(sDataDir);
end
if ~exist('SpeakerName','var')
    SpeakerName = '';
end
%%%%%%%%%%%%%
maxDur = dur(end);
maxRamp = ramp(end);
alldurs = maxDur*ones(NBlock,1);
allramps = maxRamp*ones(NBlock,1);
if length(dur)>1
    minDur = dur(1);
    alldurs = [minDur*ones(NBlock,1), alldurs];
end
if length(ramp)>1
    minRamp = ramp(1);
    allramps = [minRamp*ones(NBlock,1), allramps];
end
% generate onset times
[times, ToneRep] = genTimes(alldurs,allramps,Fs,FsCam,recDur,minOnset,minOffset,muISI,seed);
% generate amplitude levels
Amp = genAmps(Amps,ToneRep,seed);
%%
if ~isempty(SpeakerName)
    [calib_V, calib_level, ~, ~, FIRcoeff] = getCalib(SpeakerName,Fs);
end
recDurN = round((recDur)*Fs);
Stim = zeros(recDurN,1);
StimTDT = zeros(recDurN,1);
FreqsInfo = cell(NBlock,1);
for i = 1:NBlock
    seqTemp = zeros(recDurN,1);
    seqTempTDT = zeros(recDurN,1);
    for j = 1:ToneRep(i)
        N = times{i}.N(j);
        ramp = times{i}.ramp(j);
        Pulse_Dur = times{i}.dur(j);
        amp = Amp{i}(j);
        % Draw frequencies but don't use pure tones!
        f1 = rand() * (MaxFreq - MinFreq) + MinFreq;
        while 1
            f2  = rand() * (MaxFreq - MinFreq) + MinFreq;
            if f2 ~= f1
                break;
            end
        end
        if f2 > f1
            FreqsInfo{i}.fmin(j,1) = f1;
            FreqsInfo{i}.fmax(j,1) = f2;
            FreqsInfo{i}.rate(j,1) = log2(f2/f1)/Pulse_Dur;            
            FreqsInfo{i}.Sdir{j,1} = 'up';
        else
            FreqsInfo{i}.fmin(j,1) = f2;
            FreqsInfo{i}.fmax(j,1) = f1;
            FreqsInfo{i}.rate(j,1) = log2(f2/f1)/Pulse_Dur;
            FreqsInfo{i}.Sdir{j,1} = 'down';
        end
        t = linspace(0,Pulse_Dur,N)';
        signal = chirp(t,f1,Pulse_Dur,f2,SweepType);
        Gate = Cos2Gate(Fs,ramp,'s'); % onset and offset gate sequences
        signal = applyGate(signal,Gate); % apply gate sequence to onset and offset of the stimulus
        sigRMS = rms(signal);
        sigVolts = sigRMS*sqrt(2);
        signal = signal/sigVolts; % normalize to 1
        if exist('calib_V','var')
            a = 10^((amp-calib_level)/20)*calib_V;
            signalTDT = a*conv(signal ,FIRcoeff,'same');
        else
            signalTDT = signal;
        end
        
        seqTemp(times{i}.inds(j):times{i}.inds(j)+N-1,1) = ...
            seqTemp(times{i}.inds(j):times{i}.inds(j)+N-1,1) + signal; % N being the samples for each tone burst
        seqTempTDT(times{i}.inds(j):times{i}.inds(j)+N-1,1) = ...
            seqTempTDT(times{i}.inds(j):times{i}.inds(j)+N-1,1) + signalTDT; % N being the samples for each tone burst
    end
    Stim = Stim+seqTemp;
    StimTDT = StimTDT+seqTempTDT;
end

%% -----------------------------------------------
% Plotting
% -----------------------------------------------
% plot estimated stimulus spectrogram
cmap = load('selfmade3.txt'); cmap = cmap/255;
TimeCam   = (1:round(FsCam*recDur))/FsCam;    % camera time vector
M = max(Amps); m = min(Amps);
if m == M
    m = M - 30;
end
M = M+10;
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
imagesc(spTime, spFreq, spPdB,[m M])
axis xy; colormap(cmap); colorbar;
xlabel('Time (s)', 'FontSize', 16), xlim([0 recDur])
ylabel('Frequency (Hz)', 'FontSize', 16), ylim([MinFreq MaxFreq])
set(gca(), 'FontSize', 14)
if saveFlag
%     saveas(gcf(), fullfile(sDataDir, 'Spectrogram_Estimated'), 'pdf');
    figName = 'Spectrogram_Estimated';
    dpi = 600;
    pdfFilePath = sprintf('%s\\%s.pdf',sDataDir,figName);
%     save2pdf(pdfFilePath,gcf,dpi);
end
% Estimated Frequency-wise autocorrelation function
MaxLagSec = 0.5;
nMaxLag    = ceil(MaxLagSec * nFsSpec);
Corr      = zeros(size(spP,1), 2*nMaxLag+1);
tLags = (-nMaxLag:nMaxLag) / nFsSpec;
for iFreq=1:size(spP,1)
    Corr(iFreq,:) = xcorr(abs(spP(iFreq,:)), abs(spP(iFreq,:)), nMaxLag);
    Corr(iFreq,:) = Corr(iFreq,:) / max(Corr(iFreq,:));
end
figure('Name','Frequency-wise stimulus autocorrelation');
imagesc(tLags, spFreq, Corr,[0 1]);
axis xy; colormap(cmap); colorbar;
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
StimLen  = recDur-minOnset-minOffset;        % stimulus length is s
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
spec.FsCam        = FsCam;       % Camera sampling frequency in Hz
spec.SpeakerName  = SpeakerName; % speaker name
spec.Amp          = Amp;         % Mean sound pressure level
spec.times        = times;
spec.ToneRep      = ToneRep;
spec.seed         = seedOrig;    % seed value for reproducability
spec.NBlock       = NBlock;
spec.SweepType    = SweepType;
spec.FreqsInfo    = FreqsInfo;
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
    save(matFile,'Stim','StimTDT','spec','Speaker','Spectrogram','EstCorr','-v7.3');
%     audiowrite(sprintf('%s\\Stim.wav', sDataDir), StimTDT, round(Fs));
end



