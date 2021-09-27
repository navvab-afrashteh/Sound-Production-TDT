function [VocalOrder, TimeOrder] = run_button_Callback_Vocal_Tinnitus(varargin)
clc
RP = TDTRP('Ephus_MiceVocalization_GUI.rcx', 'RX6');

Nvocal_h = findobj('tag','Nvocal');
Amp_h = findobj('tag','Amp');
ISI_h = findobj('tag','ISI');
jitter_percentage_h = findobj('tag','jitter_percentage');
Ntrials_h = findobj('tag','Ntrials');
NCam_h = findobj('tag','NCam');
Randomized_h = findobj('tag','Randomized');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');
FsCam_h = findobj('tag','FsCam');
nFrames_h = findobj('tag','nFrames');
nSponB_h = findobj('tag','nSponB');
nSponE_h = findobj('tag','nSponE');
% Get the specifications from GUI
Nvocal = get(Nvocal_h,'string');
Nvocal = str2double(Nvocal);
Nvocal(isnan(Nvocal)) = 10; % if the box is empty then this is default value.
Amp = get(Amp_h,'string');
Amp = str2double(Amp);
Amp(isnan(Amp)) = 85; % if the box is empty then this is default value.
ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 5000; % if the box is empty then this is default value.
jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;
nSponB = get(nSponB_h,'string');
nSponB = str2double(nSponB);
nSponB(isnan(nSponB)) = 2; % if the box is empty then this is default value.
nSponE = get(nSponE_h,'string');
nSponE = str2double(nSponE);
nSponE(isnan(nSponE)) = 0; % if the box is empty then this is default value.
Ntrials = get(Ntrials_h,'string');
Ntrials = str2double(Ntrials);
Ntrials(isnan(Ntrials)) = 25; % if the box is empty then this is default value.
NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 7; % if the box is empty then this is default value.
Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);

% get device sampling frequency
fs = RP.FS;
% fs = 24414.0625*8;
[calib_V, calib_level, ~, ~, FIRcoeff] = getCalib(SpeakerName,fs);

% Read vocalization files and resample
k = 0;
for n = 1:Nvocal
    k = k+1;
    % PathName = sprintf('%s\\Sine Waves',pwd);
    % filePath = sprintf('%s\\sineWave_%dkHz.wav',PathName,k+1);
    PathName = sprintf('%s\\Mice Vocalization',pwd);
    filePath = sprintf('%s\\P100_%d.wav',PathName,n);

    [d, fsold] = audioread(filePath);
    if abs(fsold-fs)<1; fsold = fs; end
    [p,q] = rat(fs/fsold);
    signals(k).Signal = resample(d,p,q);
    signals(k).fs = fs;
    signals(k).N = length(signals(k).Signal);
    signals(k).T = signals(k).N/signals(k).fs;
    signals(k).t = linspace(0,signals(k).T,signals(k).N);
    
    % apply calibration filter and adjust amplitudes for required dB SPL
    sigRMS = rms(signals(k).Signal);
    sigVolts = sigRMS*sqrt(2);
    a = 10^((Amp-calib_level)/20)*calib_V;
    signals(k).Signal = signals(k).Signal*a/sigVolts;
    signals(k).Signal = conv(signals(k).Signal,FIRcoeff,'same');
    
    sigRMS = rms(signals(k).Signal);
    sigVolts = sigRMS*sqrt(2);
    
    M = max(abs(signals(k).Signal));
    if M>10
%         signals(k).Signal = 10*signals(k).Signal/M;
    end
end

nEvoked = NCam-nSponB-nSponE;
[VocalOrder, TimeOrder] = VocalTorder(Nvocal, nEvoked, ISI, jitter, Ntrials, Randomized,FsCam, nFrames);

[RecNames{1:nEvoked}] = deal('Vocal');
[SponBNames{1:nSponB}] = deal('Spon');
RecNames = [SponBNames, RecNames];
[SponENames{1:nSponE}] = deal('Spon');
RecNames = [RecNames, SponENames];

%% First spontaneous recording
for spB = 1:nSponB
    recDur = nFrames/FsCam;
    % waiting for trigger from other hardware
    fprintf('\nSpontaneous Recording:\n');
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~RP.read('Trigger')
    end
    start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fprintf('Trigger received. Recording in progress...\n\n');
    pause(recDur);
end

%% play out the sounds
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
for nc = 1:nEvoked
    idx = 1+(nc-1)*Nvocal*Ntrials/nEvoked : nc*Nvocal*Ntrials/nEvoked;
    VocalOrdertemp = VocalOrder(idx);
    delays = TimeOrder(idx);
    tplay = 0;
    
    % waiting for trigger from other hardware
    disp('TDT is waiting for the trigger from A-M system')
    while ~RP.read('Trigger')
    end
    tot = tic;
    for n = 1:length(idx)
        % load up entire buffer
        tic
        k = VocalOrdertemp(n);
        RP.write('datain', signals(k).Signal);
        RP.write('PulseSamp', signals(k).N);
        RP.write('SourceSize', signals(k).N);
        disp(['Vocalization #', num2str(k),'; ',num2str((nc-1)*Nvocal*Ntrials/nEvoked+n),'/',num2str(Nvocal*Ntrials)]);
        tload = toc;
        pause(delays(n)-tplay-tload)
        
        % start playing
        RP.trg(1);
        tic
        pause(0.5*signals(k).T)
        curindex = RP.read('index');
        endInd = signals(k).N;
        while(curindex > 0) && (curindex < endInd)
            curindex = RP.read('index');
        end
        tplay = toc;
    end
    toc(tot)
    fprintf('\n');
end

pause(ISI/1000);

%% End spontanious recordings
for spE = 1:nSponE
    recDur = nFrames/FsCam;
    % waiting for trigger from other hardware
    fprintf('\n Spontaneous Recording:\n');
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~RP.read('Trigger')
    end
    start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fprintf('Trigger received. Recording in progress...\n\n');
    pause(recDur);
end
%% stop playing
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
RP.halt;
disp('Done Dude!')

%% Saving Stimulus Set Information
StimSet.StimType    = 'Vocalization';
StimSet.NCam = NCam;
StimSet.RecNames = RecNames;
StimSet.nEvoked = nEvoked;
StimSet.nSponB = nSponB;
StimSet.nSponE = nSponE;
StimSet.signals = signals;
StimSet.Nvocal = Nvocal;
StimSet.Amp = Amp;
StimSet.ISI = ISI;
StimSet.jitter_percentage = jitter_percentage;
StimSet.Randomized = Randomized;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;
StimSet.nFrames = nFrames;
StimSet.distinctTones = Nvocal;

tot_trials = length(VocalOrder);
for n = 1:tot_trials
	StimSet.trials(n).VocalSignal = VocalOrder(n);
	StimSet.trials(n).Time = TimeOrder(n);
end
fileno = str2double(fileno); fileno = round(fileno);
Amp_str = [];
for k = 1:length(Amp)
	Amp_str = [Amp_str, num2str(Amp(k)),'_'];
end
if nSponB>0 && nSponE>0
    filename = ['0',num2str(fileno),'_Vocalization_',num2str(nSponB),'Spon_',...
        num2str(nEvoked),'Evoked_',num2str(nSponE),'Spon',...
        '_Nvocal_',num2str(Nvocal),'_Amp_',Amp_str,...
        'ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
        '_Ntrials_',num2str(Ntrials)];
elseif nSponB>0 && nSponE==0
    filename = ['0',num2str(fileno),'_Vocalization_',num2str(nSponB),'Spon_',...
        num2str(nEvoked),'Evoked',...
        '_Nvocal_',num2str(Nvocal),'_Amp_',Amp_str,...
        'ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
        '_Ntrials_',num2str(Ntrials)];
elseif nSponB==0 && nSponE>0
    filename = ['0',num2str(fileno),'_Vocalization_',...
        num2str(nEvoked),'Evoked_',num2str(nSponE),'Spon',...
        '_Nvocal_',num2str(Nvocal),'_Amp_',Amp_str,...
        'ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
        '_Ntrials_',num2str(Ntrials)];
end

if fileno <= 9.5
	eval(['StimSet_0',num2str(fileno), '=', 'StimSet']);
	save([filename,'.mat'], ['StimSet_0',num2str(fileno)]);
	
	if fileno >= 8.5
		set(fileno_h,'string',num2str(fileno+1));
	else
		set(fileno_h,'string',['0',num2str(fileno+1)]);
	end
else
	eval(['StimSet_',num2str(fileno), '=', 'StimSet']);
	save([filename(2:end),'.mat'],['StimSet_',num2str(fileno)]);
	set(fileno_h,'string',num2str(fileno+1));
end



