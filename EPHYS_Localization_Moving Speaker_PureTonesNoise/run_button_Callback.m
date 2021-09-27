function run_button_Callback(varargin)
clc
RP = TDTRP('Ephus_Localization_MovingSpeaker_GUI.rcx', 'RX6');

min_Freq_h = findobj('tag','min_Freq');
max_Freq_h = findobj('tag','max_Freq');
Freq_per_Oct_h = findobj('tag','Freq_per_Oct');
Amp_Array_h = findobj('tag','Amp_Array');
ramp_h = findobj('tag','ramp');
Pulse_Dur_h = findobj('tag','Pulse_Dur');
ISI_h = findobj('tag','ISI');
jitter_percentage_h = findobj('tag','jitter_percentage');
Ntrials_h = findobj('tag','Ntrials');
NCam_h = findobj('tag','NCam');
Randomized_h = findobj('tag','Randomized');
Noise_h = findobj('tag','Noise');
seed_h = findobj('tag','seed');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');
FsCam_h = findobj('tag','FsCam');
nFrames_h = findobj('tag','nFrames');
nContStimB_h = findobj('tag','nContStimB');
nContStimE_h = findobj('tag','nContStimE');
% Get the specifications from GUI
min_Freq = get(min_Freq_h,'string');
min_Freq = str2double(min_Freq);
min_Freq(isnan(min_Freq)) = 4; % if the box is empty then this is default value.
max_Freq = get(max_Freq_h,'string');
max_Freq = str2double(max_Freq);
max_Freq(isnan(max_Freq)) = 32; % if the box is empty then this is default value.
Freq_per_Oct = get(Freq_per_Oct_h,'string'); %Freqency per Octave
Freq_per_Oct =  str2double(Freq_per_Oct);
Freq_per_Oct(isnan(Freq_per_Oct)) = 0.667; % if the box is empty then this is default value.
Fspec = [min_Freq, max_Freq, Freq_per_Oct];
Amp_Array = get(Amp_Array_h,'string');
Amp_Array = str2num(Amp_Array);
Amp_Array(isnan(Amp_Array)) = 75; % if the box is empty then this is default value.
ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 5; % if the box is empty then this is default value.
Pulse_Dur = get(Pulse_Dur_h,'string');
Pulse_Dur = str2double(Pulse_Dur);
Pulse_Dur(isnan(Pulse_Dur)) = 100; % if the box is empty then this is default value.
ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 4000; % if the box is empty then this is default value.
jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;
Ntrials = get(Ntrials_h,'string');
Ntrials = str2double(Ntrials);
Ntrials(isnan(Ntrials)) = 10; % if the box is empty then this is default value.
NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 7; % if the box is empty then this is default value.
Randomized = get(Randomized_h, 'value');
Noise = get(Noise_h, 'value');
seed = get(seed_h,'string');
seed = str2num(seed);
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);
nContStimB = get(nContStimB_h,'string');
nContStimB = str2double(nContStimB);
nContStimB(isnan(nContStimB)) = 2; % if the box is empty then this is default value.
nContStimE = get(nContStimE_h,'string');
nContStimE = str2double(nContStimE);
nContStimE(isnan(nContStimE)) = 0; % if the box is empty then this is default value.

% load Tone-Azimuth-Time orders
FolderPath = [pwd,'\\Tone-Azimuth-Time Orders'];
fileNameMat = sprintf('Sound Localization_Moving Speaker_Seed_%d.mat',seed);
filePathMat = sprintf('%s\\%s',FolderPath,fileNameMat);
load(filePathMat,'NTones','NperTone','NAz','ToneOrder','AzOrder','TimeOrder')

% nContStimB = 0
NncCam = NCam-nContStimB-nContStimE; % number of camera recordings for non-cont stim
Ntotal = NTones*NperTone*NAz; % number of all epochs of stims for non-cont recordins
NperCam = Ntotal/NncCam; % number of stims per camera recordings

% get device sampling frequency
Fs = RP.FS;
% Fs = 24414.0625*8;
% Fs = round(Fs);
[calib_V, calib_level, ~, ~, FIRcoeff] = getCalib(SpeakerName,Fs);
F = octavefreq(Fspec);
F = 0.1*round(10*F);
Amp_TDT = getTDTAmp(F, F, repmat(Amp_Array,1,length(F)), SpeakerName);
F = F*1000;
recDur = nFrames/FsCam;
% Make sounds for continuous stims
rng(seed);
delayBetweenConStim = 2;
Dur = round((recDur-2*delayBetweenConStim)/2); % in sec
% Dur = 115.4;
Dur = Dur*1000;  % in msec
for k = 1:length(F)
    f = F(k);
    tone = gTone(Fs,f,Dur,ramp);
    cSigs(k).Signal = Amp_TDT(k)*tone;
    cSigs(k).N = length(tone);
    cSigs(k).dur = length(tone)/Fs;
    cSigs(k).name = sprintf('%dkHz',round(f/1000));
end
if Noise
    k = length(F)+1;
    [noiseG, N] = WideBandNoise(Fs,Dur,ramp);
    % apply calibration filter and adjust amplitudes for required dB SPL
    Amp = Amp_Array;
    sigRMS = rms(noiseG);
    sigVolts = sigRMS*sqrt(2);
    a = 10^((Amp-calib_level)/20)*calib_V;
    noiseG = noiseG*a/sigVolts;
    % apply calibration filter and adjust amplitudes for required dB SPL
    noiseGC = conv(noiseG,FIRcoeff,'same');
    cSigs(k).Signal = noiseGC;
    cSigs(k).N = N;
    cSigs(k).dur = N/Fs;
    cSigs(k).name = 'WBNoise';
end
% Make sounds for non-continuous stims
for k = 1:length(F)
    f = F(k);
    tone = gTone(Fs,f,Pulse_Dur,ramp);
    ncSigs(k).Signal = Amp_TDT(k)*tone;
    ncSigs(k).N = length(tone);
    ncSigs(k).dur = length(tone)/Fs;
    ncSigs(k).name = sprintf('%dkHz',round(f/1000));
    RecNames{k} = ncSigs(k).name;
end
if Noise
    k = length(F)+1;
    for n = 1:NperTone*NAz
        [noiseG, N] = WideBandNoise(Fs,Pulse_Dur,ramp);
        % apply calibration filter and adjust amplitudes for required dB SPL
        Amp = Amp_Array;
        sigRMS = rms(noiseG);
        sigVolts = sigRMS*sqrt(2);
        a = 10^((Amp-calib_level)/20)*calib_V;
        noiseG = noiseG*a/sigVolts;
        % apply calibration filter and adjust amplitudes for required dB SPL
        noiseGC = conv(noiseG,FIRcoeff,'same');
        ncSigs(k).Signal{n} = noiseGC;
        ncSigs(k).N = N;
        ncSigs(k).dur = N/Fs;
        ncSigs(k).name = 'WBNoise';
        RecNames{k} = ncSigs(k).name;
    end
end
%% First continuous recording
npts = 200000;
RP.write('SourceSize', npts);
bufpts = npts/2;
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
for csB = 1:nContStimB
    fprintf('\nContinuous Recording #%d\n',csB);
    % load stim sequence
    k = 2*csB-1;
    signal = 0.5*cSigs(k).Signal;
    N = cSigs(k).N;
    RP.write('N', N);
    % load up entire buffer with segments A and B
    index = 1;
    top = min(index+npts, N);
    if top == N
        top = top+1;
    end
    RP.write('datain', signal(index:top-1));
    index = top;
    % waiting for trigger from other hardware
    fprintf('TDT is waiting for the trigger from RPi\n')
    while ~RP.read('PiTrigger')
    end
    RP.trg(1);
    fprintf('Trigger received. Recording in progress...\n\n');
    try
        stimName = ncSigs(k).name;
    catch
        stimName = ncSigs(k).name{1};
    end
    fprintf('%s\n',stimName)
    curindex = RP.read('index');
    while index < N
        % wait until done playing A
        while(curindex < bufpts)
            curindex = RP.read('index');
            pause(.05);
        end
        % load the next signal segment
        top = min(index+bufpts, N);
        if top == N
            top = top+1;
        end
        RP.write('datain', signal(index:top-1));
        index = top;
        % checks to see if the data transfer rate is fast enough
        curindex = RP.read('index');
        if(curindex < bufpts)
            warning('Transfer rate is too slow');
        end
        % wait until start playing A
        while(curindex > bufpts)
            curindex = RP.read('index');
            pause(.05);
        end
        % load segment B
        if index > N
            break;
        end
        top = min(index+bufpts, N);
        if top == N
            top = top+1;
        end
        RP.write('datain', signal(index:top-1), 'OFFSET', bufpts);
        index = top;
        % make sure we're still playing A
        curindex = RP.read('index');
        if(curindex > bufpts)
            warning('Transfer rate too slow');
        end
    end
%     pause(npts/Fs); % make sure that buffer is palyed out

    % second part 
    % waiting for trigger from other hardware
    fprintf('TDT is waiting for the trigger from RPi\n')
    while ~RP.read('PiTrigger')
    end
    fprintf('Trigger received from Pi. Carriage is stopped...\n\n');
    pause(0.2)
    % load stim sequence
    k = 2*csB;
    signal = cSigs(k).Signal;
    N = cSigs(k).N;
    RP.write('N', N);
    % load up entire buffer with segments A and B
    index = 1;
    top = min(index+npts, N);
    if top == N
        top = top+1;
    end
    RP.write('datain', signal(index:top-1));
    index = top;
    fprintf('TDT is waiting for the trigger from RPi\n')
    while ~RP.read('PiTrigger')
    end
    RP.trg(1);
    fprintf('Trigger received from Pi. Carriage is moving. Second Part of Stim in progress...\n\n');
    try
        stimName = ncSigs(k).name;
    catch
        stimName = ncSigs(k).name{1};
    end
    fprintf('%s\n',stimName)
    curindex = RP.read('index');
    while index < N
        % wait until done playing A
        while(curindex < bufpts)
            curindex = RP.read('index');
            pause(.05);
        end
        % load the next signal segment
        top = min(index+bufpts, N);
        if top == N
            top = top+1;
        end
        RP.write('datain', signal(index:top-1));
        index = top;
        % checks to see if the data transfer rate is fast enough
        curindex = RP.read('index');
        if(curindex < bufpts)
            warning('Transfer rate is too slow');
        end
        % wait until start playing A
        while(curindex > bufpts)
            curindex = RP.read('index');
            pause(.05);
        end
        % load segment B
        if index > N
            break;
        end
        top = min(index+bufpts, N);
        if top == N
            top = top+1;
        end
        RP.write('datain', signal(index:top-1), 'OFFSET', bufpts);
        index = top;
        % make sure we're still playing A
        curindex = RP.read('index');
        if(curindex > bufpts)
            warning('Transfer rate too slow');
        end
    end
%     pause(npts/Fs); % make sure that buffer is palyed out
    % waiting for trigger from other hardware
    fprintf('TDT is waiting for the trigger from RPi\n')
    while ~RP.read('PiTrigger')
    end
    fprintf('Trigger received from Pi. Carriage is stopped...\n\n');
    pause(0.2)
end
%%
% play out the sounds
Nnoise = 0; % a counter for noise realizations
for nc = 1:NncCam
    idx = 1+(nc-1)*NperCam : nc*NperCam;
    ToneOrdertemp = ToneOrder(idx);
    delays = TimeOrder(idx);
    tplay = 0;
    % waiting for trigger from other hardware
%     disp('TDT is waiting for the trigger from A-M system')
%     while ~RP.read('Trigger')
%     end
    tot = tic;
    for n = 1:length(idx)
        % load up entire buffer
        tic
        k = ToneOrdertemp(n);
        if k == length(F)+1
            Nnoise = Nnoise+1;
            signal = ncSigs(k).Signal{Nnoise};
        else
            signal = ncSigs(k).Signal;
            
        end
        N = ncSigs(k).N;
        dur = ncSigs(k).dur;
        stimName = ncSigs(k).name;
        RP.write('datain', signal);
        RP.write('N', N);
        RP.write('SourceSize', N);
        fprintf('%8s; %s/%s\n',stimName,num2str((nc-1)*NperCam+n),num2str(Ntotal))
        % wait for trigger from raspberry pi. This means that speaker is
        % located at the next azimuth
        while ~RP.read('PiTrigger')
        end
%         tPi = toc;
        % pause(delays(k)-tplay-tPi-dur)
        
        % start playing
        RP.trg(1);
        tic
        pause(0.5*dur)
        curindex = RP.read('index');
        while(curindex > 0) && (curindex < N)
            curindex = RP.read('index');
        end
        tplay = toc;
    end
    toc(tot)
    fprintf('\n');
end

pause(ISI/1000);
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
RP.halt;
disp('Done Dude!')

%% Saving Stimulus Set Information
StimSet.seed = seed;
StimSet.RecNames = RecNames;

StimSet.Amp = Amp_Array;
StimSet.NCam = NCam;
StimSet.nContStimB = nContStimB;
StimSet.nContStimE = nContStimE;

StimSet.NTones = NTones;
StimSet.NperTone = NperTone;
StimSet.NAz = NAz;
StimSet.ToneOrder = ToneOrder;
StimSet.AzOrder = AzOrder;
StimSet.TimeOrder = TimeOrder;
StimSet.NperCam = NperCam;


StimSet.ISI = ISI;
StimSet.jitter_percentage = jitter_percentage;
StimSet.Randomized = Randomized;
StimSet.Noise = Noise;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;
StimSet.FsTDT = Fs;
StimSet.nFrames = nFrames;

StimSet.distinctTones = NTones;

tot_trials = length(ToneOrder);
for n = 1:tot_trials
    k = ToneOrder(n);
	StimSet.trials(n).StimName = ncSigs(k).name;
	StimSet.trials(n).StimNum = ToneOrder(n);
	StimSet.trials(n).Azimuth = AzOrder(n);
	StimSet.trials(n).Time = TimeOrder(n);
end

fileno = str2double(fileno); fileno = round(fileno);

Amp_str = [];
for k = 1:length(Amp_Array)
	Amp_str = [Amp_str, num2str(Amp_Array(k)),'_'];
end

filename = ['0',num2str(fileno),'_MovingSpeakerLocalization',...
    '_NTones_',num2str(NTones),...
	'_Amp_',Amp_str,...
	'ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
	'_Ntrials_',num2str(Ntrials)];

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



