function [SweepTimeOrder, SweepDirOrder, rampTimeOrder, AmpOrder, TimeOrder] = run_button_Callback(varargin)
clc
RP = TDTRP('Ephus_FMSweeps_GUI.rcx', 'RX6');

min_Freq_h = findobj('tag','min_Freq');
max_Freq_h = findobj('tag','max_Freq');
Amp_Array_h = findobj('tag','Amp_Array');
min_Pulse_Dur_h = findobj('tag','min_Pulse_Dur');
max_Pulse_Dur_h = findobj('tag','max_Pulse_Dur');
Pulse_Dur_Levels_h = findobj('tag','Pulse_Dur_Levels');
min_ramp_h = findobj('tag','min_ramp');
max_ramp_h = findobj('tag','max_ramp');
ISI_h = findobj('tag','ISI');
jitter_percentage_h = findobj('tag','jitter_percentage');
Ntrials_h = findobj('tag','Ntrials');
Randomized_h = findobj('tag','Randomized');
UpDown_h = findobj('tag','UpDown');
Log_Levels_h = findobj('tag','Log_Levels');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');
FsCam_h = findobj('tag','FsCam');
nFrames_h = findobj('tag','nFrames');
NCam_h = findobj('tag','NCam');

% Get the specifications from GUI
min_Freq = get(min_Freq_h,'string');
min_Freq = str2double(min_Freq);
min_Freq(isnan(min_Freq)) = 4; % if the box is empty then this is default value.

max_Freq = get(max_Freq_h,'string');
max_Freq = str2double(max_Freq);
max_Freq(isnan(max_Freq)) = 32; % if the box is empty then this is default value.

Fspec = [min_Freq, max_Freq];

Amp_Array = get(Amp_Array_h,'string');
Amp_Array = str2num(Amp_Array);
if isempty(Amp_Array)
	Amp_Array = 60; % if the box is empty then this is default value.
end

min_Pulse_Dur = get(min_Pulse_Dur_h,'string'); %Freqency per Octave
min_Pulse_Dur =  str2double(min_Pulse_Dur);
min_Pulse_Dur(isnan(min_Pulse_Dur)) = 7.5; % if the box is empty then this is default value.

max_Pulse_Dur = get(max_Pulse_Dur_h,'string'); %Freqency per Octave
max_Pulse_Dur =  str2double(max_Pulse_Dur);
max_Pulse_Dur(isnan(max_Pulse_Dur)) = 600; % if the box is empty then this is default value.

Pulse_Dur_Levels = get(Pulse_Dur_Levels_h,'string'); %Freqency per Octave
Pulse_Dur_Levels =  str2double(Pulse_Dur_Levels);
Pulse_Dur_Levels(isnan(Pulse_Dur_Levels)) = 5; % if the box is empty then this is default value.

Pulse_Dur_spec = [min_Pulse_Dur, max_Pulse_Dur, Pulse_Dur_Levels];

min_ramp = get(min_ramp_h,'string'); % min ramp in msec
min_ramp =  str2double(min_ramp);
min_ramp(isnan(min_ramp)) = 1; % if the box is empty then this is default value.

max_ramp = get(max_ramp_h,'string'); % max ramp in msec
max_ramp =  str2double(max_ramp);
max_ramp(isnan(max_ramp)) = 5; % if the box is empty then this is default value.

ramp_spec = [min_ramp, max_ramp, Pulse_Dur_Levels];

ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 5000; % if the box is empty then this is default value.

jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;

Ntrials = get(Ntrials_h,'string');
Ntrials = str2double(Ntrials);
Ntrials(isnan(Ntrials)) = 25; % if the box is empty then this is default value.

NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 5; % if the box is empty then this is default value.

Randomized = get(Randomized_h, 'value');
UpDown  = get(UpDown_h, 'value');
Log_Levels  = get(Log_Levels_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);

% get device sampling frequency
fs = RP.GetSFreq();
% fs = 24414.0625*8;
[calib_V, calib_level, ~, ~, FIRcoeff] = getCalib(SpeakerName,fs);

[SweepTimes, SweepTimeOrder, SweepDirOrder, rampTimeOrder, AmpOrder, TimeOrder] = ...
    SATorder(Pulse_Dur_spec, ramp_spec, Amp_Array, ISI, jitter, Ntrials, ...
             Randomized, UpDown, Log_Levels, FsCam, nFrames, NCam);

SweepOct = log2(max_Freq/min_Freq);
SweepRateOrder = 1000*SweepOct./SweepTimeOrder;
% play out the sounds
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
nS = length(SweepTimes);
nSD = length(unique(SweepDirOrder));
nA = length(Amp_Array);
NTones = nS*nSD*nA;
totTrials = NTones*Ntrials;
NperCam = totTrials/NCam;
for nc = 1:NCam
    idx = 1+(nc-1)*NperCam : nc*NperCam;
    STOrdertemp = SweepTimeOrder(idx);
    SROrdertemp = SweepRateOrder(idx);
    SDOrdertemp = SweepDirOrder(idx);
    rOrdertemp = rampTimeOrder(idx);
    AOrdertemp = AmpOrder(idx);
    delays = TimeOrder(idx);
    tplay = 0;
    
    % waiting for trigger from other hardware
    disp('TDT is waiting for the trigger from A-M system')
    while ~RP.GetTagVal('Trigger')
    end
    totCam = tic;
    for n = 1:length(idx)
        Ttrial = tic;
        tic; % load timer
        % Generate FM sweep        
        Sdir = SDOrdertemp{n};
        Pulse_Dur = STOrdertemp(n);
        ramp = rOrdertemp(n);
        [FMSweepG, N] = FMchirp(fs,min_Freq,max_Freq,Sdir,Pulse_Dur,ramp);

        % apply calibration filter and adjust amplitudes for required dB SPL
        Amp = AOrdertemp(n);
        sigRMS = rms(FMSweepG);
        sigVolts = sigRMS*sqrt(2);
        a = 10^((Amp-calib_level)/20)*calib_V;
        FMSweepG = FMSweepG*a/sigVolts;
        FMSweepGC = conv(FMSweepG,FIRcoeff,'same');
        
        % load up entire buffer
        RP.WriteTagVEX('datain', 0, 'F32', FMSweepGC);
        RP.SetTagVal('PulseSamp', N);
        RP.SetTagVal('SourceSize', N);
        Rate = SROrdertemp(n);
        disp(sprintf('%s/%s. FM Sweep: %s oct/s; %sward; %skHz - %skHz; %s dB SPL',...
            num2str((nc-1)*NperCam+n), num2str(totTrials), num2str(round(Rate)), Sdir,...
            num2str(0.1*round(10*min_Freq)), num2str(0.1*round(10*max_Freq)), num2str(Amp)));
        tload = toc;
        pause(delays(n)-tload-tplay)
        
        % start playing
        RP.SoftTrg(1);
        tic
        pause(0.5*Pulse_Dur/1000)
        curindex = RP.GetTagVal('index');
        while(curindex > 0) && (curindex < N)
            curindex = RP.GetTagVal('index');
        end
        tplay = toc;
        TimeOrder((nc-1)*NperCam+n) = toc(Ttrial);
    end
    TCam = toc(totCam);
    disp(['Total time: ',num2str(TCam), 's']);
    fprintf('\n');
end

pause(ISI/1000);
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
RP.Halt;
disp('Done Dude!')

%% Saving Stimulus Set Information 
if UpDown
    SweepDirs = {'up', 'down'};
else
    SweepDirs = {'up'};
end
rampTimes = unique(rampTimeOrder');

StimSet.StimType = 'FM Sweep';
StimSet.Fspec = Fspec;
StimSet.SweepOct = SweepOct;
StimSet.SweepTimes = SweepTimes;
StimSet.rampTimes = rampTimes;
StimSet.SweepDirs = SweepDirs;
StimSet.Amp_Array = Amp_Array;
StimSet.ISI = ISI;
StimSet.jitter_percentage = jitter_percentage;
StimSet.Randomized = Randomized;
StimSet.UpDown = UpDown;
StimSet.LogarithmicLevels = Log_Levels;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;
StimSet.FsTDT = fs;
StimSet.nFrames = nFrames;
StimSet.distinctTones = NTones;
StimSet.NperCam = NperCam;
StimSet.NCam = NCam;

for n = 1:totTrials
	StimSet.trials(n).SweepTime = SweepTimeOrder(n);
	StimSet.trials(n).SweepDir = SweepDirOrder(n);
	StimSet.trials(n).SweepRate = SweepRateOrder(n);
	StimSet.trials(n).rampTime = rampTimeOrder(n);
	StimSet.trials(n).Amp = AmpOrder(n);
	StimSet.trials(n).Time = TimeOrder(n);
end

fileno = str2double(fileno); fileno = round(fileno);

Amp_str = [];
for k = 1:length(Amp_Array)
	Amp_str = [Amp_str, num2str(Amp_Array(k)),'_'];
end

filename = ['0',num2str(fileno),'_FMSweep', '_fmin_',num2str(min_Freq),...
	'_fmax_',num2str(max_Freq), '_Amp_',Amp_str,...
	'mindur_',num2str(min_Pulse_Dur), '_maxdur_',num2str(max_Pulse_Dur),...
    '_NLevels_',num2str(Pulse_Dur_Levels),...
    '_minramp_',num2str(min_ramp), '_maxramp_',num2str(max_ramp),...
	'_ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
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


