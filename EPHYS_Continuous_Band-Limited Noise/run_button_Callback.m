function [FreqOrder, AmpOrder, TimeOrder] = run_button_Callback(varargin)
clc
RP = TDTRP('Ephus_BandLimitedNoise_GUI.rcx', 'RX6');

min_Freq_h = findobj('tag','min_Freq');
max_Freq_h = findobj('tag','max_Freq');
Freq_per_Oct_h = findobj('tag','Freq_per_Oct');
BW_Oct_h = findobj('tag','BW_Oct');
Mod_Freq_h = findobj('tag','Mod_Freq');
Amp_Array_h = findobj('tag','Amp_Array');
ramp_h = findobj('tag','ramp');
Pulse_Dur_h = findobj('tag','Pulse_Dur');
ISI_h = findobj('tag','ISI');
jitter_percentage_h = findobj('tag','jitter_percentage');
Ntrials_h = findobj('tag','Ntrials');
Randomized_h = findobj('tag','Randomized');
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

Freq_per_Oct = get(Freq_per_Oct_h,'string'); %Freqency per Octave
Freq_per_Oct =  str2double(Freq_per_Oct);
Freq_per_Oct(isnan(Freq_per_Oct)) = 1; % if the box is empty then this is default value.

Fspec = [min_Freq, max_Freq, Freq_per_Oct];

BW_Oct = get(BW_Oct_h,'string'); % Bandwidth in Octave
BW_Oct =  str2double(BW_Oct);
BW_Oct(isnan(BW_Oct)) = 1; % if the box is empty then this is default value.

Mod_Freq = get(Mod_Freq_h,'string'); % Modulation frequency for amplitude modulation (in Hz)
Mod_Freq =  str2double(Mod_Freq);
Mod_Freq(isnan(Mod_Freq)) = 10; % if the box is empty then this is default value.

Amp_Array = get(Amp_Array_h,'string');
Amp_Array = str2num(Amp_Array);
if isempty(Amp_Array)
	Amp_Array = 40:20:80; % if the box is empty then this is default value.
end
ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 2; % if the box is empty then this is default value.

Pulse_Dur = get(Pulse_Dur_h,'string');
Pulse_Dur = str2double(Pulse_Dur);
Pulse_Dur(isnan(Pulse_Dur)) = 300; % if the box is empty then this is default value.

ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 5000; % if the box is empty then this is default value.

jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;

Ntrials = get(Ntrials_h,'string');
Ntrials = str2double(Ntrials);
Ntrials(isnan(Ntrials)) = 20; % if the box is empty then this is default value.

NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 5; % if the box is empty then this is default value.

Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);

% get device sampling frequency
fs = RP.GetSFreq();
% fs = 24414.0625*8;
[calib_V, calib_level, ~, ~, FIRcoeff] = getCalib(SpeakerName,fs);

[F, FreqOrder, AmpOrder, TimeOrder] = FATorder(Fspec, Amp_Array,...
	ISI, jitter, Ntrials, Randomized, FsCam, nFrames, NCam); 

% play out the sounds
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
nF = length(F);
nA = length(Amp_Array);
NTones = nF*nA;
totTrials = NTones*Ntrials;
NperCam = totTrials/NCam;
for nc = 1:NCam
    idx = 1+(nc-1)*NperCam : nc*NperCam;
    FOrdertemp = FreqOrder(idx);
    AOrdertemp = AmpOrder(idx);
    delays = TimeOrder(idx);
    tplay = 0;
    
    % waiting for trigger from other hardware
    disp('TDT is waiting for the trigger from A-M system')
%     while ~RP.GetTagVal('Trigger')
%     end
    totCam = tic;
    for n = 1:length(idx)
        Ttrial = tic;
        tic; % load timer
        % Generate band-limited noise
        f = FOrdertemp(n)*1000; % convert to Hz
        [BLnoiseGAM, N, f1, f2] = BandLimNoise(fs,f,BW_Oct,Pulse_Dur,ramp,Mod_Freq);
        % apply calibration filter and adjust amplitudes for required dB SPL
        Amp = AOrdertemp(n);
        sigRMS = rms(BLnoiseGAM);
        sigVolts = sigRMS*sqrt(2);
        a = 10^((Amp-calib_level)/20)*calib_V;
        BLnoiseGAM = BLnoiseGAM*a/sigVolts;
        BLnoiseGAMC = conv(BLnoiseGAM,FIRcoeff,'same');
        
        % load up entire buffer
        RP.WriteTagVEX('datain', 0, 'F32', BLnoiseGAMC);
        RP.SetTagVal('PulseSamp', N);
        RP.SetTagVal('SourceSize', N);
        disp(sprintf('%s/%s. Band-Limited Noise: %skHz - %skHz; %s dB SPL', num2str((nc-1)*NperCam+n), num2str(totTrials),num2str(0.1*round(10*f1/1000)), num2str(0.1*round(10*f2/1000)), num2str(Amp)))
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
StimSet.Fspec = Fspec;
StimSet.BW_Oct = BW_Oct;
StimSet.Mod_Freq = Mod_Freq;
StimSet.Amp_Array = Amp_Array;
StimSet.ramp = ramp;
StimSet.Pulse_Dur = Pulse_Dur;
StimSet.NCam = NCam;
StimSet.ISI = ISI;
StimSet.jitter_percentage = jitter_percentage;
StimSet.Randomized = Randomized;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;
StimSet.nFrames = nFrames;
StimSet.distinctTones = NTones;

for n = 1:totTrials
	StimSet.trials(n).Freq = FreqOrder(n);
	StimSet.trials(n).Amp = AmpOrder(n);
	StimSet.trials(n).Time = TimeOrder(n);
end

fileno = str2double(fileno); fileno = round(fileno);

Amp_str = [];
for k = 1:length(Amp_Array)
	Amp_str = [Amp_str, num2str(Amp_Array(k)),'_'];
end

filename = ['0',num2str(fileno),'_BandLimNoise', '_fmin_',num2str(min_Freq),...
	'_fmax_',num2str(max_Freq),'_FreqPerOct_',num2str(Freq_per_Oct),'_BWoct_',num2str(BW_Oct),...
	'_Amp_',Amp_str,'dur_',num2str(Pulse_Dur),'_ramp_',num2str(ramp),'_AMfreq_',num2str(Mod_Freq),...
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


