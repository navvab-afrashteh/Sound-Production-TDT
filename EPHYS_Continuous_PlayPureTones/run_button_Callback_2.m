function [FreqOrder, AmpOrder, TimeOrder] = run_button_Callback_2(varargin)
clc
RP = TDTRP('Ephus_PlayPureTones_GUI_2.rcx', 'RX6');
% RP.Halt;

min_Freq_h = findobj('tag','min_Freq');
max_Freq_h = findobj('tag','max_Freq');
Freq_per_Oct_h = findobj('tag','Freq_per_Oct');
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

% Get the specifications from GUI
min_Freq = get(min_Freq_h,'string');
min_Freq = str2double(min_Freq);
min_Freq(isnan(min_Freq)) = 2; % if the box is empty then this is default value.

max_Freq = get(max_Freq_h,'string');
max_Freq = str2double(max_Freq);
max_Freq(isnan(max_Freq)) = 64; % if the box is empty then this is default value.

Freq_per_Oct = get(Freq_per_Oct_h,'string'); %Freqency per Octave
Freq_per_Oct =  str2double(Freq_per_Oct);
Freq_per_Oct(isnan(Freq_per_Oct)) = 2; % if the box is empty then this is default value.

Fspec = [min_Freq, max_Freq, Freq_per_Oct];

Amp_Array = get(Amp_Array_h,'string');
Amp_Array = str2num(Amp_Array);
if isempty(Amp_Array)
	Amp_Array = 10:10:80; % if the box is empty then this is default value.
end

ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 5; % if the box is empty then this is default value.

Pulse_Dur = get(Pulse_Dur_h,'string');
Pulse_Dur = str2double(Pulse_Dur);
Pulse_Dur(isnan(Pulse_Dur)) = 50; % if the box is empty then this is default value.

ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 2500; % if the box is empty then this is default value.

jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;

Ntrials = get(Ntrials_h,'string');
Ntrials = str2double(Ntrials);
Ntrials(isnan(Ntrials)) = 5; % if the box is empty then this is default value.

Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);

[F, FreqOrder, AmpOrder, TimeOrder] = FATorder(Fspec, Amp_Array, ramp,...
	Pulse_Dur, ISI, jitter,...
	Ntrials, Randomized, FsCam);

Amp_TDT = getTDTAmp(F, FreqOrder, AmpOrder, SpeakerName);
FreqOrder = FreqOrder*1000;

tot_trials = size(FreqOrder,1);
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
NTones = tot_trials/Ntrials;

FsTDT = RP.GetSFreq();
RP.SetTagVal('Fs',FsTDT);
RP.SetTagVal('NTones',NTones);
RP.SetTagVal('PulseDurTg',Pulse_Dur-ramp);

% TimeOrder = (1:tot_trials)*1;
for trial = 1:Ntrials
    idx = 1+(trial-1)*NTones : trial*NTones;
    delays = TimeOrder(idx)';
    amps = AmpOrder(idx);
    ampsTDT = Amp_TDT(idx)';
    freqs = FreqOrder(idx)';  

    % set the first silence period and first (amp,freq)
    baseline_delay = delays(1);
    delays = [delays(2:end) baseline_delay];
    RP.SetTagVal('Baseline',baseline_delay*1000);
    RP.WriteTagV('Delays', 0, delays);
    
    ampsTDT = [ampsTDT(end), ampsTDT(1:end-1)];
    freqs = [freqs(end), freqs(1:end-1)];
    RP.WriteTagV('Amps', 0, ampsTDT);
    RP.WriteTagV('Freqs', 0, freqs);
    
    % waiting for trigger from other hardware
    disp('TDT is waiting for the trigger from A-M system')
    fprintf('\n');
    while ~RP.GetTagVal('Trigger')
    end
    trial
    
    while RP.GetTagVal('Count')<NTones
    end
    pause(2*Pulse_Dur*0.001);
end
pause(1);
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');

RP.Halt;
disp('Done Dude!')

%% Saving Stimulus Set Information

StimSet.Fspec = Fspec;
StimSet.Amp_Array = Amp_Array;
StimSet.ramp = ramp;
StimSet.Pulse_Dur = Pulse_Dur;
StimSet.ISI = ISI;
StimSet.jitter_percentage = jitter_percentage;
StimSet.Randomized = Randomized;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;

Ntones = length(FreqOrder)/Ntrials;
StimSet.distinctTones = Ntones;

for n = 1:tot_trials
	StimSet.trials(n).Freq = FreqOrder(n)/1000;
	StimSet.trials(n).Amp = AmpOrder(n);
	StimSet.trials(n).Time = TimeOrder(n);
end

fileno = str2double(fileno); fileno = round(fileno);

Amp_str = [];
for k = 1:length(Amp_Array)
	Amp_str = [Amp_str, num2str(Amp_Array(k)),'_'];
end
filename = ['0',num2str(fileno),'_fmin_',num2str(min_Freq),...
	'_fmax_',num2str(max_Freq),'_FreqPerOct_',num2str(Freq_per_Oct),...
	'_Amp_',Amp_str,'dur_',num2str(Pulse_Dur),'_ramp_',num2str(ramp),...
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


