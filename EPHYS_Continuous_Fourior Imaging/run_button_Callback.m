function [FreqOrder, AmpOrder] = run_button_Callback(varargin)
clc
RP = TDTRP('Ephus_FI_GUI.rcx', 'RX6');
% RP.Halt;

min_Freq_h = findobj('tag','min_Freq');
max_Freq_h = findobj('tag','max_Freq');
Amp_h = findobj('tag','Amp');
ramp_h = findobj('tag','ramp');
Pulse_Dur_h = findobj('tag','Pulse_Dur');
ISI_h = findobj('tag','ISI');
sweep_Dur_h = findobj('tag','sweep_Dur');
base_Dur_h = findobj('tag','base_Dur');
Mode_h = findobj('tag','Mode');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');
FsCam_h = findobj('tag','FsCam');
nFrames_h = findobj('tag','nFrames');

% Get the specifications from GUI
min_Freq = get(min_Freq_h,'string');
min_Freq = str2double(min_Freq);
min_Freq(isnan(min_Freq)) = 2; % if the box is empty then this is default value.

max_Freq = get(max_Freq_h,'string');
max_Freq = str2double(max_Freq);
max_Freq(isnan(max_Freq)) = 64; % if the box is empty then this is default value.

sweep_Dur = get(sweep_Dur_h,'string'); %sweep Duration
sweep_Dur =  str2double(sweep_Dur);
sweep_Dur(isnan(sweep_Dur)) = 8000; % if the box is empty then this is default value.

base_Dur = get(base_Dur_h,'string'); %sweep Duration
base_Dur =  str2double(base_Dur);
base_Dur(isnan(base_Dur)) = 0; % if the box is empty then this is default value.


Amp = get(Amp_h,'string');
Amp = str2double(Amp);
Amp(isnan(Amp)) = 80; % if the box is empty then this is default value.

ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 5; % if the box is empty then this is default value.

Pulse_Dur = get(Pulse_Dur_h,'string');
Pulse_Dur = str2double(Pulse_Dur);
Pulse_Dur(isnan(Pulse_Dur)) = 50; % if the box is empty then this is default value.

ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 250; % if the box is empty then this is default value.

Mode = get(Mode_h,'string');
Mode = Mode{get(Mode_h,'value')};

fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);

Rec_Dur = 1000*nFrames/FsCam; % total recording duration by the camera

[F, Fspec, FreqOrder, AmpOrder, startDur, NTrials, NTones] = FI_FATotder(Amp, min_Freq, max_Freq, base_Dur, ISI, sweep_Dur, Mode, Rec_Dur);
Amp_TDT = getTDTAmp(F, FreqOrder, AmpOrder, SpeakerName);
FreqOrder = FreqOrder*1000;
Amp_TDT = [Amp_TDT(end), Amp_TDT(1:end-1)];
FreqOrder = [FreqOrder(end), FreqOrder(1:end-1)];

RP.SetTagVal('NTones', NTones);
RP.SetTagVal('NTrials', NTrials);
FsTDT = RP.GetSFreq();
sweepSamp = ceil(FsTDT*(sweep_Dur+base_Dur)/1000);
RP.SetTagVal('sweepSamp', sweepSamp);
RP.SetTagVal('PulseDur', Pulse_Dur);
RP.SetTagVal('ramp', ramp);
Tsilence = ISI-Pulse_Dur;
RP.SetTagVal('Tsilence', Tsilence);
RP.SetTagVal('startDur', startDur);
RP.WriteTagV('Amps', 0, Amp_TDT);
RP.WriteTagV('Freqs', 0, FreqOrder);
    
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
% waiting for trigger from other hardware
disp('TDT is waiting for the trigger from A-M system')
fprintf('\n');

while ~RP.GetTagVal('Trigger')
end

trial = 0;
playing = 1;
while playing
     trialNum = RP.GetTagVal('trialNum');
    if ~(trial-1 == trialNum)
        trial = trial+1;
        disp(['trial ', num2str(trial), ' of ', num2str(NTrials)])
        fprintf('\n');
        pause((sweep_Dur+base_Dur)/1000)
    end
    if (trial == NTrials)
        playing = 0;
    end
end

stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
RP.Halt;
disp('Done Dude!')

%% Saving Stimulus Set Information
StimSet.Fspec = Fspec;
StimSet.Amp_Array = Amp;
StimSet.ramp = ramp;
StimSet.Pulse_Dur = Pulse_Dur;
StimSet.ISI = ISI;
StimSet.sweep_Dur = sweep_Dur;
StimSet.base_Dur = base_Dur;
StimSet.Rec_Dur = Rec_Dur;
StimSet.Mode = Mode;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;
StimSet.distinctTones = NTones;

fileno = str2double(fileno); fileno = round(fileno);

Amp_str = [];
for k = 1:length(Amp)
	Amp_str = [Amp_str, num2str(Amp(k)),'_'];
end
sweep_Dur = round(sweep_Dur/1000);
base_Dur = round(base_Dur/1000);
Rec_Dur = round(Rec_Dur/1000);
filename = ['0',num2str(fileno),'_',Mode,'_fmin_',num2str(min_Freq),...
	'_fmax_',num2str(max_Freq),...
	'_Amp_',Amp_str,'dur_',num2str(Pulse_Dur),'_ramp_',num2str(ramp),...
	'_ISI_',num2str(ISI),'_sweepDur_',num2str(sweep_Dur),...
	'_baseDur_', num2str(base_Dur), '_RecDur_',num2str(Rec_Dur)];
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


