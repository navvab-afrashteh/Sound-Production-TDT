function [FreqOrder, AmpOrder, TimeOrder] = run_button_Callback_Tinnitus_Tones(varargin)
clc
FolderPath = pwd;
TDTFileName = 'Ephus_PureTones_Tinnitus_GUI.rcx';
TDTFilePath = sprintf('%s\\%s',FolderPath,TDTFileName);
RP = TDTRP(TDTFilePath, 'RX6');
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
nFrames_h = findobj('tag','nFrames');
NCam_h = findobj('tag','NCam');
nSponB_h = findobj('tag','nSponB');
nSponE_h = findobj('tag','nSponE');
% Get the specifications from GUI
NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 7; % if the box is empty then this is default value.
nSponB = get(nSponB_h,'string');
nSponB = str2double(nSponB);
nSponB(isnan(nSponB)) = 2; % if the box is empty then this is default value.
nSponE = get(nSponE_h,'string');
nSponE = str2double(nSponE);
nSponE(isnan(nSponE)) = 0; % if the box is empty then this is default value.
min_Freq = get(min_Freq_h,'string');
min_Freq = str2double(min_Freq);
min_Freq(isnan(min_Freq)) = 3; % if the box is empty then this is default value.
max_Freq = get(max_Freq_h,'string');
max_Freq = str2double(max_Freq);
max_Freq(isnan(max_Freq)) = 48; % if the box is empty then this is default value.
Freq_per_Oct = get(Freq_per_Oct_h,'string'); %Freqency per Octave
Freq_per_Oct =  str2double(Freq_per_Oct);
Freq_per_Oct(isnan(Freq_per_Oct)) = 1; % if the box is empty then this is default value.
Fspec = [min_Freq, max_Freq, Freq_per_Oct];
Amp_Array = get(Amp_Array_h,'string');
Amp_Array = str2num(Amp_Array);
if isempty(Amp_Array)
	Amp_Array = 20:15:80; % if the box is empty then this is default value.
end
ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 5; % if the box is empty then this is default value.
Pulse_Dur = get(Pulse_Dur_h,'string');
Pulse_Dur = str2double(Pulse_Dur);
Pulse_Dur(isnan(Pulse_Dur)) = 50; % if the box is empty then this is default value.
ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 3000; % if the box is empty then this is default value.
jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;
Ntrials = get(Ntrials_h,'string');
Ntrials = str2double(Ntrials);
Ntrials(isnan(Ntrials)) = 16; % if the box is empty then this is default value.
Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);
nEvoked = NCam-nSponB-nSponE;
[F, FreqOrder, AmpOrder, TimeOrder] = FATorder(Fspec, Amp_Array,...
	ISI, jitter,Ntrials, Randomized, FsCam,nFrames,nEvoked);
Amp_TDT = getTDTAmp(F, FreqOrder, AmpOrder, SpeakerName);
FreqOrder = FreqOrder*1000;
RP.write('PulseDur', Pulse_Dur-ramp);
tot_trials = size(FreqOrder,1);
NTones = tot_trials/nEvoked;
t1(1:tot_trials,1) = 0;

[RecNames{1:nEvoked}] = deal('Tone');
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

%%
for nc = 1:nEvoked
    % waiting for trigger from other hardware
    disp('TDT is waiting for the trigger from A-M system')
    fprintf('\n');
    while ~RP.read('Trigger')
    end
    idx = 1+(nc-1)*NTones : nc*NTones;
    delays = TimeOrder(idx)';
    ampOrder = AmpOrder(idx);
    ampsTDT = Amp_TDT(idx)';
    freqs = FreqOrder(idx)';  
    tplay = 0;
    totCam = tic;
    for tone = 1:NTones
        Ttrial = tic;
        % pause between tones
        pause(delays(tone)-tplay)
        % loading amplitude and ferquency to TDT
        Tplay = tic;
        a = ampsTDT(tone);
        f = freqs(tone);
        if a > 9.75
            a = 9.75;
            disp(['Clipping @',num2str(f/1000), ' kHz; ', num2str(ampOrder(tone)), ' dB'])
        end
        RP.write('Amp',a);
        RP.write('Freq',f);
        % playing the tone
        RP.trg(1);
        pause((3*ramp + Pulse_Dur)/1000)
        disp([num2str(f/1000), ' kHz, ',...
            num2str(ampOrder(tone)), ' dB, ',...
            num2str(tone+(nc-1)*NTones),'/',num2str(tot_trials)])
        fprintf('\n');
        tplay = toc(Tplay);
        t1((nc-1)*NTones +tone) = toc(Ttrial); % This is true time order
    end
    TCam = toc(totCam);
    disp(['Total time: ',num2str(TCam), 's']);
    fprintf('\n');
end
pause(ISI/1000);

%% End spontanious recordings
for spE = 1:nSponE
    recDur = nFrames/FsCam;
    % waiting for trigger from other hardware
    fprintf('\n Spontaneous Recording:\n');
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~TDT.read('Trigger')
    end
    start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fprintf('Trigger received. Recording in progress...\n\n');
    pause(recDur);
end
%% stop playing
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
TimeOrder = t1;
RP.halt;
disp('Done Dude!')
%% Saving Stimulus Set Information
StimSet.StimType    = 'Tone';
StimSet.NCam = NCam;
StimSet.RecNames = RecNames;
StimSet.nEvoked = nEvoked;
StimSet.nSponB = nSponB;
StimSet.nSponE = nSponE;
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
StimSet.distinctTones = tot_trials/Ntrials;
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
if nSponB>0 && nSponE>0
    filename = ['0',num2str(fileno),'_Tone_',num2str(nSponB),'Spon_',...
        num2str(nEvoked),'Evoked_',num2str(nSponE),'Spon','_fmin_',num2str(min_Freq),...
        '_fmax_',num2str(max_Freq),'_FreqPerOct_',num2str(Freq_per_Oct),...
        '_Amp_',Amp_str,'dur_',num2str(Pulse_Dur),'_ramp_',num2str(ramp),...
        '_ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
        '_Ntrials_',num2str(NCam)];
elseif nSponB>0 && nSponE==0
    filename = ['0',num2str(fileno),'_Tone_',num2str(nSponB),'Spon_',...
        num2str(nEvoked),'Evoked','_fmin_',num2str(min_Freq),...
        '_fmax_',num2str(max_Freq),'_FreqPerOct_',num2str(Freq_per_Oct),...
        '_Amp_',Amp_str,'dur_',num2str(Pulse_Dur),'_ramp_',num2str(ramp),...
        '_ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
        '_Ntrials_',num2str(NCam)];
elseif nSponB==0 && nSponE>0
    filename = ['0',num2str(fileno),'_Tone_',...
        num2str(nEvoked),'Evoked_',num2str(nSponE),'Spon','_fmin_',num2str(min_Freq),...
        '_fmax_',num2str(max_Freq),'_FreqPerOct_',num2str(Freq_per_Oct),...
        '_Amp_',Amp_str,'dur_',num2str(Pulse_Dur),'_ramp_',num2str(ramp),...
        '_ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
        '_Ntrials_',num2str(NCam)];
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


