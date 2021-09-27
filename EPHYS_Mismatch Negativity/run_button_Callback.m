function run_button_Callback(varargin)
clc
FolderPath = pwd;
TDTFileName = 'Ephus_PureTones_MMN_GUI.rcx';
TDTFilePath = sprintf('%s\\%s',FolderPath,TDTFileName);
RP = TDTRP(TDTFilePath, 'RX6');

MMPercent_h = findobj('tag','MMPercent');
FI_h = findobj('tag','FI');
FII_h = findobj('tag','FII');
FC_h = findobj('tag','FC');
Amp_Array_h = findobj('tag','Amp_Array');
ramp_h = findobj('tag','ramp');
Pulse_Dur_h = findobj('tag','Pulse_Dur');
ISI_h = findobj('tag','ISI');
jitter_percentage_h = findobj('tag','jitter_percentage');
NperCam_h = findobj('tag','NperCam');
Randomized_h = findobj('tag','Randomized');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');
FsCam_h = findobj('tag','FsCam');
nFrames_h = findobj('tag','nFrames');
NCam_h = findobj('tag','NCam');
nSponB_h = findobj('tag','nSponB');
nSponE_h = findobj('tag','nSponE');
seed_h = findobj('tag','seed');
% Get the specifications from GUI
MMPercent = get(MMPercent_h,'string');
MMPercent = str2num(MMPercent);
MMPercent(isnan(MMPercent)) = 12.5; % if the box is empty then this is default value.
FI = get(FI_h,'string');
FI = str2double(FI);
FI(isnan(FI)) = 3.6; % if the box is empty then this is default value.
FII = get(FII_h,'string');
FII = str2double(FII);
FII(isnan(FII)) = 5.1; % if the box is empty then this is default value.
FC = get(FC_h,'string'); %Freqency per Octave
FC =  str2num(FC);
if isempty(FC)
	FC = [3.3:0.3:5.4]; % if the box is empty then this is default value.
end
AmpArray = get(Amp_Array_h,'string');
AmpArray = str2num(AmpArray);
if isempty(AmpArray)
	AmpArray = 80; % if the box is empty then this is default value.
end
ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 5; % if the box is empty then this is default value.
PulseDur = get(Pulse_Dur_h,'string');
PulseDur = str2double(PulseDur);
PulseDur(isnan(PulseDur)) = 100; % if the box is empty then this is default value.
ISI = get(ISI_h,'string');
ISI = str2double(ISI);
ISI(isnan(ISI)) = 1000; % if the box is empty then this is default value.
jitter_percentage = get(jitter_percentage_h,'string');
jitter_percentage = str2double(jitter_percentage);
jitter_percentage(isnan(jitter_percentage)) = 0; % if the box is empty then this is default value.
jitter = (jitter_percentage/100)*ISI;
NperCam = get(NperCam_h,'string');
NperCam = str2double(NperCam);
NperCam(isnan(NperCam)) = 400; % if the box is empty then this is default value.
nSponB = get(nSponB_h,'string');
nSponB = str2double(nSponB);
nSponB(isnan(nSponB)) = 0; % if the box is empty then this is default value.
nSponE = get(nSponE_h,'string');
nSponE = str2double(nSponE);
nSponE(isnan(nSponE)) = 0; % if the box is empty then this is default value.
seed = get(seed_h,'string');
seed = str2double(seed);
seed(isnan(seed)) = 1; % if the box is empty then this is default value.
Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 3; % if the box is empty then this is default value.
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RP.write('PulseDur', PulseDur-ramp);
minDist = 3;
recDur = nFrames/FsCam;
[TimeOrderAB, FreqOrderAB, idxMMN, FmainAB] = blockAB(FI,FII,minDist, MMPercent, ...
    recDur, PulseDur/1000, ISI/1000, jitter,NperCam,seed);
AmpOrderAB = AmpArray*ones(size(FreqOrderAB));
Amp_TDT_AB = getTDTAmp(FmainAB, FreqOrderAB, AmpOrderAB, SpeakerName);
FreqOrderAB = FreqOrderAB*1000;
[TimeOrderC, FreqOrderC] = blockC(FC, recDur, PulseDur/1000, ISI/1000, jitter,NperCam,seed);
AmpOrderC = AmpArray*ones(size(FreqOrderC));
Amp_TDT_C = getTDTAmp(FC, FreqOrderC, AmpOrderC, SpeakerName);
FreqOrderC = FreqOrderC*1000;

FreqOrder = [FreqOrderC,FreqOrderAB]; FreqOrder = round(FreqOrder);
AmpOrder = [AmpOrderC,AmpOrderAB];
Amp_TDT = [Amp_TDT_C,Amp_TDT_AB];
TimeOrder = [TimeOrderC,TimeOrderAB];

nEvoked = NCam-nSponB-nSponE;
t1 = zeros(NperCam,nEvoked);

[RecNames{1:nEvoked}] = deal('Tone');
[SponBNames{1:nSponB}] = deal('Spon');
RecNames = [SponBNames, RecNames];
[SponENames{1:nSponE}] = deal('Spon');
RecNames = [RecNames, SponENames];

%% First spontaneous recording
for spB = 1:nSponB
    % waiting for trigger from other hardware
    fprintf('\nSpontaneous Recording:\n');
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~RP.read('Trigger')
    end
    start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fprintf('Trigger received. Recording in progress...\n\n');
    pause(recDur);
end

%% Block C and then A and B (A and B are randomized)
for nc = 1:nEvoked
    delays = TimeOrder(:,nc);
    ampOrder = AmpOrder(:,nc);
    ampsTDT = Amp_TDT(:,nc);
    freqs = FreqOrder(:,nc);  
    tplay = 0;
    % waiting for trigger from other hardware
    disp('TDT is waiting for the trigger from A-M system')
    fprintf('\n');
    while ~RP.read('Trigger')
    end
    if nSponB == 0
        start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    end
    totCam = tic;
    for tone = 1:NperCam
        Ttrial = tic;
        % pause between tones
        pause(delays(tone)-PulseDur/1000)
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
        pause((ramp + PulseDur)/1000)
        disp([num2str(f/1000), ' kHz, ',...
            num2str(ampOrder(tone)), ' dB, ',...
            num2str(tone+(nc-1)*NperCam),'/',num2str(NperCam*nEvoked)])
        fprintf('\n');
        tplay = toc(Tplay);
        t1(tone,nc) = toc(Ttrial); % This is true time order
    end
    TCam = toc(totCam);
    disp(['Total time: ',num2str(TCam), 's']);
    fprintf('\n');
end
pause(ISI/1000);

%% End spontanious recordings
for spE = 1:nSponE
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
StimSet.StimType    = 'MMN Tone';
StimSet.NCam = NCam;
StimSet.RecNames = RecNames;
StimSet.nEvoked = nEvoked;
StimSet.nSponB = nSponB;
StimSet.nSponE = nSponE;
StimSet.Freqs = unique(FreqOrder);
StimSet.FreqOrderAB = FreqOrderAB;
StimSet.FreqOrderC = FreqOrderC;
StimSet.idxMMN = idxMMN;
StimSet.FmainAB = FmainAB;
StimSet.Amp_Array = AmpArray;
StimSet.ramp = ramp;
StimSet.Pulse_Dur = PulseDur;
StimSet.ISI = ISI;
StimSet.jitter_percentage = jitter_percentage;
StimSet.Randomized = Randomized;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam = FsCam;

StimSet.distinctTones = length(unique(FreqOrder))*length(unique(AmpOrder));
FreqOrder = FreqOrder(:);
AmpOrder = AmpOrder(:);
TimeOrder = TimeOrder(:);
Ntot = length(FreqOrder);
for n = 1:Ntot
	StimSet.trials(n).Freq = FreqOrder(n)/1000;
	StimSet.trials(n).Amp = AmpOrder(n);
	StimSet.trials(n).Time = TimeOrder(n);
end
fileno = str2double(fileno); fileno = round(fileno);
Amp_str = [];
for k = 1:length(AmpArray)
	Amp_str = [Amp_str, num2str(AmpArray(k)),'_'];
end
filename = ['0',num2str(fileno),'_MMNtone_',num2str(nSponB),'Spon_',...
    num2str(nEvoked),'Evoked_',num2str(nSponE),'Spon','_FI_',num2str(FI),...
    '_FII_',num2str(FII),'_Amp_',Amp_str,'dur_',...
    num2str(PulseDur),'_ramp_',num2str(ramp),...
    '_ISI_',num2str(ISI),'_jitter_',num2str(round(jitter)),...
    '_Ntrials_',num2str(NCam)];

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


