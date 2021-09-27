function [FreqTone, AmpTone] = run_button_Callback_ABR_Tones(varargin)
clc
FolderPath = pwd;
TDTFileName = 'Ephus_ContinuousPlay_ABR_GUI.rcx';
TDTFilePath = sprintf('%s\\%s',FolderPath,TDTFileName);
TDT = TDTRP(TDTFilePath, 'RX6');
% TDT.halt;
min_Freq_h = findobj('tag','min_Freq');
max_Freq_h = findobj('tag','max_Freq');
Freq_per_Oct_h = findobj('tag','Freq_per_Oct');
Amp_Array_h = findobj('tag','Amp_Array');
ramp_h = findobj('tag','ramp');
Pulse_Dur_h = findobj('tag','Pulse_Dur');
ClickDur_h = findobj('tag','ClickDur');
Rate_h = findobj('tag','Rate');
Ntone_h = findobj('tag','Ntone');
Nclick_h = findobj('tag','Nclick');
Click_h = findobj('tag','Click');
Randomized_h = findobj('tag','Randomized');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');

% Get the specifications from GUI
min_Freq = get(min_Freq_h,'string');
min_Freq = str2double(min_Freq);
min_Freq(isnan(min_Freq)) = 3; % if the box is empty then this is default value.
max_Freq = get(max_Freq_h,'string');
max_Freq = str2double(max_Freq);
max_Freq(isnan(max_Freq)) = 96; % if the box is empty then this is default value.
Freq_per_Oct = get(Freq_per_Oct_h,'string'); %Freqency per Octave
Freq_per_Oct =  str2double(Freq_per_Oct);
Freq_per_Oct(isnan(Freq_per_Oct)) = 1; % if the box is empty then this is default value.
Fspec = [min_Freq, max_Freq, Freq_per_Oct];
Amp_Array = get(Amp_Array_h,'string');
Amp_Array = str2num(Amp_Array);
if isempty(Amp_Array)
    Amp_Array = 20:10:90; % if the box is empty then this is default value.
end
ramp = get(ramp_h,'string');
ramp = str2double(ramp);
ramp(isnan(ramp)) = 1.5; % if the box is empty then this is default value.
Pulse_Dur = get(Pulse_Dur_h,'string');
Pulse_Dur = str2double(Pulse_Dur);
Pulse_Dur(isnan(Pulse_Dur)) = 6; % if the box is empty then this is default value.
ClickDur = get(ClickDur_h,'string');
ClickDur = str2double(ClickDur);
ClickDur(isnan(ClickDur)) = 1; % if the box is empty then this is default value.
Rate = get(Rate_h,'string');
Rate = str2double(Rate);
Rate(isnan(Rate)) = 21; % if the box is empty then this is default value.
ISI = 1000/Rate; % ISI in msec
Ntone = get(Ntone_h,'string');
Ntone = str2double(Ntone);
Ntone(isnan(Ntone)) = 1024; % if the box is empty then this is default value.
Nclick = get(Nclick_h,'string');
Nclick = str2double(Nclick);
Nclick(isnan(Nclick)) = 512; % if the box is empty then this is default value.
Click = get(Click_h, 'value');
Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');

%% circuit info
% size of the entire serial buffer
npts = TDT.RP.GetTagSize('datain');
% npts = 200000;
% serial buffer will be divided into two buffers A & B
bufpts = npts/2;
Fs = TDT.FS; % TDT sampling frequency
% Fs = 195312.500000;

%% Generate click
Amp_Array_Click = min(Amp_Array):5:max(Amp_Array);
[C,ClickOrder,AmpClick,PhaseClick] = ClickAPorder(Amp_Array_Click,Nclick,Randomized);
AmpTDT_Click = getTDTAmpClick(AmpClick, SpeakerName);
N = round(ClickDur*Fs/1000); N = 2*round(N/2);
click = ones(N,1);
nC = length(C);
N_ISI = round(Fs*ISI/1000);
N_ISI = N_ISI-N;
N_ISI = 2*round(N_ISI/2);
click = [zeros(N_ISI/2,nC); click; zeros(N_ISI/2,nC)];
N_ISI = size(click,1);
NClicks = length(ClickOrder);
NsC = 2; % divide clicks to this number of slots. For each slot we need one iteration in ephus.
ClickrecDurN = NClicks*N_ISI/NsC;

%% Generate tones
[F,FreqTone,AmpTone,PhaseTone] = FATorder(Fspec,Amp_Array,Ntone,Randomized);
AmpTDT_Tone = getTDTAmp(F, FreqTone, AmpTone, SpeakerName);
FreqTone = FreqTone*1000;
F = F*1000;
N = round(Pulse_Dur*Fs/1000); N = 2*round(N/2);
nF = length(F);
tone = zeros(N,nF);
FreqIndex = FreqTone;
for f = 1:nF
    phi = 0;
    tone(:,f) = Cos2Gate(Fs,F(f),phi,Pulse_Dur,ramp);
    FreqIndex(FreqIndex==F(f)) = f;
end
N_ISI = round(Fs*ISI/1000);
N_ISI = N_ISI-N;
N_ISI = 2*round(N_ISI/2);
tone = [zeros(N_ISI/2,nF); tone; zeros(N_ISI/2,nF)];
N_ISI = size(tone,1);
NTones = length(FreqTone);
TonerecDurN = NTones*N_ISI/nF;

%% Play out clicks
start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
if Click
    fprintf('Click Stimuli\n');
    TDT.write('recDurN', ClickrecDurN);
    for n = 1:NsC
        ii = 1+(n-1)*NClicks/NsC:n*NClicks/NsC;
        % load up entire buffer with segments A and B
        n2=0;idx2=0;
        [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,npts,N_ISI);
        StimTDT = bufferStim(ClickOrder(ii),AmpTDT_Click(ii),PhaseClick(ii),click,n1,n2,idx1,idx2);
        TDT.write('datain', StimTDT);
        index = length(StimTDT);
        % waiting for trigger from other hardware
        fprintf('TDT is waiting for the trigger from A-M system\n')
        while ~TDT.read('Trigger')
        end
        TDT.trg(1);
        fprintf('Trigger received. Recording in progress...\n\n');
        curindex = TDT.read('index');
        while index < ClickrecDurN
            % segment A stim
            [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,bufpts,N_ISI);
            StimTDT = bufferStim(ClickOrder(ii),AmpTDT_Click(ii),PhaseClick(ii),click,n1,n2,idx1,idx2);
            % wait until done playing A
            while(curindex < bufpts)
                curindex = TDT.read('index');
                pause(.05);
            end
            % load the next signal segment
            TDT.write('datain', StimTDT);
            index = index + length(StimTDT);
            % checks to see if the data transfer rate is fast enough
            curindex = TDT.read('index');
            if(curindex < bufpts)
                warning('Transfer rate is too slow');
            end
            % segment B stim
            if index < ClickrecDurN
                [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,bufpts,N_ISI);
                StimTDT = bufferStim(ClickOrder(ii),AmpTDT_Click(ii),PhaseClick(ii),click,n1,n2,idx1,idx2);
            end
            % wait until start playing A
            while(curindex > bufpts)
                curindex = TDT.read('index');
                pause(.05);
            end
            % load segment B
            if index >= ClickrecDurN
                break;
            end
            TDT.write('datain', StimTDT, 'OFFSET', bufpts);
            index = index + length(StimTDT);
            % make sure we're still playing A
            curindex = TDT.read('index');
            if(curindex > bufpts)
                warning('Transfer rate too slow');
            end
        end
        pause(npts/Fs); % make sure that buffer is palyed out
        curindex = TDT.read('index');
        fprintf(['Current index: ' num2str(curindex),'\n']);
    end
end

%% Play out tones
fprintf('Tone Stimuli\n');
TDT.write('recDurN', TonerecDurN);
for f = 1:nF
    ii = 1+(f-1)*NTones/nF:f*NTones/nF;
    % load up entire buffer with segments A and B
    n2=0;idx2=0;
    [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,npts,N_ISI);
    StimTDT = bufferStim(FreqIndex(ii),AmpTDT_Tone(ii),PhaseTone(ii),tone,n1,n2,idx1,idx2);
    TDT.write('datain', StimTDT);
    index = length(StimTDT);
    % waiting for trigger from other hardware
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~TDT.read('Trigger')
    end
    TDT.trg(1);
    fprintf('Trigger received. Recording in progress...\n\n');
    curindex = TDT.read('index');
    while index < TonerecDurN
        % segment A stim
        [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,bufpts,N_ISI);
        StimTDT = bufferStim(FreqIndex(ii),AmpTDT_Tone(ii),PhaseTone(ii),tone,n1,n2,idx1,idx2);
        % wait until done playing A
        while(curindex < bufpts)
            curindex = TDT.read('index');
            pause(.05);
        end
        % load the next signal segment
        TDT.write('datain', StimTDT);
        index = index + length(StimTDT);
        % checks to see if the data transfer rate is fast enough
        curindex = TDT.read('index');
        if(curindex < bufpts)
            warning('Transfer rate is too slow');
        end
        % segment B stim
        if index < TonerecDurN
            [n1,n2,idx1,idx2] = bufferInfo(n2,idx2,bufpts,N_ISI);
            StimTDT = bufferStim(FreqIndex(ii),AmpTDT_Tone(ii),PhaseTone(ii),tone,n1,n2,idx1,idx2);
        end
        % wait until start playing A
        while(curindex > bufpts)
            curindex = TDT.read('index');
            pause(.05);
        end
        % load segment B
        if index >= TonerecDurN
            break;
        end
        TDT.write('datain', StimTDT, 'OFFSET', bufpts);
        index = index + length(StimTDT);
        % make sure we're still playing A
        curindex = TDT.read('index');
        if(curindex > bufpts)
            warning('Transfer rate too slow');
        end
    end
    pause(npts/Fs); % make sure that buffer is palyed out
    curindex = TDT.read('index');
    fprintf(['Current index: ' num2str(curindex),'\n']);
end

%% stop playing
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
TDT.halt;
disp('Done Dude!')

%% Saving Stimulus Set Information
StimSet.StimType = 'ABR';
StimSet.Tone = 1;
StimSet.Fspec = Fspec;
StimSet.Amp_Array = Amp_Array;
StimSet.Amp_Array_Click = Amp_Array_Click;
StimSet.distinctTones = NTones/Ntone;
StimSet.Pulse_Dur = Pulse_Dur;
StimSet.ramp = ramp;
StimSet.Ntone = Ntone;
StimSet.ISI = ISI;
StimSet.Rate = Rate;
StimSet.Click = Click;
StimSet.distinctTones = NClicks/Nclick;
StimSet.ClickDur = ClickDur;
StimSet.Nclick = Nclick;
StimSet.NsC = NsC;
StimSet.Randomized = Randomized;
StimSet.start_date = start_date;
StimSet.stop_date = stop_date;
StimSet.Notes = get(Notes_h,'string');
StimSet.SpeakerName = SpeakerName;

k = 0;
if Click
    for n = 1:NClicks
        k = k+1;
        StimSet.trials(k).Freq = Inf;
        StimSet.trials(k).Amp = AmpClick(n);
        StimSet.trials(k).Phase = PhaseClick(n);
    end
end
for n = 1:NTones
    k = k+1;
    StimSet.trials(k).Freq = FreqTone(n)/1000;
    StimSet.trials(k).Amp = AmpTone(n);
    StimSet.trials(k).Phase = PhaseTone(n);
end
fileno = str2double(fileno); fileno = round(fileno);
Amp_str = sprintf('%d_%d_%d',Amp_Array(1),Amp_Array(end),Amp_Array(2)-Amp_Array(1));

filename = ['0',num2str(fileno),'_ABR','_F_',num2str(min_Freq),...
    '_',num2str(max_Freq),'_',num2str(Freq_per_Oct),'_Amp_',Amp_str,...
    '_Click_',num2str(Click)];

if fileno <= 9.5
    eval(['StimSet_0',num2str(fileno), '=', 'StimSet']);
    save([filename,'.mat'], 'StimSet');
    
    if fileno >= 8.5
        set(fileno_h,'string',num2str(fileno+1));
    else
        set(fileno_h,'string',['0',num2str(fileno+1)]);
    end
else
    eval(['StimSet_',num2str(fileno), '=', 'StimSet']);
    save([filename(2:end),'.mat'],'StimSet');
    set(fileno_h,'string',num2str(fileno+1));
end


