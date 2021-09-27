function run_button_Callback(varargin)

clc;
seed_h = findobj('tag','seed');
DRCs_h = findobj('tag','DRCs');
DMRs_h = findobj('tag','DMRs');
FM_h = findobj('tag','FM');
Vocal_h = findobj('tag','Vocal');
Randomized_h = findobj('tag','Randomized');
fileno_h = findobj('tag','fileno');
Notes_h = findobj('tag','Notes');
Speaker_h = findobj('tag','Speaker');
FsCam_h = findobj('tag','FsCam');
nFrames_h = findobj('tag','nFrames');
NCam_h = findobj('tag','NCam');
% Get the specifications from GUI
seed = get(seed_h,'string');
seed = str2double(seed);
seed(isnan(seed)) = 1; % if the box is empty then this is default value.
NCam = get(NCam_h,'string');
NCam = str2double(NCam);
NCam(isnan(NCam)) = 5; % if the box is empty then this is default value.
Randomized = get(Randomized_h, 'value');
fileno = get(fileno_h,'string');
SpeakerName = get(Speaker_h,'string');
FsCam = get(FsCam_h,'string'); FsCam = str2double(FsCam);
nFrames = get(nFrames_h,'string'); nFrames = str2double(nFrames);
RecNames = {};
if get(DRCs_h, 'value')
    RecNames = [RecNames, {'DRCs'}];
end
if get(DMRs_h, 'value')
    RecNames = [RecNames, {'DMRs'}];
end
if get(FM_h, 'value')
    RecNames = [RecNames, {'FM Sweeps'}];
end
if get(Vocal_h, 'value')
    RecNames = [RecNames, {'Vocal_STRFs'}];
end
NStim = length(RecNames);
if Randomized
    idx = randperm(NStim);
    RecNames = RecNames(idx);
end
if NStim > NCam
    RecNames = RecNames(1:NCam);
end
StimNames = RecNames;
NStim = length(RecNames);

NSpon = NCam-NStim;
NSponB = ceil(NSpon/2);
NSponE = NSpon - NSponB;
[SponBNames{1:NSponB}] = deal('Spon');
RecNames = [SponBNames, RecNames];
[SponENames{1:NSponE}] = deal('Spon');
RecNames = [RecNames, SponENames];

%% Load TDT protocol and get params
FolderPath = pwd;
TDTFileName = 'Continuous_Play_Modified.rcx';
TDTFilePath = sprintf('%s\\%s',FolderPath,TDTFileName);
TDT = TDTRP(TDTFilePath, 'RX6');
% size of the entire serial buffer
npts = TDT.RP.GetTagSize('datain');
% serial buffer will be divided into two buffers A & B
bufpts = npts/2;
Fs = TDT.FS; % TDT sampling frequency
Spec = cell(1,NCam);
%% First spontaneous recording
s = 0;
for spB = 1:NSponB
    s = s+1;
    Spec{s}.FsCam = FsCam;
    Spec{s}.recDur = nFrames/FsCam;
    Spec{s}.recType = RecNames{spB};
    % waiting for trigger from other hardware
    fprintf('\nSpontaneous Recording:\n');
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~TDT.read('Trigger')
    end
    start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fprintf('Trigger received. Recording in progress...\n\n');
    pause(Spec{s}.recDur);
end

%% Stim sequence names
for stimN = NSponB+1:NSponB+NStim
    fprintf('%s STRF:\n',RecNames{stimN});
    % load stim sequence
    sDataDir = sprintf('%s\\%s\\seed#%d',FolderPath,RecNames{stimN},seed);
    matFile = sprintf('%s\\Stim.mat',sDataDir);
    load(matFile,'StimTDT','spec');
    recDurN = length(StimTDT);
    TDT.write('recDurN', recDurN);
    % load up entire buffer with segments A and B
    index = 1;
    top = min(index+npts, recDurN);
    if top == recDurN
        top = top+1;
    end
    TDT.write('datain', StimTDT(index:top-1));
    index = top;
    % waiting for trigger from other hardware
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~TDT.read('Trigger')
    end
    fprintf('Trigger received. Recording in progress...\n\n');
    curindex = TDT.read('index');
    while index < recDurN
        % wait until done playing A
        while(curindex < bufpts)
            curindex = TDT.read('index');
            pause(.05);
        end
        % load the next signal segment
        top = min(index+bufpts, recDurN);
        if top == recDurN
            top = top+1;
        end
        TDT.write('datain', StimTDT(index:top-1));
        index = top;
        % checks to see if the data transfer rate is fast enough
        curindex = TDT.read('index');
        if(curindex < bufpts)
            warning('Transfer rate is too slow');
        end
        % wait until start playing A
        while(curindex > bufpts)
            curindex = TDT.read('index');
            pause(.05);
        end
        % load segment B
        if index > recDurN
            break;
        end
        top = min(index+bufpts, recDurN);
        if top == recDurN
            top = top+1;
        end
        TDT.write('datain', StimTDT(index:top-1), 'OFFSET', bufpts);
        index = top;
        % make sure we're still playing A
        curindex = TDT.read('index');
        if(curindex > bufpts)
            warning('Transfer rate too slow');
        end
    end
    pause(npts/Fs); % make sure that buffer is palyed out
    curindex = TDT.read('index');
    fprintf(['Current index: ' num2str(curindex),'\n']);
    % save specs for Stim
    s = s+1;
    Spec{s} = spec;
    Spec{s}.recType = RecNames{stimN};
end
%% End spontanious recordings
for spE = 1:NSponE
    s = s+1;
    Spec{s}.FsCam = FsCam;
    Spec{s}.recDur = nFrames/FsCam;
    Spec{s}.recType = RecNames{spE};
    % waiting for trigger from other hardware
    fprintf('\n Spontaneous Recording:\n');
    fprintf('TDT is waiting for the trigger from A-M system\n')
    while ~TDT.read('Trigger')
    end
    start_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fprintf('Trigger received. Recording in progress...\n\n');
    pause(Spec{s}.recDur);
end
% stop playing
stop_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
TDT.halt;
fprintf('Done Dude!\n')

%% Saving Stimulus Set Information
StimSet.StimType    = 'STRFs';
StimSet.StimNames   = StimNames;
StimSet.RecNames    = RecNames;
StimSet.NStim       = NStim;
StimSet.NSponB      = NSponB;
StimSet.NSponE      = NSponE;
StimSet.Spec        = Spec;
StimSet.SpeakerName = SpeakerName;
StimSet.FsCam       = FsCam;
StimSet.FsTDT       = Fs;
StimSet.nFrames     = nFrames;
StimSet.NCam        = NCam;
StimSet.seed        = seed;
StimSet.Randomized  = Randomized;
StimSet.Notes       = get(Notes_h,'string');
StimSet.start_date  = start_date;
StimSet.stop_date   = stop_date;

% generate file name and save
allStim = [];
for idx = 1:length(StimNames)
    temp = StimNames{idx};
    if temp(end) == 's'
        temp = temp(1:end-1);
    end
    idxu = strfind(StimNames{idx},'_');
    if ~isempty(idxu)
        temp = temp(1:idxu(1)-1);
    end
    idxs = strfind(StimNames{idx},' ');
    if ~isempty(idxs)
        temp = temp(1:idxs(1)-1);
    end
    allStim = [allStim,'_',temp];
end
if NSponB
    allStim = ['_Spon', allStim];
end
if NSponE
    allStim = [allStim,'_Spon'];
end
fileno = str2double(fileno); fileno = round(fileno);
filename = ['0',num2str(fileno),'_STRFs', allStim,'_seed',num2str(seed)];
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


