function [calib_V, calib_level, calib_freq, calib_amp, varargout] = getCalib(SpeakerName,varargin)
Fs = NaN;
if nargin>1
    Fs = varargin{1};
end
SpeakerFolderName = 'Speaker Calibration';
PathName = pwd;
idx = strfind(PathName,'TDT_Navvab');
PathName = PathName(1:idx-1);
PathName = sprintf('%s\\%s',PathName,SpeakerFolderName);
FolderName = sprintf('%s\\%s',PathName,SpeakerName);
FileName = 'calib.csv';
FilePath = sprintf('%s\\%s',FolderName,FileName);

filterFileName = 'FIR';
outputFilePath = sprintf('%s\\%s.txt',FolderName,filterFileName);


fid = fopen(FilePath, 'rt');
x = fread(fid, inf, '*char');
fclose(fid);

x = x(:)';
lines = regexp(x, char(10), 'split');
l1 = lines{1};

parts = regexp(l1, char(44), 'split');
test_signal = parts{end};
temp = str2double(regexp(test_signal, '\d+.\d+', 'match'));

calib_V = temp(1); % calibration voltage level
calib_level = temp(2); % calibration level in dB

calib_freq = [];
calib_amp = [];
for i = 2:length(lines)-1
    parts = regexp(lines{i}, char(44), 'split');
    calib_freq = [calib_freq str2double(parts{1})];
    calib_amp = [calib_amp str2double(parts{6})];
end

if nargout>4
    if ~isnan(Fs)
        filtDur = 1;% 1ms
        ntaps = round(Fs*filtDur/1000);
        ntaps = 2^nextpow2(ntaps);
        calib_freq = calib_freq/(Fs/2);
        % make sure beginning is zero
        calib_freq = [0 calib_freq];
        calib_amp = [calib_amp(1) calib_amp];
        % make sure end is 1;
        calib_freq = [calib_freq 1];
        calib_amp = [calib_amp calib_amp(end)];
        
        FIRcoeff = fir2(ntaps,calib_freq,10.^(calib_amp/20));
        
        fid = fopen(outputFilePath, 'wt+');
        fprintf(fid, '%6f\n', FIRcoeff);
        fclose(fid);

    else
        FIRcoeff = [];
    end
    varargout{1} = FIRcoeff;
end
