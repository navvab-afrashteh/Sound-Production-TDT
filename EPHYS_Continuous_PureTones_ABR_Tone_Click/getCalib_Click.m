function [calib_V, calib_dB] = getCalib_Click(SpeakerName,varargin)

SpeakerFolderName = 'Speaker Calibration';
PathName = pwd;
idx = strfind(PathName,'TDT_Navvab');
PathName = PathName(1:idx-1);
PathName = sprintf('%s\\%s',PathName,SpeakerFolderName);
FolderName = sprintf('%s\\%s',PathName,SpeakerName);
FileName = 'calib_click.csv';
FilePath = sprintf('%s\\%s',FolderName,FileName);

fid = fopen(FilePath, 'rt');
x = fread(fid, inf, '*char');
fclose(fid);

x = x(:)';
lines = regexp(x, newline, 'split');

calib_V = [];
calib_dB = [];
for i = 2:length(lines)-1
    parts = regexp(lines{i}, char(44), 'split');
    calib_V = [calib_V str2double(parts{1})];
    calib_dB = [calib_dB str2double(parts{2})];
end

