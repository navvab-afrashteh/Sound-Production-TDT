% One-channel continuous play example using a serial buffer
% This program writes to the buffer once it has cyled halfway through
% (double-buffering)

close all; 
clear all; clc;

RP = TDTRP('C:\Users\majidlab\Desktop\Navvab\TDT_Navvab\Simultaneous Continuous Play and Continuous Acquire\Simultaneous_Continuous_Play_Continuous_Acquire.rcx', 'RX6');

% size of the entire serial buffer
npts = RP.GetTagSize('datain');  

% serial buffer will be divided into two buffers A & B
bufpts = npts/2; 

% get device sampling frequency
fs = RP.GetSFreq();

% create one 5-second chirp
duration = 5;
t = 0:1/fs:5;
s1 = 2*chirp(t, 1000, 5, 40000);
% s1 = 3*sin(2*pi*30000*t);

% load up entire buffer with segments A and B
RP.WriteTagVEX('datain', 0, 'F32', s1(1:npts));
index_play = npts;

% acquisition prep
filePath = 'C:\Users\majidlab\Desktop\Navvab\TDT_Navvab\Simultaneous Continuous Play and Continuous Acquire\';
filePath = strcat(filePath, 'fnoise.F32');
fnoise = fopen(filePath,'w');

% start playing
tic
RP.SoftTrg(1);
curindex_play = RP.GetTagVal('index_play');
disp(['Current buffer index (play): ' num2str(curindex_play)]);

curindex_acquire = RP.GetTagVal('index_acquire');
disp(['Current buffer index (acquire): ' num2str(curindex_acquire)]);
index_acquire = 0;

% main looping section
sz = length(s1);

while index_play+1 < sz || index_acquire+1 < sz
	index_play
    index_acquire
	% write into play buffer *********************
    % wait until done playing A
    if index_play+1 < sz
        while(curindex_play < bufpts)
            curindex_play = RP.GetTagVal('index_play');
        end
        
        % load the next signal segment
        top = min(index_play+bufpts, sz);
        RP.WriteTagVEX('datain', 0, 'F32', s1(index_play+1:top));
        index_play = top;
        
        % checks to see if the data transfer rate is fast enough
        curindex_play = RP.GetTagVal('index_play');
        disp(['Current buffer index (play): ' num2str(curindex_play)]);
        if(curindex_play < bufpts)
            warning('Transfer rate is too slow');
        end
    end
    
	% read from acquire buffer *******************
	% wait until done writing A
    if index_acquire+1 < sz
        curindex_acquire = RP.GetTagVal('index_acquire');
        disp(['Current buffer index (acquire): ' num2str(curindex_acquire)]);
        while(curindex_acquire < bufpts)
            curindex_acquire = RP.GetTagVal('index_acquire');
        end
        % read segment A
        noise = RP.ReadTagVEX('dataout', 0, bufpts, 'F32', 'F32', 1);
        disp(['Wrote ' num2str(fwrite(fnoise,noise,'float32')) ' points to file']);
        index_acquire = index_acquire+length(noise);
        
        % checks to see if the data transfer rate is fast enough
        curindex_acquire = RP.GetTagVal('index_acquire');
        disp(['Current buffer index (acquire): ' num2str(curindex_acquire)]);
        if(curindex_acquire < bufpts)
            warning('Transfer rate is too slow');
        end
    end	
	
	% write into play buffer *********************
    % wait until start playing A 
    if index_play+1 < sz
        while(curindex_play > bufpts)
            curindex_play = RP.GetTagVal('index_play');
        end
        
        % load segment B
        top = min(index_play+bufpts, sz);
        RP.WriteTagVEX('datain', bufpts, 'F32', s1(index_play+1:top));
        index_play = top;
        
        % make sure we're still playing A
        curindex_play = RP.GetTagVal('index_play');
        disp(['Current index (play): ' num2str(curindex_play)]);
        if(curindex_play > bufpts)
            warning('Transfer rate too slow');
        end
    end 
    
	% read from acquire buffer
	% wait until start writing A 
    if index_acquire+1 < sz
        curindex_acquire = RP.GetTagVal('index_acquire');
        while(curindex_acquire > bufpts)
            curindex_acquire = RP.GetTagVal('index_acquire');
        end
        
        % read segment B
        noise = RP.ReadTagVEX('dataout', bufpts, bufpts, 'F32', 'F32', 1);
        disp(['Wrote ' num2str(fwrite(fnoise,noise,'float32')) ' points to file']);
        index_acquire = index_acquire+length(noise);
        
        % make sure we're still playing A
        curindex_acquire = RP.GetTagVal('index_acquire');
        disp(['Current index (acquire): ' num2str(curindex_acquire)]);
        if(curindex_acquire > bufpts)
            warning('Transfer rate too slow');
        end
    end
	
end

fclose(fnoise);

% stop playing
RP.SoftTrg(2);

RP.Halt;
toc

fileID = fopen('fnoise.F32');
s2 = fread(fileID,[1,length(s1)],'float32');

figure; plot(s1/max(s1));hold on; plot(s2/max(s2),'r')

[pks, locs]=findpeaks(s2);
[pks1, locs1]=findpeaks(s1(1:length(s2)));


length(locs)
length(locs1)

disp('Done Dude!')

y1 = hilbert(s1(200:end));
y2 = hilbert(s2(200:end));

insfreq1 = fs/(2*pi)*diff(unwrap(angle(y1)));
insfreq2 = fs/(2*pi)*diff(unwrap(angle(y2)));

figure; plot(insfreq1,'linewidth',2);hold on; plot(insfreq2,'r')

