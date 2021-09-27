clear;
clc

NTones = 4;
NperTone = 10;
azList = 0:30:180;
NAz = length(azList);
azNum = 1:NAz;
Ntotal = NTones*NperTone*NAz;
maxTravel = (NAz-1)/2; % = 3

ISI = 4.5; % in sec
jitter = 15; % in percentage
jitter = jitter*ISI/100; % in sec
NCam = 5;
nFrames = 7000;
FsCam = 30;
TCam = 1/FsCam;
DurPerCam = nFrames/FsCam;
NperCam = Ntotal/NCam;

for seed = 1:10
%     clc
    rng(1000*(seed));
    AzOrder = repmat([1:NAz]',NperTone,NTones);
    ToneOrder = repmat(1:NTones,NperTone*NAz,1);
    p1 = randperm(Ntotal);
    ToneOrder = ToneOrder(p1);
    AzOrder = AzOrder(p1);
    permute = true;
    while permute
        travel = abs(diff(AzOrder))>maxTravel;
        dmin = 2;
        if mod(seed,20)==0
            dmin = 3;
        end
        if sum(travel)<=dmin
            permute = false;
            break;
        end
        idx1 = [1,find(~travel)+1];
        idx2 = find(travel)+1;
        p1 = randperm(length(idx2));
        idx2 = idx2(p1);
        ToneOrder = [ToneOrder(idx1),ToneOrder(idx2)];
        AzOrder = [AzOrder(idx1),AzOrder(idx2)];
    end
    % timing
    travel = abs(diff([1,AzOrder]))+rand(1,Ntotal);
    travel = travel/max(travel);
    travel = travel*0.4/mean(travel);
    minISI = ISI - jitter;
    TimeOrder = minISI + travel*2*jitter;
    TimeOrder = TCam*round(TimeOrder/TCam);
    
    if ~isnan(DurPerCam)
        for ncam = 1:NCam
            idx = 1+(ncam-1)*NperCam : ncam*NperCam;
            timeOrderTemp = TimeOrder(idx);
            s = sum(timeOrderTemp);
            if s > DurPerCam-ISI % sum of silence time be at least ISI sec less than total recording time
                d = (s+ISI-DurPerCam)/length(timeOrderTemp);
                timeOrderTemp = timeOrderTemp-d;
                TimeOrder(idx) = timeOrderTemp;
            end
        end
    end
    [min(TimeOrder),mean(TimeOrder),std(TimeOrder),max(TimeOrder),sum(TimeOrder)-1150]
    figure(1);clf;plot(TimeOrder,'g');hold on;plot(abs(diff([1,AzOrder]))+2.5,'r')
    figure(2);clf;hist([TimeOrder;abs(diff([1,AzOrder]))].',[0:0.3:6])
    pause(1)
    nnn = 0;
    travel  = abs(diff([1,AzOrder]));
    mtime = travel*4.25/(NAz-1); % motor travel time. 4.25sec for full one way travel with speed 250
    br = TimeOrder-mtime; % sum of baseline and response time
    brNorm = (br/mean(br));
    brNorm = brNorm/(3*std(brNorm))/5; % make brNorm to have the change range of 1/5 (now we have 3*std(br) = 1/5)
    btime = 0.5+rand(size(brNorm)).*brNorm/1; % baseline time
    rtime = TimeOrder-btime-mtime; % response time
    figure(3);clf;hold on;plot(btime,'r');plot(rtime)
    [mean(btime) mean(rtime) mean(mtime)]
    % saving
    FolderPath = [pwd,'\\Tone-Azimuth-Time Orders'];
    fileNameMat = sprintf('Sound Localization_Moving Speaker_Seed_%d.mat',seed);
    filePathMat = sprintf('%s\\%s',FolderPath,fileNameMat);
    save(filePathMat,'NTones','NperTone','NAz','maxTravel','ToneOrder','AzOrder',...
        'ISI','jitter','NCam','DurPerCam','TimeOrder','mtime','btime','rtime')
    fileNameTxt = sprintf('Sound Localization_Moving Speaker_Seed_%d.txt',seed);
    filePathTxt = sprintf('%s\\%s',FolderPath,fileNameTxt);
    fid = fopen(filePathTxt,'w');
    
    fprintf(fid,'%9s %9s %12s %9s %9s %9s\n',...
                   'ToneOrder','AzOrder','TimeOrder', 'mtime','btime','rtime');
    fprintf(fid,'%9d %9d %12.4f %9.4f %9.4f %9.4f\n',...
                   [ToneOrder;AzOrder;TimeOrder;mtime;btime;rtime]);
    fclose(fid);
end