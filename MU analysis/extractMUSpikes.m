function extractMUSpikes(expFolder,animal,unit,exp,probeID, threshlevel,name,copyToZ,parts,JobID)
%extract MU timestamps using an automated threshold
%this code implements the same spike detection algorithm as for the single
%units (minimum is spread out over time, there is a timeout after each
%crossing, but there is no combination across channels)
%artefacts are flagged by counting the number of concurrently occuring
%events; duplicates across channels within a set distance are also
%indicated
% input parameters:
% expFolder - base folder for experiments (string)
% animal - animal ID (string)
% unit - unit ID (string)
% exp - experiment ID (string)
% probe - probe ID (number)
% threshlevel - multiplication factor for automatic threshold
% name - name or initials of person running the script (for bookkeeping)
% copyToZ - copy id file to Z?
% parts - number of parts that the amplifier file will be split into
% jobID - job ID of raw spike file to process (number)
%
% output parameters:
% structure MUspk with fields (each is a vector with entries for each
% spike)
% - spkTimes: spike times for detection channels
% - detCh: id of detection channel
% - detChSort: id of detection channel, sorted according to Z and shank
% - NDuplicate: for each spike, number of channels with a simultaneously
% detect spike

%% global settings
settings.refrTime=1;
settings.refrCross=0.5; %timeout after threshold crossing in ms
settings.offsetSamples=400; %nr of samples to drop at start because of filtering artefact 
settings.filterLow=5000; %low pass filter setting
settings.filerHigh=250; %high pass filter setting
settings.spkTol=5; %window over which threshold crossings are considered duplicates/artefacts


%% basic info
%construct filenames
expname=[animal '_u' unit '_' exp];

%thresholds
load(fullfile(expFolder,animal,expname,[expname '_p' num2str(probeID) '_MUthresh.mat'])); %generates MUthresholding
chthresh=threshlevel*MUthresholding.thresholds;

%id
load(fullfile(expFolder,animal,expname,[expname '_id.mat'])); %generates id
nChannels=sum([id.probes.nChannels]);

%figure out length of recording and length 
filename=fullfile(expFolder,animal,expname,[expname '_amplifier.dat']);
fileinfo = dir(filename);
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file
samplesPerJob = ceil(samples/parts); % Number of samples to allocate to each of the  jobs

%make filter - 500Hz to 3kHz
[butter_b,butter_a] = butter(3,[settings.filerHigh settings.filterLow]/(id.sampleFreq/2),'bandpass');




%% apply threshold
fid = fopen(filename,'r');

firstSample = samplesPerJob*JobID - settings.offsetSamples; 
if firstSample<0
    firstSample=0;
end

fseek(fid,2*nChannels*firstSample,'bof'); % Offset from beginning of file

if JobID == parts-1 % The last job - first JobID is 0
    samplesLeft = samples - samplesPerJob*(parts-1) + settings.offsetSamples; % samplesLeft=TotalSamples-SamplesDone+Overhang
    Data = fread(fid, [nChannels samplesLeft], 'int16'); % If JobID is the last job, read all the samples left
else
    Data = fread(fid, [nChannels samplesPerJob], 'int16'); % If JobID isn't the last job, read samplesPerJob samples past the file position set by fseek
end

if length(id.probes)>1
    startidx=sum([id.probes(1:probeID-1).nChannels])+1; %0 for probe 1
    stopidx=startidx+id.probes(probeID).nChannels-1;
    Data=Data(startidx:stopidx,:);
end

Data=Data'; %dim 1 - time, dim 2 - channel
Data = filter(butter_b, butter_a, Data,[],1);

%spread out minimum
nTime=floor(settings.refrTime/1000*id.sampleFreq);
minData=movmin(Data,[nTime+1 nTime],1);

%get threshold crossings
AboveTh=Data>chthresh;

%get the transition points from below to above threshold (negative
CrossTh = AboveTh & circshift(~AboveTh,[-1 0]); 

%expand the threshold crossing out to cover the refractory period
nCross=floor(settings.refrCross/1000*id.sampleFreq);
CrossTh = movmax(CrossTh,[nCross 0],1);

%find spikes: Minimum across refrTime within refrCross after
%threshold crossing
Spikes = CrossTh & minData==Data; 

% Removes spikes detected in the first 1msec overlap at the beginning and end of each job. 
Spikes(1:floor(settings.offsetSamples/2),:)=0; 
Spikes(end-floor(settings.offsetSamples/2):end,:)=0; 

%% extract channel data
probeOrg=[id.probes(probeID).z id.probes(probeID).shaft id.probes(probeID).channels+1];
probeSort=sortrows(probeOrg,[2 1]);

MUspk = struct;
spkCount=1;
for i=1:size(Data,2)
    
    Times = find(Spikes(:,i)>0);
    Nspikes=length(Times);
    if Nspikes>0
        MUspk.spkTimes(spkCount:spkCount+Nspikes-1)=Times+firstSample;
        MUspk.detCh(spkCount:spkCount+Nspikes-1)=i;
        MUspk.detChSort(spkCount:spkCount+Nspikes-1)=find(probeSort(:,3)==i);
% 
%         wv=Data([-settings.spikeSamples:settings.spikeSamples]+Times,spikeData.channelIds);
%             
%             Ntime=2*settings.spikeSamples+1;
%             spikeData.rawWvfrms=reshape(wv,[Nspikes Ntime Nch]); %dimensions: spike x timepoints x channel
%             
%             %normalize by baseline
%             Nbase=floor(settings.spikeSamples/2);
%             spikeData.Wvfrms=spikeData.rawWvfrms-mean(spikeData.rawWvfrms(:,1:Nbase,:),2);
    end
    
    spkCount=spkCount+Nspikes;
end

%% flag duplicates
if isfield(MUspk,'spkTimes')
    
    %for each threshold crossing, find out how many other threshold
    %crossings occur within the tolerance window
    MUspk.NDuplicate=zeros(size(MUspk.spkTimes));

    tol=settings.spkTol/max(MUspk.spkTimes); %uniquetol scales by maximum

    [~,~,idxB]=uniquetol(MUspk.spkTimes,tol); %get unique events plus/minus tolerance
    countDup=accumarray(idxB,1); %count how often each duplicate shows up in spkTimesDet
    MUspk.NDuplicate=countDup(idxB); %build vector that gives N for each event
    MUspk.NDuplicate=MUspk.NDuplicate';


    %in addition to figuring out on how many channels an event occurs, also
%     %flag events that occur at the same time on nearby channels
% 
%     
%     for i=1:id.probes(probeID).nChannels
%         if sum(MUspk.detCh==i)>0
%             
%             %find spikes at detection channel, get their energy and position in
%             %spkTimes vector
%             detTimes=MUspk.spkTimes(MUspk.detCh==i);
%             timesIdxDet=timesIdx(spk.detCh==i);
%             enDet=spk.EnDet(spk.detCh==i);
%             
%             %spread out each event over neighboring samples (to give interval
%             %for detection)
%             detTimesConv=detTimes+[spkWindow(1):spkWindow(2)]';
%             detTimesConv=detTimesConv(:);
%             
%             %also need an index for grouping later
%             detIdxConv=repmat([1:length(detTimes)],spkWindow(2)-spkWindow(1)+1,1);
%             detIdxConv=detIdxConv(:);
%             
%             %get the other channels
%             chList=spikeData(i).channelIds; 
%             chList=chList(2:end); %need to remove center channel
%             
%             %find overlapping events
%             tf=ismember(spk.detCh,chList) & ismember(spk.spkTimesDet,detTimesConv);
%             
%             %only continue if there actually are any
%             if sum(tf)>0
%                 %get energy of overlapping events
%                 enCh=spk.EnDet(tf);
%                 
%                 %get their index in the spkTimesDet vector
%                 timesIdxCh=timesIdx(tf);
%                 
%                 %get index of original event to use as grouping par
%                 [~,locB]=ismember(spk.spkTimesDet(tf),detTimesConv);
%                 idxCh=detIdxConv(locB);
%                 
%                 %generate one big array for accumarray
%                 enAll=[enDet(unique(idxCh)) enCh]'; %we only want the relevant original spikes
%                 idxAll=[unique(idxCh);idxCh];
%                 timesIdxAll=[timesIdxDet(unique(idxCh)) timesIdxCh];
%                 
%                 %get maximum per group
%                 maxData=accumarray(idxAll,enAll,[],@max);
%                 
%                 %flag which of the overlapping events are not the maximum (i.e.
%                 %duplicates)
%                 flagD = enAll~=maxData(idxAll);
%                 
%                 %put back into large vector
%                 idxD=timesIdxAll(flagD);
%                 spk.flagDuplicate(idxD)=1; %set to zero outside loop; so events that get flagged multiple times are still ok
%             end
%         end %if sum
%     end %for ch

end



%% save
save(fullfile(expFolder,animal,expname,'SpikeFiles',[expname '_j' num2str(JobID) '_p' num2str(probeID) '_MUspk.mat']),'MUspk','settings','expname');

%% documentation
%for job 0, add info to id file for bookkeeping
if JobID==0
    id.extractMUspk.date{probeID}=date;
    id.extractMUspk.name{probeID}=name;
   
    jobVec=[0:parts-1];
    startSample = samplesPerJob*jobVec - settings.offsetSamples; 
    startSample(startSample<0)=0;
    stopSample=startSample+samplesPerJob;
    stopSample(end)=samples;
    edgeSample=startSample+settings.offsetSamples/2; %boundaries between samples
    edgeSample(end+1)=samples; %to finish the last bin
    
    id.extractMUspk.jobStart{probeID}=startSample;
    id.extractMUspk.jobStop{probeID}=stopSample;
    id.extractMUspk.jobEdges{probeID}=edgeSample;
    
    save(fullfile(expFolder,animal,expname,[expname '_id.mat']),'id'); 
    
    if copyToZ==1
        zbase='Z:\EphysNew\processedSpikes';
        save(fullfile(zbase,animal,expname,[expname '_id.mat']),'id'); 
    end
    
end

disp(['extractSpikes job ID ' num2str(JobID) ' done.'])

    