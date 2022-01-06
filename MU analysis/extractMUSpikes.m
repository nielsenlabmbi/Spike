function extractMUSpikes(physpath,animal,unit,exp,probeID, threshlevel,threshlength,name,copyToZ,parts,JobID)
%extract MU timestamps using an automated threshold


%% global settings
settings.refrTime=1;
settings.refrCross=0.5; %timeout after threshold crossing in ms
settings.offsetSamples=200; %nr of samples to drop at start because of filtering artefact 
settings.filterLow=5000; %low pass filter setting
settings.filerHigh=250; %high pass filter setting


%% basic inforead data
%construct filenames
basename=fullfile(physpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp]);


%need id file for number of channels and sampling rate
load([basename '_id.mat']); %generates id
nChannels=sum([id.probes.nChannels]);

%figure out length of recording and length  
fileinfo = dir([basename '_amplifier.dat']);
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file
samplesPerJob = ceil(samples/parts); % Number of samples to allocate to each of the  jobs
partsOverlapSamples = floor((2/1000)*id.sampleFreq); % get 2msec overlap between samples

%make filter - 500Hz to 3kHz
[butter_b,butter_a] = butter(3,[settings.filerHigh settings.filterLow]/(id.sampleFreq/2),'bandpass');


%% get threshold
%samples and time points for threshold - threshold is fixed throughout recording period, so
%use 3 segments to figure out level, one at start one in middle, one at end
baseSample=round(threshlength*id.sampleFreq);
startSample(1)=0;
startSample(2)=samples/2;
startSample(3)=samples-baseSample;


fid = fopen([basename '_amplifier.dat'],'r');

%compute threshold
frewind(fid);
for s=1:3
    fseek(fid,2*startSample(s)*nChannels,'bof');
    Data = fread(fid, [nChannels baseSample], 'int16');
    
    if length(id.probes)>1
        startidx=sum([id.probes(1:probeID-1).nChannels])+1; %0 for probe 1
        stopidx=startidx+id.probes(probeID).nChannels-1;
        Data=Data(startidx:stopidx,:);
    end
    
    Data=Data'; %dim 1 - time, dim 2 - channel
    Data = filter(butter_b, butter_a, Data,[],1);

    %there is a filter artefact at the start that needs to be avoided
    Data=Data(settings.offsetSamples:end,:);
    
    %compute threshold
    chthresh(s,:) = squeeze(round(1.4826 * median(abs(Data - median(Data,1)),1)));
end
chthresh=-threshlevel*mean(chthresh,1);

%% apply threshold
%now read actual data and threshold
firstSample = samplesPerJob*JobID - partsOverlapSamples; 
if firstSample<0
    firstSample=0;
end

fseek(fid,2*nChannels*firstSample,'bof'); % Offset from beginning of file

if JobID == parts-1 % The last job - first JobID is 0
    samplesLeft = samples - samplesPerJob*(parts-1) + partsOverlapSamples; % samplesLeft=TotalSamples-SamplesDone+Overhang
    Data = fread(DataFile, [nChannels samplesLeft], 'int16'); % If JobID is the last job, read all the samples left
else
    Data = fread(DataFile, [nChannels samplesPerJob], 'int16'); % If JobID isn't the last job, read samplesPerJob samples past the file position set by fseek
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
Spikes(1:floor(partsOverlapSamples/2),:)=0; 
Spikes(end-floor(partsOverlapSamples/2):end,:)=0; 

%% extract channel data
MUSpikesData = struct;
for i=1:size(Data,2)
    
    Times = find(Spikes(:,i)>0);
    if ~isempty(Times)
        MUSpikesData(i).spikeTimes=Times+firstSample;
    else
        MUSpikesData(i).spikeTimes=NaN;
    end
    MUSpikesData(i).threshold=chthresh(i);
end

MUSpikesinfo.xpos=id.probes(probeID).x;
MUSpikesinfo.zpos=id.probes(probeID).z;
MUSpikesinfo.shaft=id.probes(probeID).shaft;
MUSpikesinfo.zshaft=id.probes(probeID).z+100*(id.probes(probeID).shaft-1);%buffer of 100um between shafts

expname=[animalID '_u' unitID '_' expID];
save([basename '_p' num2str(probeID) '_MUSpikes.mat'],'MUSpikesData','MUSpikesinfo','settings','expname');

%% documentation
%for job 0, add info to id file for bookkeeping
if JobID==0
    id.extractMUSpikes.date=date;
    id.extractMUSpikes.name=name;
   
    jobVec=[0:parts-1];
    startSample = samplesPerJob*jobVec - partsOverlapSamples; 
    startSample(startSample<0)=0;
    stopSample=startSample+samplesPerJob;
    stopSample(end)=samples;
    edgeSample=startSample+partsOverlapSamples/2; %boundaries between samples
    edgeSample(end+1)=samples; %to finish the last bin
    
    id.extractMUSpikes.jobStart=startSample;
    id.extractMUSpikes.jobStop=stopSample;
    id.extractMUSpikes.jobEdges=edgeSample;
    
    save([basename '_id.mat'],'id'); 
    
    if copyToZ==1
        zbase='Z:\EphysNew\processedSpikes';
        save(fullfile(zbase,animalID,expname,[expname '_id.mat']),'id'); 
    end
    
end

disp(['extractSpikes job ID ' num2str(JobID) ' done.'])

    