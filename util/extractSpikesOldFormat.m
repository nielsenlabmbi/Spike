function extractSpikesOldFormat(expFolder,animalID,unitID,expID,probeID,name,copyToZ,parts,JobID)
% extractSpikesOldFormat extracts spike waveforms based on timestamps
% generated with the old spike sorter, using the data stored in spkSort
% input parameters:
% expFolder - experiment folder
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)
% parts - number of segments to divide the data file into
% JobID - current segment to process; starts with 0
% name - initials
% copyToZ - copy id file to Z

% output (saved in matlab file)
% spikeData - structure with fields (one per detection channel)
%    spikeTimes - time of minimum at detection channel (in samples, not seconds)
%    channelIds - channels for which waveforms are extracted (within set
%    distance from detection channel); IDs are numbered as in acquisition
%    file, not according to location
%    rawWvfrms - raw waveforms at the selected channel in a window around
%    the timestamp in spikeTimes; spikes x time x channel
%    Wvfrms - same waveforms as in rawWvfrms, but baseline correct using
%    start of each waveform

% requires SpikeFiles directory to be present in the experiment folder
% also updates the id file

%% global settings
settings.spikeSamples=15; %number of sample points per spike before and after the spike time
settings.spikeRadius=100; %distance radius over which to extract spike waveforms

%% generate basic info
%load spkSort and id data
expname=[animalID '_u' unitID '_' expID];

%test whether spikefiles folder exists (usually generated using threshold gui)
if exist(fullfile(expfolder,animalID,expname,'SpikeFiles'),'dir')~=7
    disp('SpikeFiles folder does not exist! Cannot proceed.')
    return;
end

load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_spkSort'])); %generates spkSort
load(fullfile(expFolder,animalID,expname,[expname '_id.mat'])); %generates id

%filter settings
hp=250;
lp=5000;
[b1,a1]=butter(3,[hp/id.sampleFreq,lp/id.sampleFreq]*2,'bandpass');

%compute total channel number
nChannels=sum([id.probes.nChannels]);

%get file size
filename=fullfile(expFolder,animalID,expname,[expname '_amplifier.dat']);
fileinfo = dir(filename);
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file
samplesPerJob = ceil(samples/parts); % Number of samples to allocate to each of the 200 jobs

%padding
partsOverlapSamples = floor((2/1000)*id.sampleFreq); % get 2msec overlap between samples


%% read data - we could just read snippets, but easier to read all at once
%file boundaries (in terms of samples) for time stamps
startFile=samplesPerJob*JobID; %first job is job 0
stopFile=samplesPerJob*(JobID+1);

%file boundaries for reading (with overlap)
firstSample = samplesPerJob*JobID - partsOverlapSamples; % Sets first sample to process; 
if firstSample<0
    firstSample=0;
end

DataFile = fopen(filename,'r'); % Load amplifier data
fseek(DataFile,2*nChannels*firstSample,'bof'); % Offset from beginning of file

if JobID == parts-1 % The last job - first JobID is 0
    samplesLeft = samples - samplesPerJob*(parts-1) + partsOverlapSamples; % samplesLeft=TotalSamples-SamplesDone
    Data = fread(DataFile, [nChannels samplesLeft], 'int16'); % If JobID is the last job, read all the samples left
else
    Data = fread(DataFile, [nChannels samplesPerJob+2*partsOverlapSamples], 'int16'); % If JobID isn't the last job, read samplesPerJob samples past the file position set by fseek
end

%extract only data for the relevant channels (only necessary if there are 2
%probes)
if length(id.probes)>1
    startidx=sum([id.probes(1:probeID-1).nChannels])+1; %0 for probe 1
    stopidx=startidx+id.probes(probeID).nChannels-1;
    Data=Data(startidx:stopidx,:);
end
    
%transpose data - matlab is faster with the longer dimension as rows rather
%than columns
Data=Data';

% Filter 
Data = filter(b1, a1, Data,[],1);


%% extract the actual waveforms

%output file - we're saving each job separately so that things can run in parallel
outname=fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(JobID) '_p' num2str(probeID)  '_spike.mat']); 
matOut=matfile(outname,'Writable',true);

%add settings and original file name for record keeping
matOut.settings=settings;
matOut.expname=expname;

for i=1:size(Data,2)

    %initialize output - we're collecting things in a structure here, to
    %add it to the matfile object later; we're only extracting spike
    %data here, rest will happen in next file to make it easier to add
    %new properties
    spikeData = struct;

    spkidx=find(spkSort.spktimes>=startFile & spkSort.spktimes<=stopFile & spkSort.detCh==i);

    Nspikes=length(spkidx);
    %continue only if there are spikes
    if Nspikes>0

        %get spike times
        spikeData.spikeTimes=spkSort.spktimes(spkidx);
        Times=spikeData.spikeTimes'-firstSample+1;

        %extract waveforms - we are ignoring the shank here, since shanks might
        %be close enough to pick up the same waveforms
        %the number of channels in this radius will be variable across
        %channels, but should be very similar for neighboring channels and will
        %be constant for each channel
        %we're reorganizing things according to distance to the
        %detection channel, which also makes the detection channel
        %the first entry in the waveform matrices
        distCh=sqrt((id.probes(probeID).x-id.probes(probeID).x(i)).^2+(id.probes(probeID).z-id.probes(probeID).z(i)).^2);
        [distOrg,distIdx]=sort(distCh);
        spikeData.channelIds=distIdx(distOrg<=settings.spikeRadius); %add offset back to get to correct channels
        Nch=length(spikeData.channelIds);

        if min(-settings.spikeSamples+Times)<1
            disp('Index too small')
        end
        if max(settings.spikeSamples+Times)>size(Data,1)
            disp('Index too large')
            %disp(Times)
        end
        
        wv=Data([-settings.spikeSamples:settings.spikeSamples]+Times,spikeData.channelIds);

        Ntime=2*settings.spikeSamples+1;
        spikeData.rawWvfrms=reshape(wv,[Nspikes Ntime Nch]); %dimensions: spike x timepoints x channel

        %normalize by baseline
        Nbase=floor(settings.spikeSamples/2);
        spikeData.Wvfrms=spikeData.rawWvfrms-mean(spikeData.rawWvfrms(:,1:Nbase,:),2);

        %save
        matOut.spikeData(1,i)=spikeData;
    else
        spikeData.spikeTimes=NaN;
        spikeData.channelIds=NaN;
        spikeData.rawWvfrms=NaN;
        spikeData.Wvfrms=NaN;
        matOut.spikeData(1,i)=spikeData;

    end
    
end


%for job 0, add info to id file for bookkeeping
if JobID==0
    id.extractSpikes.date=date;
    id.extractSpikes.name=name;
   
    jobVec=[0:parts-1];
    startSample = samplesPerJob*jobVec; 
    startSample(startSample<0)=0;
    stopSample=startSample+samplesPerJob;
    stopSample(end)=samples;
    edgeSample=startSample; %boundaries between samples
    edgeSample(end+1)=samples; %to finish the last bin
    
    id.extractSpikes.jobStart=startSample;
    id.extractSpikes.jobStop=stopSample;
    id.extractSpikes.jobEdges=edgeSample;
    
    save(fullfile(expFolder,animalID,expname,[expname '_id.mat']),'id'); 
    
    if copyToZ==1
        zbase='Z:\EphysNew\processedSpikes';
        save(fullfile(zbase,animalID,expname,[expname '_id.mat']),'id'); 
    end
    
end

disp(['extractSpikes job ID ' num2str(JobID) ' done.'])

        

