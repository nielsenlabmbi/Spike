function extractSpikes(expFolder,animalID,unitID,expID,probeID,name,copyToZ,MUflag,legacyFlag,parts,JobID,varargin)
% extractSpikes computes spike waveforms
% input parameters:
% expFolder - experiment folder
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)
% name - initials
% copyToZ - copy id file to Z
% MU flag - 1 = use MUthresholds rather than thresholds, save as MUspk
% legacyflag - implement old (less conservative) way to determine amount of
% overlap between files
% parts - number of segments to divide the data file into
% JobID - current segment to process; starts with 0
% varargin -  2 options: 'id'+id (avoids loading the id file for batch processing); 
%             'tSuffix'+suffix for threshold file, as in threshold_XXX (so multiple versions can
%             exist)

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
% generates extractSpk file with info

%% global settings
settings.refrTime=1; %timeout before and after large minima in ms
settings.refrCross=0.5; %timeout after threshold crossing in ms
settings.spikeSamplesT=[0.5 0.83]; %in ms, corresponds to [15 25] setting before, assuming 30k rate
settings.spikeRadius=50; %distance radius over which to extract spike waveforms
settings.offsetSamples=800; %this used to be partsOverlapSamples; overlap between files (increased to avoid filtering artefact)
settings.legacyFlag=legacyFlag; %for bookkeeping

%% deal with varargin
% possibilities for varargin: addition to threshold file, id 
% provided
loadId=1;
tSuffix='';
tname='threshold';
spkFolder='SpikeFiles';

if ~isempty(varargin)
    idx=find(strcmp(varargin,'tSuffix'));
    if ~isempty(idx)
        tSuffix=varargin{idx+1};
        tname=[tname '_' tSuffix];
        spkFolder=[spkFolder '_' tSuffix];
    end
    idx=find(strcmp(varargin,'id'));
    if ~isempty(idx)
        loadId=0;
        id=varargin{idx+1};
    end
end

%% generate basic info

%load threshold and id data
expname=[animalID '_u' unitID '_' expID];

if MUflag==0
    threshname=fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_' tname '.mat']);
    load(threshname); %generates thresholding
else
    threshname=fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_MU' tname '.mat']);
    load(threshname); %generates MUthresholding
    thresholding=MUthresholding;
end
if loadId==1
    load(fullfile(expFolder,animalID,expname,[expname '_id.mat'])); %generates id
end

%compute total channel number
nChannels=sum([id.probes.nChannels]);

%get file size for amplifier file
filename=fullfile(expFolder,animalID,expname,[expname '_amplifier.dat']);
fileinfo = dir(filename);
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file
samplesPerJob = ceil(samples/parts); % Number of samples to allocate to each of the 200 jobs

if legacyFlag==1
    %offset fixed at 2s previously
    settings.offsetSamples=floor((2/1000)*id.sampleFreq);
    %need to adjust the spikeSamples so that they are not in conflict with
    %the shorter window
    settings.spikeSamples(settings.spikeSamples>floor(settings.offsetSamples/2))=floor(settings.offsetSamples/2);
end

%convert spikeSamples to number
settings.spikeSamples=round(settings.spikeSamplesT/1000*id.sampleFreq); 


%% read data
firstSample = samplesPerJob*JobID - settings.offsetSamples; % Sets first sample to process; each job hasoverlap with previous job and next job
if firstSample<0
    firstSample=0;
end

DataFile = fopen(filename,'r'); % Load amplifier data
fseek(DataFile,2*nChannels*firstSample,'bof'); % Offset from beginning of file

if JobID == parts-1 % The last job - first JobID is 0
    % If JobID is the last job, read all the samples left
    samplesLeft = samples - samplesPerJob*(parts-1) + settings.offsetSamples; % samplesLeft=TotalSamples-SamplesDone+Overhang
    Data = fread(DataFile, [nChannels samplesLeft], 'int16'); 
elseif JobID == 0
    Data = fread(DataFile, [nChannels samplesPerJob], 'int16');
else
    % If JobID isn't the first or last job, read samplesPerJob+offset samples past the file position set by fseek
    Data = fread(DataFile, [nChannels samplesPerJob+settings.offsetSamples], 'int16'); 
end

fclose(DataFile);

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
Data = filter(thresholding.butter.b1, thresholding.butter.a1, Data,[],1);

% set bad channels to NaN (could also remove them)
Data(:,logical(thresholding.badChannels))=NaN; %this propagates the choice to the output



%% implement artificial threshold across time
% this is implemented by spreading out minimum peaks of voltage across time. we are 
% keeping the minimum value in a sliding window of the refractory
%period length. the window extends forwards and backwards. effectively, this causes
%a 'timeout' period before and after each strong minimum, in which other smaller minima are erased

%this code is consistent with Augusto's version, with the exception of the
%first and last data points [settings for Augusto's code: shift [31 30]]
nTime=floor(settings.refrTime/1000*id.sampleFreq);
minData=movmin(Data,[nTime+1 nTime],1);


%% Detect threshold crossings within .3msec
%detect threshold crossings
if size(thresholding.thresholds,1)==1
    AboveTh=Data>thresholding.thresholds;
else
    AboveTh=Data>thresholding.thresholds';
end
    
%get the transition points from below to above threshold (negative
%circshift shifts backwards, so this marks the last 1 (above threshold)
%before a 0 (below threshold)
CrossTh = AboveTh & circshift(~AboveTh,[-1 0]); 

%expand the threshold crossing out to cover the refractory period
nCross=floor(settings.refrCross/1000*id.sampleFreq);
CrossTh = movmax(CrossTh,[nCross 0],1);


%% Final spikes detection

%find spikes: Minimum across refrTime and refrSpace within refrCross after
%threshold crossing
%this sets the occurence of the minimum of a waveform to 1
Spikes = CrossTh & minData==Data; 


% Removes spikes detected in the overlap at the beginning and end of each job. 
%This is important as some of these may go beyond recording to get waveform.
Spikes(1:floor(settings.offsetSamples/2),:)=0; 
Spikes(end-floor(settings.offsetSamples/2)+1:end,:)=0; 


%% extract the actual waveforms

%output file - we're saving each job separately so that things can run in parallel
if MUflag==0
    outname=fullfile(expFolder,animalID,expname,spkFolder,[expname '_j' num2str(JobID) '_p' num2str(probeID)  '_spike.mat']);
else
    outname=fullfile(expFolder,animalID,expname,spkFolder,[expname '_j' num2str(JobID) '_p' num2str(probeID)  '_MUspike.mat']); 
end

for i=1:size(Data,2)
    
    if ~thresholding.badChannels(i)
        
        %initialize output - we're only extracting spike
        %data here, rest will happen in next file to make it easier to add
        %new properties
        spikeDataOut = struct;
        
        Times = find(Spikes(:,i)>0); % Find coordinates where spikes occurred
        spikeDataOut.spikeTimes=Times+firstSample;
        
        Nspikes=length(Times);
        %continue only if there are spikes
        if Nspikes>0
           
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
            spikeDataOut.channelIds=distIdx(distOrg<=settings.spikeRadius); %add offset back to get to correct channels
            Nch=length(spikeDataOut.channelIds);

            wv=Data([-settings.spikeSamples(1):settings.spikeSamples(2)]+Times,spikeDataOut.channelIds);
            
            Ntime=sum(settings.spikeSamples)+1;
            spikeDataOut.rawWvfrms=reshape(wv,[Nspikes Ntime Nch]); %dimensions: spike x timepoints x channel
            
            %normalize by baseline
            spikeDataOut.Wvfrms=spikeDataOut.rawWvfrms-mean(spikeDataOut.rawWvfrms(:,1:settings.spikeSamples(1),:),2);
            
            %save
            spikeData(1,i)=spikeDataOut;
        else
            spikeDataOut.spikeTimes=NaN;
            spikeDataOut.channelIds=NaN;
            spikeDataOut.rawWvfrms=NaN;
            spikeDataOut.Wvfrms=NaN;
            spikeData(1,i)=spikeDataOut;
            
        end
    else
        spikeDataOut.spikeTimes=NaN;
        spikeDataOut.channelIds=NaN;
        spikeDataOut.rawWvfrms=NaN;
        spikeDataOut.Wvfrms=NaN;
        spikeData(1,i)=spikeDataOut;
    end
end

save(outname,'spikeData','settings','expname','-v7.3','-nocompression');


%for job 0, save info file
if JobID==0
    jobVec=[0:parts-1];
    startSample = samplesPerJob*jobVec - settings.offsetSamples;
    startSample(startSample<0)=0;
    stopSample=startSample+samplesPerJob+settings.offsetSamples;
    stopSample(1)=samplesPerJob;
    stopSample(end)=samples;
    edgeSample=startSample+settings.offsetSamples/2; %boundaries between samples
    edgeSample(end+1)=samples; %to finish the last bin

    extractSpk.jobStart=startSample;
    extractSpk.jobStop=stopSample;
    extractSpk.jobEdges=edgeSample;
    extractSpk.date=date;
    extractSpk.name=name;
    extractSpk.settings=settings;
    extractSpk.probeNr=probeID;
    extractSpk.exptId=expname;
    extractSpk.MUflag=MUflag;

    %threshold information as available
    extractSpk.threshold.filename=tname;
    if isfield(thresholding,'date')
        extractSpk.threshold.date=thresholding.date;
    else
        fn=dir(threshname);
        fn2=strsplit(fn.date,' '); %to retain only the month
        extractSpk.threshold.date=fn2{1};
    end
    if isfield(thresholding,'name')
        extractSpk.threshold.name=thresholding.name;
    end
    
    if MUflag==0
        infoname=[expname '_p' num2str(probeID) '_extractSpk'];
    else
        infoname=[expname '_p' num2str(probeID) '_MUextractSpk'];
    end
    if ~isempty(tSuffix)
        infoname=[infoname '_' tSuffix '.mat'];
    else
        infoname=[infoname '.mat'];
    end
   

    save(fullfile(expFolder,animalID,expname,infoname),'extractSpk'); 
    
    if copyToZ==1
        zbase='Z:\EphysNew\processedSpikes';
        save(fullfile(zbase,animalID,expname,infoname),'extractSpk'); 
    end
    
end

disp(['extractSpikes job ID ' num2str(JobID) ' done.'])

        

