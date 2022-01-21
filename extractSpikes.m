function extractSpikes(expFolder,animalID,unitID,expID,probeID,name,copyToZ,MUflag,legacyFlag,parts,JobID)
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
settings.refrTime=1; %timeout before and after large minima in ms
settings.refrCross=0.5; %timeout after threshold crossing in ms
settings.spikeSamples=[15 25]; %number of sample points per spike before and after the minimum, used to be [15 15]
settings.spikeRadius=100; %distance radius over which to extract spike waveforms
settings.offsetSamples=400; %this used to be partsOverlapSamples; overlap between files (increased to avoid filtering artefact)
settings.legacyFlag=legacyFlag; %for bookkeeping

%% generate basic info
%load threshold and id data
expname=[animalID '_u' unitID '_' expID];
if MUflag==0
    load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_threshold.mat'])); %generates thresholding
else
    load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_MUthreshold.mat'])); %generates MUthresholding
    thresholding=MUthresholding;
end
load(fullfile(expFolder,animalID,expname,[expname '_id.mat'])); %generates id

%compute total channel number
nChannels=sum([id.probes.nChannels]);

%get file size
filename=fullfile(expFolder,animalID,expname,[expname '_amplifier.dat']);
fileinfo = dir(filename);
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file
samplesPerJob = ceil(samples/parts); % Number of samples to allocate to each of the 200 jobs

if legacyFlag==1
    settings.offsetSamples=floor((2/1000)*id.sampleFreq);
    %need to adjust the spikeSamples so that they are not in conflict with
    %the shorter window
    settings.spikeSamples(settings.spikeSamples>floor(settings.offsetSamples/2))=floor(settings.offsetSamples/2);
end


%% read data
firstSample = samplesPerJob*JobID - settings.offsetSamples; % Sets first sample to process; each job has 1st 2msec overlap with previous job and last 2msec overlap with next job
if firstSample<0
    firstSample=0;
end

DataFile = fopen(filename,'r'); % Load amplifier data
fseek(DataFile,2*nChannels*firstSample,'bof'); % Offset from beginning of file

if JobID == parts-1 % The last job - first JobID is 0
    samplesLeft = samples - samplesPerJob*(parts-1) + settings.offsetSamples; % samplesLeft=TotalSamples-SamplesDone+Overhang
    Data = fread(DataFile, [nChannels samplesLeft], 'int16'); % If JobID is the last job, read all the samples left
else
    Data = fread(DataFile, [nChannels samplesPerJob], 'int16'); % If JobID isn't the last job, read samplesPerJob samples past the file position set by fseek
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
% this is implemented by spreading out minimum peaks of voltage across time and
% space (channels, see next part). we are keeping the minimum value in a sliding window of the refractory
%period length. the window extends forwards and backwards. effectively, this causes
%a 'timeout' period before and after each strong minimum, in which other smaller minima are erased

%this code is consistent with Augusto's version, with the exception of the
%first and last data points [settings for Augusto's code: shift [31 30]]
nTime=floor(settings.refrTime/1000*id.sampleFreq);
minData=movmin(Data,[nTime+1 nTime],1);


%% implement artificial threshold across space
%similar to the logic for the artifical threshold across time, this spreads
%out a minimum across channels
%the faster way would be to use movmin after sorting channels according to
%position; not used here because the movmin approach makes assumptions
%about regular spacing between probe sites that often are violated

%base code fragment for movmin
%for i=1:id.probes.nShanks
%    minTmp=minData(:,id.probes.config(id.probes.shaft==i));
%    minData(:,1+(i-1)*id.probes.nChannels/id.probes.nShanks:i*id.probes.nChannels/id.probes.nShanks)=movmin(minTmp,9,2);
%end

% if settings.useRefrSpace
%     %this should only be executed for the probes; we want to keep the signals
%     %from tetrode wires independent
%     if length(id.probes)>1 || ~strcmp(id.probes.type,'single') && ~strcmp(id.probes.type,'tetrode')
%         
%         minDataTmp=minData;
%         
%         for p=1:length(id.probes)
%             for i=1:id.probes(p).nChannels
%                 
%                 offsetCh=sum([id.probes(1:p-1).nChannels]); %0 for p=1
%                 
%                 if ~thresholding.badChannels(i+offsetCh)
%                     distCh=sqrt((id.probes(p).x-id.probes(p).x(i)).^2+(id.probes(p).z-id.probes(p).z(i)).^2);
%                     mask=(distCh<=settings.refrSpace);
%                     
%                     maskFull=false(nChannels,1);
%                     maskFull(offsetCh+1:offsetCh+length(mask),1)=mask;
%                     
%                     minData(:,i+offsetCh)=min(minDataTmp(:,maskFull),[],2);
%                 end
%             end
%         end
%         
%         clear minDataTmp;
%     end
% end


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


%% Make sure there is no repeated value of max during artificial refractory period (necessary for raw data with low bit depth)
%double - unlikely this will actually ever occur, so removed for now
% RepeatedMax = zeros(size(Data));
% % Is the first max in 15 samples across neighbor channels
% for ch = -8:8 % channel window
%     for sm = -30:0 % sample window
%         if ~(sm==0 && ch==0)
%         RepeatedMax = RepeatedMax | Data == circshift(Data,[sm ch]);
%         end
%     end
% end

%% Final spikes detection

%find spikes: Minimum across refrTime and refrSpace within refrCross after
%threshold crossing
%this sets the occurence of the minimum of a waveform to 1
Spikes = CrossTh & minData==Data; 

% Removes spikes detected in the first 1msec overlap at the beginning and end of each job. 
%This is important as some of these may go beyond recording to get waveform.
Spikes(1:floor(settings.offsetSamples/2),:)=0; 
Spikes(end-floor(settings.offsetSamples/2):end,:)=0; 


%% extract the actual waveforms

%output file - we're saving each job separately so that things can run in parallel
if MUflag==0
    outname=fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(JobID) '_p' num2str(probeID)  '_spike.mat']);
else
    outname=fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(JobID) '_p' num2str(probeID)  '_MUspike.mat']); 
end
matOut=matfile(outname,'Writable',true);

%add settings and original file name for record keeping
matOut.settings=settings;
matOut.expname=expname;


for i=1:size(Data,2)
    
    if ~thresholding.badChannels(i)
        
        %initialize output - we're collecting things in a structure here, to
        %add it to the matfile object later; we're only extracting spike
        %data here, rest will happen in next file to make it easier to add
        %new properties
        spikeData = struct;
        
        Times = find(Spikes(:,i)>0); % Find coordinates where spikes occurred
        spikeData.spikeTimes=Times+firstSample;
        
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
            spikeData.channelIds=distIdx(distOrg<=settings.spikeRadius); %add offset back to get to correct channels
            Nch=length(spikeData.channelIds);

            wv=Data([-settings.spikeSamples(1):settings.spikeSamples(2)]+Times,spikeData.channelIds);
            
            Ntime=sum(settings.spikeSamples)+1;
            spikeData.rawWvfrms=reshape(wv,[Nspikes Ntime Nch]); %dimensions: spike x timepoints x channel
            
            %normalize by baseline
            spikeData.Wvfrms=spikeData.rawWvfrms-mean(spikeData.rawWvfrms(:,1:settings.spikeSamples(1),:),2);
            
            %save
            matOut.spikeData(1,i)=spikeData;
        else
            spikeData.spikeTimes=NaN;
            spikeData.channelIds=NaN;
            spikeData.rawWvfrms=NaN;
            spikeData.Wvfrms=NaN;
            matOut.spikeData(1,i)=spikeData;
            
        end
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
    jobVec=[0:parts-1];
    startSample = samplesPerJob*jobVec - settings.offsetSamples;
    startSample(startSample<0)=0;
    stopSample=startSample+samplesPerJob;
    stopSample(end)=samples;
    edgeSample=startSample+settings.offsetSamples/2; %boundaries between samples
    edgeSample(end+1)=samples; %to finish the last bin

  


    if MUflag==0
        %need to clean up previous versions that didn't index according to
        %probe
        if isfield(id,'extractSpikes')
            if ~iscell(id.extractSpikes.date) %name, date, extractSpikes go together, so only do this once
                %if there is only one probe, or only one probe has ever been thresholded
                %simply delete since it's now obsolete
                if length(id.probes)==1 || sum(id.threshold.processedProbe)==1
                    id=rmfield(id,'extractSpikes');
                else %try to keep information, using date info
                    t1=datetime([id.threshold.date]);
                    t2=datetime(id.extractSpikes.date);
                    [~,oldProbe]=min(t2-t1); %only returns 1 value

                    tmpId=id.extractSpikes;
                    id=rmfield(id,'extractSpikes');

                    id.extractSpikes.date{oldProbe}=tmpId.date;
                    id.extractSpikes.name{oldProbe}=tmpId.name;
                    id.extractSpikes.jobStart{oldProbe}=tmpId.jobStart;
                    id.extractSpikes.jobStop{oldProbe}=tmpId.jobStop;
                    id.extractSpikes.jobEdges{oldProbe}=tmpId.jobEdges;
                end
            end
        end

        id.extractSpikes.date{probeID}=date;
        id.extractSpikes.name{probeID}=name;

        id.extractSpikes.jobStart{probeID}=startSample;
        id.extractSpikes.jobStop{probeID}=stopSample;
        id.extractSpikes.jobEdges{probeID}=edgeSample;

        id.extractSpikes.legacyFlag(probeID)=legacyFlag;
    else
        id.MUextractSpikes.date{probeID}=date;
        id.MUextractSpikes.name{probeID}=name;

        id.MUextractSpikes.jobStart{probeID}=startSample;
        id.MUextractSpikes.jobStop{probeID}=stopSample;
        id.MUextractSpikes.jobEdges{probeID}=edgeSample;
    end
    
    save(fullfile(expFolder,animalID,expname,[expname '_id.mat']),'id'); 
    
    if copyToZ==1
        zbase='Z:\EphysNew\processedSpikes';
        save(fullfile(zbase,animalID,expname,[expname '_id.mat']),'id'); 
    end
    
end

disp(['extractSpikes job ID ' num2str(JobID) ' done.'])

        

