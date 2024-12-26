function computeMUThreshold(physpath,animal,unit,exp,probeID,threshlength,threshlevel,name,copyToZ,varargin)
%compute the threshold for automatic MU extraction, save in separate file
%will add sampleFreq to id if not done yet; will add SpikeFiles folder if
%not done yet
%input:
% physpath - experiment folder
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)
% threshlength - length of intervals used for threshold computations (in
% sec); if left empty set to a quarter of the file
% threshlevel - scale factor for threshold (thresholds will be set to
% -threshlevel x automatically determined level)
% name - initials
% copyToZ - copy id file to Z
% varargin: allows to add an addition to the file name (to distinguish
% different versions)



%construct file name
basename=fullfile(physpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp]);

%need id file for number of channels and sampling rate
load([basename '_id.mat']); %generates id
nChannels=sum([id.probes.nChannels]);

%id may not have sampleFreq yet - add here if yes
if ~isfield(id,'sampleFreq')
    headerName=[basename '_info.rhd'];
    header=read_Intan_Header(headerName);
    id.sampleFreq = header.sample_rate;
end

%figure out length of recording and length  
fileinfo = dir([basename '_amplifier.dat']);
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file

%construct filter
MUthresholding.butter.lowPass=5000; %low pass filter setting
MUthresholding.butter.highPass=250; %high pass filter setting

[butter_b,butter_a] = butter(3,[MUthresholding.butter.highPass MUthresholding.butter.lowPass]/(id.sampleFreq/2),'bandpass');
MUthresholding.butter.b1=butter_b;
MUthresholding.butter.a1=butter_a;

%for compatibility
MUthresholding.badChannels=zeros(1,nChannels);


%samples and time points for threshold - threshold is fixed throughout recording period, so
%use 3 segments to figure out level, one at start one in middle, one at end
if isempty(threshlength)
    MUthresholding.threshlength=(samples/id.sampleFreq)*.25;
else
    MUthresholding.threshlength=threshlength;
end
MUthresholding.offsetSamples=400; %nr of samples to drop at start because of filtering artefact 

baseSample=round(MUthresholding.threshlength*id.sampleFreq);
startSample(1)=0;
startSample(2)=samples/2;
startSample(3)=samples-baseSample;

%compute threshold
fid = fopen([basename '_amplifier.dat'],'r');
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
    Data=Data(MUthresholding.offsetSamples:end,:);
    
    %compute threshold
    chthresh(s,:) = squeeze(round(1.4826 * median(abs(Data - median(Data,1)),1)));
end
MUthresholding.thresholds=-threshlevel*mean(chthresh,1);

%add documentation to MUthresholding
MUthresholding.threshlevel=threshlevel;
MUthresholding.threshlength=threshlength;
MUthresholding.date=date;
MUthresholding.name=name;
MUthresholding.probeNr=probeID;
MUthresholding.exptId=[animal '_u' unit '_' exp];
if ~isempty(varargin)
    MUthresholding.suffix=varargin{1};
end


if ~isempty(varargin)
    save([basename '_p' num2str(probeID) '_MUthreshold_' varargin{1} '.mat'],'MUthresholding');
else
    save([basename '_p' num2str(probeID) '_MUthreshold.mat'],'MUthresholding');
end

    
if copyToZ==1
    expname=[animal '_u' unit '_' exp];
    zbase='Z:\EphysNew\processedSpikes';

    if ~isempty(varargin)
        save(fullfile(zbase,animal,expname,[expname '_p' num2str(probeID) '_MUthreshold_' varargin{1} '.mat']),'MUthresholding');
    else
        save(fullfile(zbase,animal,expname,[expname '_p' num2str(probeID) '_MUthreshold.mat']),'MUthresholding');
    end
    
end

%add SpikeFiles folder if necessary
if isempty(varargin)
    spkbase=fullfile(physpath,animal,[animal '_u' unit '_' exp],'SpikeFiles');
else
    spkbase=fullfile(physpath,animal,[animal '_u' unit '_' exp],['SpikeFiles_' varargin{1}]);
end

if ~exist(spkbase,'dir')
    mkdir(spkbase)
end

