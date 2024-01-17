function computeMUThreshold(physpath,animal,unit,exp,probeID,threshlength,threshlevel,name,copyToZ)
%compute the threshold for automatic MU extraction, save in separate file
%will add sampleFreq to id if not done yet; will add SpikeFiles folder if
%not done yet
%input:
% expFolder - experiment folder
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)
% threshlength - length of intervals used for threshold computations (in
% sec)
% threshlevel - scale factor for threshold (thresholds will be set to
% -threshlevel x automatically determined level)
% name - initials
% copyToZ - copy id file to Z



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
MUthresholding.threshlength=threshlength;
MUthresholding.offsetSamples=400; %nr of samples to drop at start because of filtering artefact 

baseSample=round(threshlength*id.sampleFreq);
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
MUthresholding.threshlevel=threshlevel;

save([basename '_p' num2str(probeID) '_MUthreshold.mat'],'MUthresholding');

%documentation
id.MUthreshold.date{probeID}=date;
id.MUthreshold.name{probeID}=name;
id.MUthreshold.settings{probeID}.scaleFactor = threshlevel;

save([basename '_id.mat'],'id'); %this will also add id.sampleFreq if neeed
    
if copyToZ==1
    expname=[animal '_u' unit '_' exp];
    zbase='Z:\EphysNew\processedSpikes';
    save(fullfile(zbase,animal,expname,[expname '_id.mat']),'id');
    save(fullfile(zbase,animal,expname,[expname '_p' num2str(probeID) '_MUthreshold.mat']),'MUthresholding');
end

%add SpikeFiles folder if necessary
spkbase=fullfile(physpath,animal,[animal '_u' unit '_' exp],'SpikeFiles');
if ~exist(spkbase,'dir')
    mkdir(spkbase)
end

