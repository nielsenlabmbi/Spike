function MUEnvTrialData(physpath,animal,unit,exp,probeID,eventType,eId,baseTime,stimTime)
%compute MUEnv (Nikos method) for every trial

%input:
%physpath: path to amplifier file
%animal: animal id
%unit: unit id (string)
%exp: exp id (string)
%probeId: probe id (number)
%eventType: id - decimal number for event; ch - 'on' for 1 selected channel
%eId: either decimal number for trigger event, or channel number
%baseTime: time before event to include (in s)
%stimTime: time after event to include (in s)
%
%output:
%MUEnv: matrix with dimensions time x channel x event
%zMUEnv: z transformed MUEnv, using mean and standard deviation from
%baseline period
%MUEnvInfo: parameters for MU computation

basename=fullfile(physpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp]);

%load id file (for recording and probe info)
load([basename '_id.mat']); %generates id
nrChTotal=sum([id.probes.nChannels]);

%load trialinfo file (for trial information)
load([basename '_trialInfo.mat']); %generates trialInfo

%make filter - 500Hz to 3kHz
[butter_b,butter_a] = butter(3,[500 3000]/(id.sampleFreq/2),'bandpass');

%translate time windows into samples
baseSample=round(baseTime*id.sampleFreq);
stimSample=round(stimTime*id.sampleFreq);
nSamples=baseSample+stimSample;

%find events
if strcmp(eventType,'id')
    eventIdx=find(trialInfo.eventId==eId);
else
    eventIdx=find(trialInfo.eventCh(:,eId)==1);
end


%open amplifier file
dataFileId = fopen([basename '_amplifier.dat'],'r');


for i=1:length(eventIdx)
    eTime=trialInfo.eventTimes(eventIdx(i)); %in samples
    startSample=eTime-baseSample;
    
    %read all data
    frewind(dataFileId);
    fseek(dataFileId,2*startSample*nrChTotal,'bof');
    Data = fread(dataFileId, [nrChTotal nSamples], 'int16');
   
    %only keep the correct probe and reshape
    if length(id.probes)>1
        startidx=sum([id.probes(1:probeID-1).nChannels])+1; %0 for probe 1
        stopidx=startidx+id.probes(probeID).nChannels-1;
        Data=Data(startidx:stopidx,:);
    end
    
    %transpose for speed - dimension 1 is samples, dimension 2 channels
    Data=Data';
    
    %filter
    Data = filter(butter_b, butter_a, Data,[],1);
    
    %rectify
    Data=abs(Data);

    %decimate to 1kHz - need to do this for every channel separately
    for c=1:id.probes(probeID).nChannels
        decData(:,c)=decimate(Data(:,c),id.sampleFreq/1000);
    end
    
    
    %also compute Z transform
    baseSampleDec=baseSample/(id.sampleFreq/1000);
    meanBase=mean(decData(1:baseSampleDec,:),1);
    stdBase=std(decData(1:baseSampleDec,:),0,1);
    zDecData=(decData-meanBase)./stdBase;
    
    %collect
    MUEnv(:,:,i)=decData;
    zMUEnv(:,:,i)=zDecData;
end


%document settings
MUEnvInfo.eventType=eventType;
MUEnvInfo.eventId=eId;
MUEnvInfo.baseTime=baseTime;
MUEnvInfo.stimTime=stimTime;
MUEnvInfo.filter=[butter_a,butter_b];


%save
save([basename '_p' num2str(probeID) '_MUEnvTrial.mat'],'MUEnv','zMUEnv','MUEnvInfo');

