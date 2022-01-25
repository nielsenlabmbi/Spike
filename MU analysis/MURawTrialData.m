function MURawTrialData(physpath,animal,unit,exp,probeID,eventType,eId,baseTime,stimTime,channel)
%compute raw, continuous MU signal for every trial, for a single channel
%only

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
%channel: which channel to extract
%
%output:
%MURaw: 2D matrix with dimensions time x event
%MURawInfo: parameters for MURaw computation 

basename=fullfile(physpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp]);

%load id file (for recording and probe info)
load([basename '_id.mat']); %generates id
nrChTotal=sum([id.probes.nChannels]);

%load trialinfo file (for trial information)
load([basename '_trialInfo.mat']); %generates trialInfo

%make filter - 500Hz to 3kHz
[butter_b,butter_a] = butter(3,[250 5000]/(id.sampleFreq/2),'bandpass');

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
   
    %only keep the correct channel
    chidx=sum([id.probes(1:probeID-1).nChannels])+channel; 

    Data=Data(chidx,:);
    
    %filter
    Data = filter(butter_b, butter_a, Data);
    
    %collect
    MURaw(:,i)=Data;

end


%document settings
MURawInfo.eventType=eventType;
MURawInfo.eventId=eId;
MURawInfo.baseTime=baseTime;
MURawInfo.stimTime=stimTime;
MURawInfo.channel=channel;
MURawInfo.filter=[butter_a,butter_b];




%save
save([basename '_c' num2str(channel) '_p' num2str(probeID) '_MURawTrial.mat'],'MURaw','MURawInfo');

