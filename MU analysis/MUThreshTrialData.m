function MUThreshTrialData(physpath,animal,unit,exp,probeId,eventType,eId,baseTime,stimTime)
%this function extracts spike times for MU spikes (threshold but not sorted) 
%for each trial, also computes mean number of spikes
%and rates in windows before/after event
%
%input:
%physpath: path to phys data (e.g., z:\ephysNew\processedSpikes)
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
%structure SU, one entry per cell
%fields:
%unitID: unit id from sort file
%unitClass: unit category from sort file (SU or MU)
%spktimes: spike times for every trial, in ms
%baseNspk: number of spikes in baseline period
%baseFrate: firing rate in Hz in baseline period
%stimNspk: number of spikes in stimulus period
%stimFrate: firing rate in Hz in stimulus period
%
%structure SUinfo
%fields
%eventType, eventId, baseTime & stimTime: copies of the input variables
%dom, domval, blankId,triallist: copies from trialInfo


basename=fullfile(physpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp]);

%load trial info
load([basename '_trialInfo.mat']);

%load SU data
load([basename '_p' num2str(probeId) '_MUspkMerge.mat']);

%get sampling rate
load([basename '_id.mat']); %generates id

%translate time windows into samples
baseSample=round(baseTime*id.sampleFreq);
stimSample=round(stimTime*id.sampleFreq);


%find events
if strcmp(eventType,'id')
    eventIdx=find(trialInfo.eventId==eId);
else
    eventIdx=find(trialInfo.eventCh(:,eId)==1);
end

%loop through channels and events to get data
%we're using the sorted channels here as overall organization to make the next scripts easier 
chidx=unique(MUspkMerge.detChSort);
for u=1:length(chidx)
    %base info
    MUThresh(u).detChSort=chidx(u);
    
    idx=find(MUspkMerge.detChSort==chidx(u),1,'first');
    MUThresh(u).detCh=MUspkMerge.detCh(idx);
    
    
    for i=1:length(eventIdx)
        eTime=trialInfo.eventTimes(eventIdx(i));

        sidx=find(MUspkMerge.detChSort==chidx(u) & MUspkMerge.spktimes>eTime-baseSample & MUspkMerge.spktimes<eTime+stimSample);
        spkTimesTmp=(MUspkMerge.spktimes(sidx)-eTime)/id.sampleFreq*1000;
        flagDupTmp=MUspkMerge.flagDuplicate(sidx);
        
        [MUThresh(u).spktimes{i},sortidx]=sort(spkTimesTmp);
        MUThresh(u).Nspk(i)=length(sidx);
        MUThresh(u).flagDuplicate{i}=flagDupTmp(sortidx);
        
        %SU(u).Nspk(i)=length(sidx);
        
%         bidx=find(MUspkMerge.detCh==chidx(u) & spkSort.spktimes>eTime-baseSample & spkSort.spktimes<eTime);
%         SU(u).baseNspk(i)=length(bidx);
%         SU(u).baseFrate(i)=length(bidx)/baseTime;
%        
%         stidx=find(MUspkMerge.detCh==chidx(u) & spkSort.spktimes>eTime & spkSort.spktimes<eTime+stimSample);
%         SU(u).stimNspk(i)=length(stidx);
%         SU(u).stimFrate(i)=length(stidx)/stimTime;
    end
end

%generate structure with general info, including condition info
MUThreshInfo.eventType=eventType;
MUThreshInfo.eventId=eId;
MUThreshInfo.baseTime=baseTime;
MUThreshInfo.stimTime=stimTime;
MUThreshInfo.dom=trialInfo.dom;
MUThreshInfo.domval=trialInfo.domval;
MUThreshInfo.blankId=trialInfo.blankId;
MUThreshInfo.triallist=trialInfo.triallist;

save([basename '_p' num2str(probeId) '_MUThreshTrial.mat'],'MUThresh','MUThreshInfo');

