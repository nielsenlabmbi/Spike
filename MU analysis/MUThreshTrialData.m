function MUThreshTrialData(physpath,animal,unit,exp,probeId,eventType,eId,baseTime,stimTime,varargin)
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
%stimTime: time after event to include (in s)
%varargin: allows to add an addition to the file name (to distinguish
%different versions)
%
%output:
%structure MUThresh, one entry per channel
%channels are organized according to position on the probe (bottom of probe
%is channel 1)
%spikes marked as duplicate are removed
%fields:
%detChSort: sorted detection channel (same as position in MUThresh_
%detCh: original detection channel
%spiketime: spikes relative to the chosen event for every trial
%Nspk: total number of spikes in a trial
%baseNspike: number of spikes in baseline period
%baseFrate: firing rate in baseline period
%stimNspike: number of spikes in stimulus period
%stimFrate: firing rate in stimulus period
%
%structure MUinfo
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

        sidx=find(MUspkMerge.detChSort==chidx(u)  & ...
            MUspkMerge.spktimes>eTime-baseSample & MUspkMerge.spktimes<eTime+stimSample);
        %spkTimesTmp=(MUspkMerge.spktimes(sidx)-eTime)/id.sampleFreq*1000;
        %flagDupTmp=MUspkMerge.flagDuplicate(sidx);
        
        MUThresh(u).spktimes{i}=(MUspkMerge.spktimes(sidx)-eTime)/id.sampleFreq*1000;
        MUThresh(u).Nspk(i)=length(sidx);
        
        bidx=find(MUspkMerge.detChSort==chidx(u)  & ...
            MUspkMerge.spktimes>eTime-baseSample & MUspkMerge.spktimes<eTime);
        MUThresh(u).baseNspk(i)=length(bidx);
        MUThresh(u).baseFrate(i)=length(bidx)/baseTime;
       
        stidx=find(MUspkMerge.detChSort==chidx(u) & ...
            MUspkMerge.spktimes>eTime & MUspkMerge.spktimes<eTime+stimSample);
        MUThresh(u).stimNspk(i)=length(stidx);
        MUThresh(u).stimFrate(i)=length(stidx)/stimTime;
    end
end

%compute means and standard error per condition (baseline corrected)
condId=unique(trialInfo.triallist);
for i=1:length(condId)
    %find trials
    idx=find(trialInfo.triallist==condId(i));

    for u=1:length(chidx)
        MUThresh(u).avgRateCond(i)=mean(MUThresh(u).stimFrate(idx)-MUThresh(u).baseFrate(idx));
        MUThresh(u).semRateCond(i)=std(MUThresh(u).stimFrate(idx)-MUThresh(u).baseFrate(idx))/sqrt(length(idx));
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

trialExclude=MUThreshFlagOutlier(MUThresh,MUThreshInfo,0);
MUThreshInfo.trialExclude = trialExclude;

if ~isempty(varargin)
    save([basename '_p' num2str(probeId) '_MUThreshTrial_' varargin{1} '.mat'],'MUThresh','MUThreshInfo');
else
    save([basename '_p' num2str(probeId) '_MUThreshTrial.mat'],'MUThresh','MUThreshInfo');
end



