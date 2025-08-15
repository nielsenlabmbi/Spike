function SUTrialData(physpath,animal,unit,exp,probeId,eventType,eId,baseTime,stimTime,copyToZ,varargin)
%this function extracts spike times for identifed SU and MU spikes (noise
%and none are excluded) for each trial, also computes mean number of spikes
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
%copyToZ: copy file to Z
%varargin: allows addition of suffix to the spkSort file, as well as
%          addition of suffix to the output file (do not have to be the same)
%          indicate using 'spkSort' + suffix and 'out' + suffix
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
sortName=[basename '_p' num2str(probeId) '_spkSort'];
if ~isempty(varargin)
    idx=find(strcmp(varargin,'spkSort'));
    if ~isempty(idx)
        sortName=[sortName '_' varargin{idx+1}];
    end
end
load(sortName);

%get sampling rate - either from header file or from id file
headerName=[basename '_info.rhd'];
idName=[basename '_id.mat'];
if exist(headerName,'file')~=0
    header=read_Intan_Header(headerName);
    sampleFrq = header.sample_rate;
elseif exist(idName,'file')~=0
    load(idName);
    sampleFrq=id.sampleFreq;
else
    disp('file with sample rate missing')
    return;
end

%translate time windows into samples
baseSample=round(baseTime*sampleFrq);
stimSample=round(stimTime*sampleFrq);

%figure out which single units we want
unitIdx=find(strcmp(spkSort.unitinfo,'SU') | strcmp(spkSort.unitinfo,'MU'));

%find events
if strcmp(eventType,'id')
    eventIdx=find(trialInfo.eventId==eId);
else
    eventIdx=find(trialInfo.eventCh(:,eId)==1);
end

%loop through units and events to get data
for u=1:length(unitIdx)
    %base info
    SU(u).unitId=unitIdx(u);
    SU(u).unitClass=spkSort.unitinfo{unitIdx(u)};
    
    
    for i=1:length(eventIdx)
        eTime=trialInfo.eventTimes(eventIdx(i));

        sidx=find(spkSort.unitid==unitIdx(u) & spkSort.spktimes>eTime-baseSample & spkSort.spktimes<eTime+stimSample);
        SU(u).spktimes{i}=sort((spkSort.spktimes(sidx)-eTime)/sampleFrq*1000);
        SU(u).Nspk(i)=length(sidx);
        
        bidx=find(spkSort.unitid==unitIdx(u) & spkSort.spktimes>eTime-baseSample & spkSort.spktimes<eTime);
        SU(u).baseNspk(i)=length(bidx);
        SU(u).baseFrate(i)=length(bidx)/baseTime;
       
        stidx=find(spkSort.unitid==unitIdx(u) & spkSort.spktimes>eTime & spkSort.spktimes<eTime+stimSample);
        SU(u).stimNspk(i)=length(stidx);
        SU(u).stimFrate(i)=length(stidx)/stimTime;
    end
end

%compute means and standard error per condition (baseline corrected)
condId=unique(trialInfo.triallist);
for i=1:length(condId)
    %find trials
    idx=find(trialInfo.triallist==condId(i));

    for u=1:length(unitIdx)
        SU(u).avgRateCond(i)=mean(SU(u).stimFrate(idx)-SU(u).baseFrate(idx));
        SU(u).semRateCond(i)=std(SU(u).stimFrate(idx)-SU(u).baseFrate(idx))/sqrt(length(idx));
    end
end


%generate structure with general info, including condition info
SUinfo.eventType=eventType;
SUinfo.eventId=eId;
SUinfo.baseTime=baseTime;
SUinfo.stimTime=stimTime;
SUinfo.dom=trialInfo.dom;
SUinfo.domval=trialInfo.domval;
SUinfo.blankId=trialInfo.blankId;
SUinfo.triallist=trialInfo.triallist;
if ~isempty(varargin)
    SUinfo.tSuffix=varargin{1};
end
SUinfo.sortfile=sortName;
SUinfo.sortdate=spkSort.info.date;
SUinfo.sortname=spkSort.info.name;

%save
fname1=[basename '_p' num2str(probeId) '_SUTrial'];
if ~isempty(varargin)
    idx=find(strcmp(varargin,'out'));
    if ~isempty(idx)
        fname1=[fname1 '_' varargin{idx+1}];
    end
end

save(fname1,'SU','SUinfo');

if copyToZ==1
    fname2=fullfile('Z:\EphysNew\processedSpikes',animal,[animal '_u' unit '_' exp],...
            [animal '_u' unit '_' exp '_p' num2str(probeId) '_SUTrial']);
    if ~isempty(varargin)
        idx=find(strcmp(varargin,'out'));
        if ~isempty(idx)
            fname2=[fname2 '_' varargin{idx+1}];
        end
end
    save(fname2,'SU','SUinfo');
end
