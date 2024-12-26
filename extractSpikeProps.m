function extractSpikeProps(expFolder,animalID,unitID,expID,probeID,name,copyToZ,MUflag,jobID,varargin)
% extract properties for spikes (one job file at a time)
% input parameters:
% expFolder - base folder for experiments (string)
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe ID (number)
% name - name or initials of person running the script (for bookkeeping)
% copyToZ - copy id file to Z?
% MUflag - 0 normal, 1 use MUspike files
% jobID - job ID of raw spike file to process (number)
% varargin -  2 options: 'id'+id (avoids loading the id file for batch processing); 
%             'tSuffix'+suffix for Spikes directory (to allow multiple
%             verions)
%             
%
% output parameters:
% structure spk with fields (each is a vector/matrix with entries for each
% spike)
% - EnAll: Energy for detection plus surrounding channels
% - EnDet: Energy for detection channel only (also contained in EnAll, this
% is to speed up later computations that only need the detection channel)
% - AmpMinAll, AmpMinDet: Minimum amplitude
% - AmxMaxAll, AmpMaxDet: Maximum amplitude for peak after minimum
% - AmxMaxBeforeAll, AmpMaxBeforeDet: Maximum amplitude for peak before
% minimum
% - PksAll, PksDet: Max-Min using maximum after minimum
% - PksBeforeAll, PksBeforeDet: Max-Min using maximum before minimum
% - WidthAll, WidthDet: Time of max - time of min
% - WidthIDet: Time of max - time of min, computed based on interpolation
% of spike waveforms
% - chListAll: List of channels included
% - comXMinDet, comZMinDet: center of mass for x and z, based on minimum
% - comXEnDet,comZEnDet: center of mass for x and z, based on energy
% - comZMinDetShaft, comZEnDetShaft: center of mass for z, organized to
% split shafts (based on detection channel)
% - spkTimesDet: spike times for detection channels
% - detCh: id of detection channel
% - detChSort: id of detection channel, sorted according to Z and shank
% - NDuplicate: for each spike, number of channels (across the entire probe) with a simultaneously
% detected spike
% - flagDuplicate: detection of one spike as separate events on
% different channels within a set radius; number indicates that the same spike was detected as
% a duplicate multiple times (spike with highest energy is taken to be the
% origin of the event)
% - duplicateMxCh: for events flagged as duplicate, the channel that was
% marked as the maximum/winning event
% - duplicateMxIdx: index into spkTimesDet etc for the winning event
%
% generates extractSpkProp file with info


%% deal with varargin
loadId=1;
tSuffix='';
spkFolder='SpikeFiles';

if ~isempty(varargin)
    idx=find(strcmp(varargin,'id'));
    if ~isempty(idx)
        loadId=0;
        id=varargin{idx+1};
    end
    idx=find(strcmp(varargin,'tSuffix'));
    if ~isempty(idx)
        tSuffix=varargin{idx+1};
        spkFolder=[spkFolder '_' tSuffix];
    end
end


%we need the id file for probe settings
expname=[animalID '_u' unitID '_' expID];
if loadId==1
    load(fullfile(expFolder,animalID,expname,[expname '_id'])); %generates id
end

%compute total channel number
nChannels=sum([id.probes.nChannels]);

%parameter choices
spkWindow=[-4 6]; %determines over how many datapoints we're looking for the minimum
spkTol=5; %window over which threshold crossings are considered duplicates/artefacts
spkInterp=0.1; %interpolation factor for width calculation
spkWindowI=15; %window in which to find minimum in interpolated data (scaling up spkWindow ends up too large)

%open matfile with spike data
%generates spikeData and settings
if MUflag==0
    fileIn=fullfile(expFolder,animalID,expname,spkFolder,[expname  '_j' num2str(jobID) '_p' num2str(probeID) '_spike']);

else
    fileIn=fullfile(expFolder,animalID,expname,spkFolder,[expname  '_j' num2str(jobID) '_p' num2str(probeID) '_MUspike']);
end
load(fileIn);

%samples per spike waveform
spikeSamples=settings.spikeSamples;
Nsample=sum(spikeSamples)+1;

%% go through spike file, compute properties for each spike, collected in arrays

spkCount=1;
spk=struct;

%sort channels according to z and shank
probeOrg=[id.probes(probeID).z id.probes(probeID).shaft id.probes(probeID).channels+1];
probeSort=sortrows(probeOrg,[2 1]);

%maximum z (for shaft separation)
maxZ=max(id.probes(probeID).z);

for i=1:id.probes(probeID).nChannels
     
    %only go through channels that actually have spikes
    if length(spikeData(i).spikeTimes)>1 | ~isnan(spikeData(i).spikeTimes)
        
        %base parameters
        Nspikes=length(spikeData(i).spikeTimes);
        Nch=length(spikeData(i).channelIds);
        
        %compute energy: sum squares across time
        En=squeeze(sqrt(sum(spikeData(i).Wvfrms.^2,2)));
        if Nspikes==1 %need to flip the vector in this case
            En=En';
        end
        
     
        %compute minimum and maximum
        TimeMin=zeros(Nspikes,Nch);
        AmpMin=zeros(Nspikes,Nch);
        TimeMax=zeros(Nspikes,Nch);
        AmpMax=zeros(Nspikes,Nch);
        TimeMaxB=zeros(Nspikes,Nch);
        AmpMaxB=zeros(Nspikes,Nch);
        WidthI=zeros(1,Nspikes);

        for s=1:Nspikes
            for c=1:Nch
                %find local minima in window around detection timepoint
                %organization of Wvfrms is spikes x times x channel
                tf=islocalmin(spikeData(i).Wvfrms(s,spikeSamples(1)+spkWindow(1):spikeSamples(1)+spkWindow(2),c),2,'FlatSelection','first');
                %to catch case without local min - set to detection event
                if sum(tf)==0
                    tf(abs(spkWindow(1))+2)=1; %spikeSamples(1) is 1 to the left of the minimum
                end
                %convert to indices into full time vector
                minIdx=find(tf==1)+spikeSamples(1)+spkWindow(1)-1;
                %find the one closest to the detection timepoint
                [~,minIdx2]=min(abs(minIdx-spikeSamples(1)-1));
                TimeMin(s,c)=minIdx(minIdx2);
                AmpMin(s,c)=spikeData(i).Wvfrms(s,TimeMin(s,c),c); %spikes x channels
            
                
                %find maxima - before and after
                tf=islocalmax(spikeData(i).Wvfrms(s,:,c),2,'FlatSelection','first');
                maxIdx=find(tf==1);

                %maximum before (last value of the ones that are smaller
                %than the detection time point)
                maxIdxB=find(maxIdx<spikeSamples(1)+1,1,'last');
                if ~isempty(maxIdxB)
                    TimeMaxB(s,c)=maxIdx(maxIdxB);
                else
                    TimeMaxB(s,c)=1;
                end

                %similar for after
                maxIdxA=find(maxIdx>spikeSamples(1)+1,1,'first');        
                if ~isempty(maxIdxA)
                    TimeMax(s,c)=maxIdx(maxIdxA);
                else
                    TimeMax(s,c)=Nsample;
                end
                
                AmpMax(s,c)=spikeData(i).Wvfrms(s,TimeMax(s,c),c);
                AmpMaxB(s,c)=spikeData(i).Wvfrms(s,TimeMaxB(s,c),c);
            end %for channel

            %compute width based on interpolating the spikes (detection
            %channel only)
            WvTmp=griddedInterpolant(squeeze(spikeData(i).Wvfrms(s,:,1)),'cubic');
            WvInterp=WvTmp([1:spkInterp:Nsample+1]);

            %minimum - search around the location of the minimum in the
            %non-interpolated version to avoid mismatches because of the
            %interpolation
            tf=islocalmin(WvInterp,'FlatSelection','first');
            minIdx=find(tf==1);
            TimeMinOrig=(TimeMin(s,1)-1)/spkInterp+1;
            [minVal,minIdx2]=min(abs(minIdx-TimeMinOrig));
            if isempty(minVal) || minVal>spkWindowI
                TimeMinI=TimeMinOrig;
            else
                TimeMinI=minIdx(minIdx2);
            end

            %maximum - again around the timepoint of the maximum detected
            %without interpolation
            tf=islocalmax(WvInterp,'FlatSelection','first');
            maxIdx=find(tf==1);
            TimeMaxOrig=(TimeMax(s,1)-1)/spkInterp+1;
            [maxVal,maxIdx2]=min(abs(maxIdx-TimeMaxOrig));
            if isempty(maxVal) || maxVal>spkWindowI
                TimeMaxI=TimeMaxOrig;
            else
                TimeMaxI=maxIdx(maxIdx2);
            end

            WidthI(s)=(TimeMaxI-TimeMinI)*spkInterp;

        end %for spikes
            
        %peak to peak value
        Pks = AmpMax-AmpMin;
        PksB = AmpMaxB-AmpMin;
        
        %width
        Width=TimeMax-TimeMin;

        %compute center of mass using minimum and energy, using the coordinates of the
        %channels
        %make sure to account for bad channels that are set to NaN
        chId=spikeData(i).channelIds;
        xpos=id.probes(probeID).x(chId);
        zpos=id.probes(probeID).z(chId);
        
        comXMin=sum(abs(AmpMin).*xpos',2,'omitnan')./sum(abs(AmpMin),2,'omitnan');
        comZMin=sum(abs(AmpMin).*zpos',2,'omitnan')./sum(abs(AmpMin),2,'omitnan');
        
        comXEn=sum(En.*xpos',2,'omitnan')./sum(En,2,'omitnan');
        comZEn=sum(En.*zpos',2,'omitnan')./sum(En,2,'omitnan');
        
        comZEnShaft=comZEn+(maxZ+100)*(id.probes(probeID).shaft(i)-1); %buffer of 100um between shafts
        comZMinShaft=comZMin+(maxZ+100)*(id.probes(probeID).shaft(i)-1);
        
        %since we are reorganizing into a structure per spike,
        %there are a couple of things we need to copy from
        %spikeData
        %time stamps
        spkTimes=spikeData(i).spikeTimes;
        
        %detection channel on probe        
        detChIdx=repmat(i,size(spkTimes));

        %detection channel sorted according to z
        detChIdxSort=find(probeSort(:,3)==i);
        detChIdxSort=repmat(detChIdxSort,size(spkTimes));
        
        %channel list 
        chList=repmat(spikeData(i).channelIds',Nspikes,1);
        
        %output
        %entries with spikes x channels - for all of these, we also
        %save the detection channel separately
        spk.EnAll(spkCount:spkCount+Nspikes-1)=num2cell(En,2); %number of channels might change
        spk.EnDet(spkCount:spkCount+Nspikes-1)=En(:,1);
        
        spk.AmpMinAll(spkCount:spkCount+Nspikes-1)=num2cell(AmpMin,2);
        spk.AmpMinDet(spkCount:spkCount+Nspikes-1)=AmpMin(:,1);
        
        spk.AmpMaxAll(spkCount:spkCount+Nspikes-1)=num2cell(AmpMax,2);
        spk.AmpMaxDet(spkCount:spkCount+Nspikes-1)=AmpMax(:,1);

        spk.AmpMaxBeforeAll(spkCount:spkCount+Nspikes-1)=num2cell(AmpMaxB,2);
        spk.AmpMaxBeforeDet(spkCount:spkCount+Nspikes-1)=AmpMaxB(:,1);
        
        spk.PksAll(spkCount:spkCount+Nspikes-1)=num2cell(Pks,2);
        spk.PksDet(spkCount:spkCount+Nspikes-1)=Pks(:,1);
        
        spk.PksBeforeAll(spkCount:spkCount+Nspikes-1)=num2cell(PksB,2);
        spk.PksBeforeDet(spkCount:spkCount+Nspikes-1)=PksB(:,1);
        
        spk.WidthAll(spkCount:spkCount+Nspikes-1)=num2cell(Width,2);
        spk.WidthDet(spkCount:spkCount+Nspikes-1)=Width(:,1);
        
        spk.WidthIDet(spkCount:spkCount+Nspikes-1)=WidthI;
        
        spk.chListAll(spkCount:spkCount+Nspikes-1)=num2cell(chList,2);
        
        %entries wit spikes x 1
        spk.comXMinDet(spkCount:spkCount+Nspikes-1)=comXMin;
        spk.comZMinDet(spkCount:spkCount+Nspikes-1)=comZMin;
        spk.comZMinShaftDet(spkCount:spkCount+Nspikes-1)=comZMinShaft;
        spk.comXEnDet(spkCount:spkCount+Nspikes-1)=comXEn;
        spk.comZEnDet(spkCount:spkCount+Nspikes-1)=comZEn;
        spk.comZEnShaftDet(spkCount:spkCount+Nspikes-1)=comZEnShaft;
        spk.spkTimesDet(spkCount:spkCount+Nspikes-1)=spkTimes;
        spk.detCh(spkCount:spkCount+Nspikes-1)=detChIdx;
        spk.detChSort(spkCount:spkCount+Nspikes-1)=detChIdxSort;
        
        spkCount=spkCount+Nspikes;
        
        
    end %if no spikes
end %for ch

if isfield(spk,'spkTimesDet')
    % now go through and flag potential duplicates - we do this at the spike
    %level (only flag events that are separately detected on multiple channels;
    %different from the amplitude extraction at neighboring channels
    %happening above)
    spk.flagDuplicate=zeros(size(spk.spkTimesDet));
    spk.NDuplicate=zeros(size(spk.spkTimesDet));
    spk.duplicateMxCh=zeros(size(spk.spkTimesDet));
    spk.duplicateMxIdx=zeros(size(spk.spkTimesDet));
    timesIdx=[1:length(spk.spkTimesDet)];
    
    %for each threshold crossing, find out how many other threshold
    %crossings occur within the tolerance window
    tol=spkTol/max(spk.spkTimesDet); %uniquetol scales by maximum
    
    [~,~,idxB]=uniquetol(spk.spkTimesDet,tol); %get unique events plus/minus tolerance, idxB is an index into spkTimesDet
    countDup=accumarray(idxB,1); %count how often each duplicate shows up in spkTimesDet
    spk.NDuplicate=countDup(idxB); %build vector that gives N for each event
    spk.NDuplicate=spk.NDuplicate';
    
    %in addition to figuring out whether an event is detected on multiple
    %channels, get the one that has the maximal amplitude and flag the rest
    %as potential duplicates; this only will be applied to channels that
    %are close enough
    
    for i=1:id.probes(probeID).nChannels
        if sum(spk.detCh==i)>0
            
            %find spikes at detection channel, get their energy and position in
            %spkTimes vector
            detTimes=spk.spkTimesDet(spk.detCh==i);
            timesIdxDet=timesIdx(spk.detCh==i);
            enDet=spk.EnDet(spk.detCh==i);
            
            %spread out each event over neighboring samples (to give interval
            %for detection)
            detTimesConv=detTimes+[spkWindow(1):spkWindow(2)]'; 
            detTimesConv=detTimesConv(:); %this contains detection times plus the interval around them
            
            %also need an index for grouping later - this indexes into the
            %triggering events, with repeating numbers for the same event
            detIdxConv=repmat([1:length(detTimes)],spkWindow(2)-spkWindow(1)+1,1);
            detIdxConv=detIdxConv(:);
            
            %get the other channels
            chList=spikeData(i).channelIds; 
            chList=chList(2:end); %need to remove center channel - chList only contains the other channels
            
            %find overlapping events - events in other channels close to
            %events in the channel of interest
            tf=ismember(spk.detCh,chList) & ismember(spk.spkTimesDet,detTimesConv);
            
            %only continue if there actually are any
            if sum(tf)>0
                %get energy of overlapping events
                enCh=spk.EnDet(tf);
                
                %get their index in the spkTimesDet vector
                timesIdxCh=timesIdx(tf);
                
                %get their channel
                detCh=spk.detCh(tf);

                %get index of triggering event to use as grouping par
                %we're doing it this way because there are different
                %numbers of overlapping events for every spike on the detection channel
                [~,locB]=ismember(spk.spkTimesDet(tf),detTimesConv); %need tf here to restrict to events in the right channels
                idxCh=detIdxConv(locB); %id of trigger event
                
                %generate one big array for accumarray (enDet is the
                %channel of interest)
                enAll=[enDet(unique(idxCh)) enCh]'; %we only want the original events that had overlaps
                idxAll=[unique(idxCh);idxCh]; %group by trigger event id
                timesIdxAll=[timesIdxDet(unique(idxCh)) timesIdxCh]'; %indexes into spkTimesDet ultimately
                detChAll=[repmat(i,1,length(unique(idxCh))) detCh]';
                
                %get maximum per group
                maxData=accumarray(idxAll,enAll,[],@max); %in maxdata, a group shows up in the row corresponding to its index value
                
                %need one loop here to figure out who is generating the
                %maximum
                maxCh=zeros(size(maxData));
                maxTime=zeros(size(maxData));
                for m=unique(idxCh)' %this only pulls maxima for existing events
                    mx=find(enAll==maxData(m),1,'first');
                    maxCh(m)=detChAll(mx);
                    maxTime(m)=timesIdxAll(mx);
                end
                maxChD = maxCh(idxAll);
                maxTimeD = maxTime(idxAll);

                %flag which of the overlapping events are not the maximum (i.e.
                %duplicates)
                %we only eliminate the spikes on the current detection
                %channel; otherwise a channel in the middle can propagate
                %events over a much larger distance than intended
                %previously: flagD = enAll~=maxData(idxAll)
                
                flagD = (enAll~=maxData(idxAll) & detChAll==i);

                if sum(flagD)>0
                    %put back into large vector
                    idxD=timesIdxAll(flagD);
                    
                    maxChD = maxChD(flagD);
                    maxTimeD = maxTimeD(flagD);
                    
                    spk.flagDuplicate(idxD)=1; 
                    spk.duplicateMxCh(idxD)=maxChD;
                    spk.duplicateMxIdx(idxD)=maxTimeD;
                end
            end
        end %if sum
    end %for ch
end
%add corresponding spike file name
spk.fileIn=fileIn;


if MUflag==0
    outname=fullfile(expFolder,animalID,expname,spkFolder,[expname  '_j' num2str(jobID) '_p' num2str(probeID) '_spkinfo']);
else
    outname=fullfile(expFolder,animalID,expname,spkFolder,[expname  '_j' num2str(jobID) '_p' num2str(probeID) '_MUspkinfo']);
end

save(outname,'spk','-v7.3','-nocompression');


%generate extractSpkProp file for job 0
if jobID==0
    extractSpkProp.date=date;
    extractSpkProp.name=name;
    extractSpkProp.probeNr=probeID;
    extractSpkProp.exptId=expname;
    extractSpkProp.suffix=tSuffix;
    extractSpkProp.MUflag=MUflag;
    extractSpkProp.spkWindow=spkWindow;
    extractSpkProp.spkTol=spkTol;
    extractSpkProp.spkInterp=spkInterp;
    extractSpkProp.spkWindowI=spkWindowI;

    if MUflag==0
        infoname=[expname '_p' num2str(probeID) '_extractSpkProp'];
    else
        infoname=[expname '_p' num2str(probeID) '_MUextractSpkProp'];
    end
    if ~isempty(tSuffix)
        infoname=[infoname '_' tSuffix '.mat'];
    else
        infoname=[infoname '.mat'];
    end

    save(fullfile(expFolder,animalID,expname,infoname),'extractSpkProp');
    
    if copyToZ==1
        zbase='Z:\EphysNew\processedSpikes';
        save(fullfile(zbase,animalID,expname,infoname),'extractSpkProp'); 
    end

end

disp(['extractSpikeProps job ID ' num2str(jobID) ' done.'])







