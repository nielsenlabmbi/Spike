function extractSpikeProps(expFolder,animalID,unitID,expID,jobID)

% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% expFolder - base folder for experiments (string)
% jobID - job ID of raw spike file to process (number)

%we need the id file for probe settings
expname=[animalID '_u' unitID '_' expID];
load(fullfile(expFolder,animalID,expname,[expname '_id'])); %generates id

%compute total channel number
nChannels=sum([id.probes.nChannels]);

%% go through spike file, compute properties for each spike, collected in arrays

%open matfile with spike data - note if too large load 1 probe only? 1
%channel at a time is too slow
%generates spikeData and settings
load(fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(jobID) '_spike'])); 
%samples per spike waveform
spikeSamples=settings.spikeSamples;

spkCount=1;
for p=1:length(id.probes)
    for i=1:id.probes(p).nChannels
        
        offsetCh=sum([id.probes(1:p-1).nChannels]); %0 for p=1
             
        %only go through channels that actually have spikes
        if length(spikeData(i+offsetCh).spikeTimes)>1 || ~isnan(spikeData(i+offsetCh).spikeTimes)
            
            %base parameters
            Nspikes=length(spikeData(i+offsetCh).spikeTimes);
            Nch=length(spikeData(i+offsetCh).channelIds);
            
            %compute energy: sum squares across time
            En=squeeze(sum(spikeData(i+offsetCh).Wvfrms.^2,2));
            if Nspikes==1 %need to flip the vector in this case
                En=En';
            end
            
            %compute minimum and maximum
            %we will base these on the derivative - we want a minimum
            %that is the local minimum around the minimum time at the
            %detection channel (for which the minimum will be at
            %spikeSamples+1), and a maximum that is the first maximum
            %after the positive going part of the waveform
            %both are more correctly detected based on the slope rather
            %than using min/max over a set interval
            %we can't use the minimum at the detection time point
            %because the minimum occurs at different time points on
            %different channels
            
            %sign of slope
            Wv1Der=sign(diff(spikeData(i+offsetCh).Wvfrms,1,2));
            
            %2nd derivative of the sign
            Wv2Der=diff(Wv1Der,1,2);
            
            %minimum: maximum in 2nd derivative around the spikeSamples+1
            %using +-5 data points for now
            [minDer,TimeMin]=max(Wv2Der(:,spikeSamples-4:spikeSamples+6,:),[],2);
            
            %because of how diff works, max of derivative is one ahead of true min, also
            %add time offset relative to original waveform
            TimeMin=squeeze(TimeMin)+spikeSamples-4; %spikes x channel
            
            %if there is no local minimum, use the minimum of the
            %detection channel
            TimeMin(squeeze(minDer)==0)=spikeSamples+1;
            
            %now get minimum values
            %need to turn TimeMin into the correct index first
            spkIdx=repmat([1:Nspikes]',1,Nch);
            chIdx=repmat([1:Nch],Nspikes,1);
            if Nspikes==1 %need to flip the time vector in this case
                TimeMin=TimeMin';
            end
            mxIdx=sub2ind([Nspikes 2*spikeSamples+1 Nch],spkIdx,TimeMin,chIdx);
            AmpMin=spikeData(i+offsetCh).Wvfrms(mxIdx); %spikes x channels
            
            %do the same for the maximum (minimum in 2nd derivative)
            %searching after occurence of the minimum for detection
            %channel
            %min automatically returns the index to the first occurence if
            %the minimum occurs more than once
            [~,TimeMax]=min(Wv2Der(:,spikeSamples+1:end,:),[],2);
            TimeMax=squeeze(TimeMax)+spikeSamples+1; %spikes x channel
            
            %now get maximum values
            if Nspikes==1 %need to flip the time vector in this case
                TimeMax=TimeMax';
            end
            mxIdx=sub2ind([Nspikes 2*spikeSamples+1 Nch],spkIdx,TimeMax,chIdx);
            AmpMax=spikeData(i+offsetCh).Wvfrms(mxIdx);
            
            %peak to peak value
            Pks = AmpMax-AmpMin;
            
            %width
            Width=TimeMax-TimeMin;
            
            %compute center of mass using minimum and energy, using the coordinates of the
            %channels; we won't separate by shank here since every
            %channel has a unique position, so it's not necessary
            xpos=id.probes(p).x(spikeData(i+offsetCh).channelIds);
            zpos=id.probes(p).z(spikeData(i+offsetCh).channelIds);
            
            comXMin=sum(abs(AmpMin).*xpos',2)./sum(abs(AmpMin),2);
            comZMin=sum(abs(AmpMin).*zpos',2)./sum(abs(AmpMin),2);
            
            comXEn=sum(En.*xpos',2)./sum(En,2);
            comZEn=sum(En.*zpos',2)./sum(En,2);
            
            
            %since we are reorganizing into a structure per spike,
            %there are a couple of things we need to copy from
            %spikeData
            %time stamps
            spkTimes=spikeData(i+offsetCh).spikeTimes;
            
            %detection channel
            detChIdx=repmat(i+offsetCh,size(spkTimes));
            
            %channel list
            chList=repmat(spikeData(i+offsetCh).channelIds',Nspikes,1);
            
            %output
            %entries with spikes x channels - for all of these, we also
            %save the detection channel separately
            spk.EnAll(spkCount:spkCount+Nspikes-1)=num2cell(En,2); %number of channels might change
            spk.EnDet(spkCount:spkCount+Nspikes-1)=En(:,1);
            
            spk.AmpMinAll(spkCount:spkCount+Nspikes-1)=num2cell(AmpMin,2); 
            spk.AmpMinDet(spkCount:spkCount+Nspikes-1)=AmpMin(:,1); 
            
            spk.AmpMaxAll(spkCount:spkCount+Nspikes-1)=num2cell(AmpMax,2);
            spk.AmpMaxDet(spkCount:spkCount+Nspikes-1)=AmpMax(:,1);
            
            spk.PksAll(spkCount:spkCount+Nspikes-1)=num2cell(Pks,2);
            spk.PksDet(spkCount:spkCount+Nspikes-1)=Pks(:,1);
            
            spk.WidthAll(spkCount:spkCount+Nspikes-1)=num2cell(Width,2); 
            spk.WidthDet(spkCount:spkCount+Nspikes-1)=Width(:,1); 
            
            spk.chListAll(spkCount:spkCount+Nspikes-1)=num2cell(chList,2);
            
            %entries wit spikes x 1
            spk.comXMinDet(spkCount:spkCount+Nspikes-1)=comXMin;
            spk.comZMinDet(spkCount:spkCount+Nspikes-1)=comZMin;
            spk.comXEnDet(spkCount:spkCount+Nspikes-1)=comXEn;
            spk.comZEnDet(spkCount:spkCount+Nspikes-1)=comZEn;
            spk.spkTimesDet(spkCount:spkCount+Nspikes-1)=spkTimes;
            spk.detCh(spkCount:spkCount+Nspikes-1)=detChIdx;
            
            
            spkCount=spkCount+Nspikes;
            
            
        end %if no spikes
        
    end %for ch
end %for probes

%generate some basic info
spk.dateProps=date;
spk.expname=expname;

outname=fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(jobID) '_spkinfo']);
save(outname,'spk');
disp(['extractSpikes job ID ' num2str(jobID) ' done.'])



            
            
            
     
            
  