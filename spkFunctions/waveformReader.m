function [wvData,wvAvg,wvError]=waveformReader(spkpath,sortpath,animalID,unitID,expID,probeID,plotUnitNr,percSpikes,plotData)

%read waveforms of a selected unit and return
% input parameters:
% spkpath (character) - base path to folder containing SpikeFiles subfolder
% sortpath (character) - base path to folder containing spkSort and id files
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)
% plotUnitNr: nr of neuron to plot
% percSpikes: percentage of spikes to be plotted (randomly selected)
% plotData: produce figures? Note: assumes intan conversion from raw signal to
% uV
%
% output parameters:
% wvData: raw waveform timecourses
% wvAvg (optional): average timecourse
% wvError (optional): SEM around timecourse


expname=[animalID '_u' unitID '_' expID];


%load sort file
load(fullfile(sortpath,animalID,expname,[expname '_p' num2str(probeID) '_spkSort.mat'])); 

%load id file (for job edges)
load(fullfile(sortpath,animalID,expname,[expname '_id.mat'])); 


%get timestamps and detection channel for unit
ts=spkSort.spktimes(spkSort.unitid==plotUnitNr);
chidx=spkSort.detCh(spkSort.unitid==plotUnitNr);

%determine how many spikes should be plotted
nrSpk=floor(length(ts)*percSpikes/100);


%randomly subselect if necessary
if percSpikes<100
    ridx=randperm(length(ts),nrSpk);
    ts=ts(ridx);
    chidx=chidx(ridx);
end

%now read waveforms - map timestamps onto job files
tsbin=discretize(ts,id.extractSpikes.jobEdges{probeID}); 
tsunique=unique(tsbin);

%need loop to open files and read waveforms - first loop over job files
count=1;
for i=1:length(tsunique)
    load(fullfile(spkpath,animalID,expname,'SpikeFiles',...
        [expname '_j' num2str(tsunique(i)-1) '_p' num2str(probeID) '_spike'])); %generates spikeData

    %within job files loop over waveforms
    tsidx=find(tsbin==tsunique(i));
    for j=1:length(tsidx)
        widx=find(spikeData(chidx(tsidx(j))).spikeTimes==ts(tsidx(j))); %find which waveform

        wv=spikeData(chidx(tsidx(j))).Wvfrms(widx,:,1); %get waveform
        wvData(count,:)=squeeze(wv);

        count=count+1;
    end
end

wvAvg=mean(wvData,1);
wvError=std(wvData,0,1)/sqrt(size(wvData,1));


%plots selected?
if plotData==1
    %test whether id file contains settings
    if isfield(id.extractSpikes,'settings')
        spikeSamples=id.extractSpikes.settings{probeID}.spikeSamples;
    else
        %use the settings from the last spike file (still in memory)
        spikeSamples=settings.spikeSamples;
    end
        
    timevec=[-spikeSamples(1):spikeSamples(2)]/id.sampleFreq*1000;

    figure
    plot(timevec,wvData'*0.195);
    xlabel('time (ms)')
    ylabel('Amplitude (uV)')

    figure
    errorbar(timevec,wvAvg*0.195,wvError*0.195)
    xlabel('time (ms)')
    ylabel('Amplitude (uV)')

end