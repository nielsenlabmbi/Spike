function detChOldFormat(expFolder,animalID,unitID,expID,probeID,detChannelMethod,spikesIn,plotChecks,rawDataFolder,nrPlotCh,nrPlotEvent)
%function to compute the raw detection channel (i.e. the Intan recording channel)
%for the spike or Spike files generated with the original spike sorter
%NOTE: MAKE SURE RESULTS APPEAR REASONABLE BY PLOTTING DATA

%input parameters:
% expFolder - experiment folder
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe number to process (number)

%spkSortIn: spkSort input file

%detChannelMethod: method used to generate detection channels
%1: use detChSort of spkSort structure (=Properties(16,:) of spikes or
%spikes File); this only works if there were no bad channels
%2: use experiment file; this only works if experiments file has not been
%overwritten (for multiple probes)
%3: use Properties(17,:) and Properties(18,:) from spikes file
%4: set to all 1 (this is for tetrode data only!)

%spikesIn: only needed for method 3, spikes or Spikes file to be loaded
%(can be empty string otherwise)

%plotCheck: 1 enables plotting spikes to check channel assignemt

%the next 3 parameters are only needed when plotting, otherwise provide
%empty strings
%rawDataFolder: base name of folder containing amplifier file for raw data
%(e.g., z:\EphysNew); animal, experiment etc are added automatically
%nrPlotCh: number of channels to plot (which channels is randomly selected)
%nrPlotEvent: number of events to plot (again, randomly selected)


%% channel conversion code
%parse file name and figure out probe number
expname=[animalID '_u' unitID '_' expID];
load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_spkSort'])); %generates spkSort
nChIn=length(unique(spkSort.detChSort));


%load id file
load(fullfile(expFolder,animalID,expname,[expname '_id.mat'])); %generates id


%get probe type
probeType=id.probes(probeID).type;
nChProbe=id.probes(probeID).nChannels;

%get channel configuration used in makeExperimentFile
switch probeType
    case '64M' %this seems wrong
        CHs = [56,10,55,11,54,12,53,13,52,14,51,15,50,16,49,17,48,18,47,19,46,20,45,21,44,22,43,23,42,24,41,25,40,26,39,27,38,28,37,29,36,30,35,31,34,32,33,9,57,8,58,7,59,6,60,5,61,4,62,3,63,2,64,1];
    case '64D'
        CHs = [17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54];
    case '64F'
        CHs = [13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49] ;
end


switch detChannelMethod
    case 1 %use detChSort as is
        %this method will depend on no channels being marked as bad - test
        %whether all channels are present
        if nChIn~=nChProbe
            disp('Channels are missing from detChSort, cannot proceed with conversion!')
            return;
        end

        %all channels present, compute detCh and add to spkSort
        spkSort.detCh=CHs(spkSort.detChSort);

    case 2 %use experiment.mat
        %load experiment file
        load(fullfile(expFolder,animalID,expname,'experiment.mat'));

        %check whether length of exeriment.mat matches length of channels in
        %spkSort
        if nChIn~=sum(~BadCh)
            disp('Channels are missing from detChSort, cannot proceed with conversion!')
            return;
        end

        %also check if second probe (if applicable)
        if probeID==2
            minEx=find(BadCh==0,1,'first');
            maxEx=find(BadCh==0,1,'last');

            if minEx<=64 || maxEx<=64
                disp('experiment.mat is for other probe, cannot proceed with conversion!')
                return;
            end
        end

        %remove bad channels
        CHs=CHs(~BadCh);

        %all channels present, compute detCh and add to spkSort
        spkSort.detCh=CHs(spkSort.detChSort);

    case 3 %use properties 17 and 18
        %load spikes or Spikes file
        load(spikesIn); %generates Properties

        %check that Properties has the right info
        if size(Properties,1)<17
            disp('Wrong Properties format, cannot proceed with conversion!')
            return;
        end

        %may need to correct x (offset in same cases)
        if min(Properties(17,:))>10000
            Properties(17,:)=Properties(17,:)-1e11;
        end

        %channel map
        chid=unique(Properties(16,:));
        chmap=zeros(max(chid),1);
        for i=chid
            idx=find(Properties(16,:)==i,1,'first');
            chidx=find(id.probes(probeID).x==Properties(17,idx) & id.probes(probeID).z==Properties(18,idx));
            if ~isempty(chidx)
                chmap(i)=chidx;
            else
                disp('cannot match properties to probe layout!')
            end
        end

        %now use channel map to get channels
        spkSort.detCh=chmap(spkSort.detChSort);

    case 4
        spkSort.detCh=ones(size(spkSort.spktimes));

end

%% plot data to check assignment
if plotChecks==1
    %filter settings
    hp=250;
    lp=5000;
    [b1,a1]=butter(3,[hp/id.sampleFreq,lp/id.sampleFreq]*2,'bandpass');

    %determine which channels to plot - linearly subsample channel range
    if nrPlotCh>nChIn
        nrPlotCh=nChIn;
    end
    plotChIdx=randsample(unique(spkSort.detCh),nrPlotCh);
 
    figure;
    t=tiledlayout(nrPlotCh,nrPlotEvent);
    t.TileSpacing='compact';
    t.Padding='compact';

    fid = fopen(fullfile(rawDataFolder,animalID,expname,[expname '_amplifier.dat']),'r');
    
    for c=1:nrPlotCh
        %random events
        evIdx=spkSort.spktimes(spkSort.detCh==plotChIdx(c));
        if nrPlotEvent>length(evIdx)
            nrPlotEvent=length(evIdx);
        end
        plotEvIdx=randsample(evIdx,nrPlotEvent);

        for i=1:nrPlotEvent
            frewind(fid);
            startSample=plotEvIdx(i)-id.sampleFreq/2; %keeping this long for filtering
                        
            if probeID==2
                cc=plotChIdx(c)+id.probes(1).nChannels;
            else
                cc=plotChIdx(c);
            end

            fseek(fid,2*sum([id.probes.nChannels])*startSample+2*(cc-1),'bof');
            tc = fread(fid, id.sampleFreq, 'int16',2*(sum([id.probes.nChannels])-1));

            %filter
            datFilt=filter(b1,a1,tc);


            %get waveform
            px=nexttile;
            plot([-20:20]./id.sampleFreq*1000,datFilt(id.sampleFreq/2-20:id.sampleFreq/2+20))
            hold on
            xline(0,'r-')
            title(['Ch: ' num2str(plotChIdx(c))])
            

        end

    end

    fclose(fid);


end


%% save
%copy detChSort into separate field (since we'll overwrite it later), and
%remove it here
spkSort.detChSortOrig=spkSort.detChSort;
spkSort=rmfield(spkSort,'detChSort');

save(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_spkSort']),'spkSort');