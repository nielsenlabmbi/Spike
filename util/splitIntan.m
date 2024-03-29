function splitIntan(expFolder,animalID,expID,probeID,name)
% split files that were sorted as one merged file
%
% input parameters:
% expFolder - base folder for experiments (string)
% animalID - animal ID (string)
% expID - experiment ID of the merged file (string; unit defaults to uMMM)
% probeID - number of probe (number)
% name - name or initials of person running the script (for bookkeeping)


%load spksort file
expname=[animalID '_uMMM_' expID];
spkSortIn=load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_spkSort.mat'])); %generates spkSort

%load merge info
load(fullfile(expFolder,animalID,expname,[expname '_mergeInfo.mat'])); %generates mergeInfo

%load id file for merged file - we need this for the extract offset and
%other
%info
idIn=load(fullfile(expFolder,animalID,expname,[expname '_id.mat'])); 

%compute edges of files (in samples)
%filesize - bytes per file merged
idx=[1:2:length(mergeInfo.filesize)*2-1];
blockLength(idx)=mergeInfo.filesize; %need to insert the buffer in between
%add buffer
for f=1:length(mergeInfo.filesize)-1
    blockLength(f*2)=mergeInfo.nChannel*mergeInfo.bufferLength*2; %2 b/c conversion to bytes
end
%now turn into edges in bytes
fileEdge=[0 cumsum(blockLength)];
%bytes to samples
fileEdge=fileEdge/(2*mergeInfo.nChannel);

%discretize samples in spkSort into the files
fileID=discretize(spkSortIn.spkSort.spktimes,fileEdge,'IncludedEdge','right');

%generate output
dataVar=fieldnames(spkSortIn.spkSort);
for f=1:length(mergeInfo.files)
    tf = fileID==(f*2)-1;
    
    spkSort=struct;

    for n=1:length(dataVar)
        if length(spkSortIn.spkSort.(dataVar{n}))==length(tf)
            spkSort.(dataVar{n})=spkSortIn.spkSort.(dataVar{n})(tf);
        else %just copy all of the info fields
            spkSort.(dataVar{n})=spkSortIn.spkSort.(dataVar{n});
        end
    end

    %add raw
    %spkSort.spktimesMerge=spkSort.spktimes;

    %need to correct the detection times (start is in odd bin of fileEdge)
    spkSort.spktimes=spkSort.spktimes-fileEdge(f*2-1);

    %usually we discard the first and last spikes for an amplifier file as
    %part of processing the first and last job
    %here, merging changes how jobs are organized, and spikes are not
    %dropped for the parts where files are connected
    offsetSamples=idIn.id.extractSpikes.settings{probeID}.offsetSamples;
    endFile=mergeInfo.filesize(f)/(mergeInfo.nChannel*2);
    tf = spkSort.spktimes>floor(offsetSamples/2);
    tf = tf & spkSort.spktimes<endFile-floor(offsetSamples/2);

    for n=1:length(dataVar)
        if length(spkSort.(dataVar{n}))==length(tf)
            spkSort.(dataVar{n})=spkSort.(dataVar{n})(tf);
        end
    end

    %save sort file
    outname=[animalID '_' mergeInfo.files{f}];
    outfile=fullfile(expFolder,animalID,outname,[outname '_p' num2str(probeID) '_' num2str(f) '_spkSort.mat']);
    %check first whether file exist - if yes ask
    if exist(outfile,'file')
        answer=questdlg('At least one spkSort file already exists and will be overwritten. Proceed?',...
            'Warning','Yes','Cancel','Yes');
        if strcmp(answer,'Cancel')
            return;
        end
    end
    save(outfile,'spkSort');

    %either generate new id files or add info to them
    idname=fullfile(expFolder,animalID,outname,[outname '_id.mat']);
    if ~exist(idname,'file')
        %generate id, only with the minimal information (since the rest was
        %not done on this file individually)
        id.exptId=outname;
        id.probes=idIn.id.probes;
        id.isBR=idIn.id.isBR;
        id.sampleFreq=idIn.id.sampleFreq;
    else
        load(idname); %generates id
    end
    id.spikeSort.name{probeID}=idIn.id.spikeSort.name{probeID};
    id.spikeSort.date{probeID}=idIn.id.spikeSort.date{probeID};
    id.spikeSort.NSingleUnit(probeID)=idIn.id.spikeSort.NSingleUnit(probeID);
    id.spikeSort.NMultiUnit(probeID)=idIn.id.spikeSort.NMultiUnit(probeID);

    %add info about merging and splitting (has not happened yet as the id file
    %is not touched/generated during merge)
    id.mergeInfo.processedMerged(probeID)=1;
    id.mergeInfo.splitSortFiles.name{probeID}=name;
    id.mergeInfo.splitSortFiles.date{probeID}=date;

    save(idname,'id');

end








%% save merge info file
save(fullfile(expFolder,animalID,outname,[outname '_mergeInfo']),'mergeInfo');