%split id files into id and history
%function splitIdFile(fileBase,animalId,unitId,expId)

fileBase='c:/kristina/ephys';
animalId='feaz5';
unitId='000';
expId='002';

expname=[animalId '_u' unitId '_' expId];
fname=fullfile(fileBase,animalId,expname,[expname '_id.mat']);

idIn=load(fname); 

%% extract id part - probes, expitId, isBR, sampleFreq
id=struct;
id.probes=idIn.id.probes;
id.exptId=idIn.id.exptId;
id.isBR=idIn.id.isBR;
id.sampleFreq=idIn.id.sampleFreq;

%move rest to old to keep record
if isfield(idIn.id,'threshold')
    id.old.threshold=idIn.id.threshold;
end

if isfield(idIn.id,'extractSpikes')
    id.old.extractSpikes=idIn.id.extractSpikes;
end

if isfield(idIn.id,'extractSpikeProps')
    id.old.extractSpikeProps=idIn.id.extractSpikeProps;
end

if isfield(idIn.id,'spikeSort')
    id.old.spikeSort=idIn.id.spikeSort;
end

%save(fnameId,'id');

%% threshold

%% extractSpikes

if isfield(idIn.id,'extractSpikes')
    %separate files for multiple probes
    for i=1:length(idIn.id.extractSpikes.date)
        extractSpk=struct;
        if iscell(idIn.id.extractSpikes.jobStart) %there are old files that don't have these stored per probe
            extractSpk.jobStart=idIn.id.extractSpikes.jobStart{i};
            extractSpk.jobStop=idIn.id.extractSpikes.jobStop{i};
            extractSpk.jobEdges=idIn.id.extractSpikes.jobEdges{i};
        else
            extractSpk.jobStart=idIn.id.extractSpikes.jobStart;
            extractSpk.jobStop=idIn.id.extractSpikes.jobStop;
            extractSpk.jobEdges=idIn.id.extractSpikes.jobEdges;
        end

        extractSpk.date=idIn.id.extractSpikes.date{i};
        extractSpk.name=idIn.id.extractSpikes.name{i};
        if isfield(idIn.id.extractSpikes,'settings')
            extractSpk.settings=idIn.id.extractSpikes.settings{i};
        else %at least legacy flag exists, so copy that
            if length(idIn.id.extractSpikes.legacyFlag)>1
                extractSpk.settings.legacyFlag=idIn.id.extractSpikes.legacyFlag(i);
            else
                extractSpk.settings.legacyFlag=idIn.id.extractSpikes.legacyFlag;
            end
        end

        %threshold
        

        extractSpk.probeNr=idIn.id.extractSpikes.i;
    end
        

end


%% extractSpikeProps

%% spkSort
% add edges to spkSort
% add jobsList
