function splitIdFile(fileBase,animalId,unitId,expId)
%split id files into separate setting files to be compatible with the
%processing pipeline allowing multiple versions
%we are ignoring all of the MU processing here as the id files often do not
%correctly reflect which threshold files were used to extract spikes;
%finally, no update of threshold file since that would just be name & date, 
%but update of spkSort file to be able to read it with sortGui
%check id file first to make sure information is actually correct

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
fnames=fieldnames(idIn.id);
for i=1:length(fnames)
    if ~any(strcmp(fnames{i},{'probes','isBR','exptId','sampleFreq'}))
        id.old.(fnames{i})=idIn.id.(fnames{i});
    end
end

save(fullfile(fileBase,animalId,expname,[expname '_id.mat']),'id');



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

        extractSpk.probeNr=i;

        %threshold - for spike files, there aren't any versions (only for
        %MU)
        extractSpk.exptId=idIn.id.exptId;
        extractSpk.threshold.filename='threshold';
        extractSpk.threshold.date=idIn.id.threshold.name{i};
        extractSpk.threshold.name=idIn.id.threshold.date{i};
        extractSpk.MUflag=0;
                   
        %save
        outname=[expname '_p' num2str(i) '_extractSpk'];
        save(fullfile(fileBase,animalId,expname,outname),'extractSpk');
        
    end
end


%% extractSpikeProps
if isfield(idIn.id,'extractSpikeProps')
    for i=1:length(idIn.id.extractSpikeProps.date)
        extractSpkProp=struct;
        extractSpkProp.date=idIn.id.extractSpikeProps.date{i};
        extractSpkProp.name=idIn.id.extractSpikeProps.name{i};
        extractSpkProp.probeNr=i;
        extractSpkProp.exptId=idIn.id.exptId;
        extractSpkProp.suffix=''; %again, no versions for spike files
        extractSpkProp.MUflag=0;
        %taking a guess at settings, but these basically were never changed
        extractSpkProp.settings.spkWindow=[-4 6];
        extractSpkProp.settings.spkTol=5;
        extractSpkProp.settings.spkWindowI=15;
        extractSpkProp.settings.spkInterp=0.1;
    
        %save
        outname=[expname '_p' num2str(i) '_extractSpkProp.mat'];
        save(fullfile(fileBase,animalId,expname,outname),'extractSpkProp');
        
    end
end


%% spkSort
if isfield(idIn.id,'spikeSort')
    for i=1:length(idIn.id.spikeSort.date)
        loadSort=0; %don't want to use return because the rest of the function should still execute

        spkSortName=fullfile(fileBase,animalId,expname,[expname '_p' num2str(i) '_spkSort.mat']);
        if isfile(spkSortName)
            loadSort=1;
        else
            sel=questdlg('No spkSort file found. Pick one?','spkSort file',...
                'Yes','Skip updating spkSort');
            if strcmp(sel,'Yes')
                [sortName,p] = uigetfile('*spkSort*.mat','Select sort file');
                spkSortName=fullfile(p,sortName);
                loadSort=1;
            end
        end

        if loadSort==1
            load(spkSortName);
            
            %double check they actually match
            updateSort=1;
            if ~strcmpi(spkSort.info.name,idIn.id.spikeSort.name{i}) || ...
                    datetime(spkSort.info.date)~=datetime(idIn.id.spikeSort.date{i})
                sel=questdlg('Mismatch between spkSort and id file info (name or date). Proceed?','Mismatched info',...
                    'Yes','No');
                if strcmp(sel,'No')
                    updateSort=0;
                end
            end

            if updateSort==1
                if ~isfield(spkSort.info,'percJobs')
                    spkSort.info.percJobs=100;
                end
                if ~isfield(spkSort.info,'jobsList')
                    spkSort.info.jobsList=[spkSort.info.jobStart:spkSort.info.jobStop];
                end
                spkSort.info.jobEdges=idIn.id.extractSpikes.jobEdges{i};
                spkSort.info.tSuffix='';

                %save
                save(spkSortName,'spkSort');
            end

        end
    end
end

%% done - warning dialog
warndlg('Updated all settings files. Please make sure to check files to make sure information is correct!');
