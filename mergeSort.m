function mergeSort(expFolder,animalID,unitID,expID,probeID,name,copyToZ)
% merging of sort files
%
% input parameters:
% expFolder - base folder for experiments (string)
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe ID (number)
% name - name or initials of person running the script (for bookkeeping)
% copyToZ - copy id file to Z?

% output:
% one file
% structure spkSort with fields unitid, spktimes, detCh,detChSort

%load id file
expname=[animalID '_u' unitID '_' expID];
load(fullfile(expFolder,animalID,expname,[expname '_id']));


%find all spkSort files in SpikeFiles folder
sortFiles=dir(fullfile(expFolder,animalID,expname,'SpikeFiles',[expname  '_j*_p' num2str(probeID) '_spkSort.mat']));

%merge
spkSort=struct;
for i=1:length(sortFiles)
    %load data
    load(fullfile(sortFiles(i).folder,sortFiles(i).name)); %creates spkSortP

    dataVar=fieldnames(spkSortP);

    for n=1:length(dataVar)
        if i==1
            spkSort.(dataVar{n})=spkSortP.(dataVar{n});
        else
            if ~strcmp(dataVar{n},'info') 
                spkSort.(dataVar{n})=[spkSort.(dataVar{n}) spkSortP.(dataVar{n})];
            end
        end

    end

    clear spkSortP;
end

%determine number of units
unitIdx=unique(spkSort.unitid(spkSort.unitid>0));
nrUnits=length(unitIdx);

%add unit assignments (based on spkSort file corresponding to history)
spkSortGui=load(fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeID) '_partSpkSort'])); %generates field spkSort
spkSort.unitinfo=spkSortGui.spkSort.unitinfo;

%display warning if things don't match
if length(spkSort.unitinfo)~=nrUnits
    errordlg('Mismatch between nr of units in merged and partially sorted file!','Error');
    return;
end


%compute unit metrics (not Lratio, since it requires metrics we don't have
%access to right now)
minISI=0.0012;
minISIsample=round(minISI*id.sampleFreq);
for i=1:nrUnits
    %isi violations
    ts=sort(spkSort.spktimes(spkSort.unitid==i));
    diffTs=diff(ts); %difference between spikes in samples
    perISIv=sum(diffTs<minISIsample)/length(ts);
    ISIv=perISIv*100;

    spkSort.unitisi(i)=ISIv;
    
    %footprint
    chidx=spkSort.detCh(spkSort.unitid==i);
    spkSort.unitFP(i)=length(unique(chidx));
end

%compute spike properties - to save time, just copy them from the previous
%spike sort
spkSort.spkProps=spkSortGui.spkSort.spkProps;
spkSort.info.docu=spkSortGui.spkSort.info.docu;
spkSort.info.docu{10}='Note: spkProps based only on properties of units in partial sort file';
 
%update id - will keep it like this, so that database GUI can deal with it
id.spikeSort.name{probeID}=name;
id.spikeSort.date{probeID}=date;
id.spikeSort.NSingleUnit(probeID)=sum(strcmp(spkSort.unitinfo,'SU'));
id.spikeSort.NMultiUnit(probeID)=sum(strcmp(spkSort.unitinfo,'MU'));
id.spikeSort.Note{probeID}='based on merged sorts';

%save spkSort and id
save(fullfile(expFolder,animalID,expname,[expname  '_p' num2str(probeID) '_spkSort']),'spkSort');
save(fullfile(expFolder,animalID,expname,[expname  '_id']),'id');


if copyToZ==1
    zbase='Z:\EphysNew\processedSpikes';
    save(fullfile(zbase,animalID,expname,[expname  '_p' num2str(probeID) '_spkSort']),'spkSort');
    save(fullfile(zbase,animalID,expname,[expname '_id']),'id');
end





