function mergeJobSort(expFolder,animalID,unitID,expID,probeID,copyToZ,varargin)
% merging of sort applied to each job file
%
% input parameters:
% expFolder - base folder for experiments (string)
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe ID (number)
% name - name or initials of person running the script (for bookkeeping)
% copyToZ - copy spkSort to Z?
% varargin - optional file suffix

% output:
% file with structure spkSort with fields unitid, spktimes, detCh,detChSort

if ~isempty(varargin)
    tSuffix=varargin{1};
else
    tSuffix='';
end

expname=[animalID '_u' unitID '_' expID];
%load id file - for sample frequency
load(fullfile(expFolder,animalID,expname,[expname '_id']));


%find all spkSort files in SpikeFiles folder
sDir='SpikeFiles';
if ~isempty(tSuffix)
    sDir=[sDir '_' tSuffix];
end
sortFiles=dir(fullfile(expFolder,animalID,expname,sDir,[expname  '_j*_p' num2str(probeID) '_spkSort.mat']));

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
partName=[expname '_p' num2str(probeID) '_partSpkSort'];
if ~isempty(tSuffix)
    partName=[partName '_' tSuffix];
end
spkSortGui=load(fullfile(expFolder,animalID,expname,partName)); %generates field spkSort
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

    %correct unit assignments if there is now an ISI violation for a SU
    if ISIv~=0 && strcmp(spkSort.unitinfo{i},'SU')
        spkSort.unitinfo{i}='MU';
    end

    %footprint
    chidx=spkSort.detCh(spkSort.unitid==i);
    spkSort.unitFP(i)=length(unique(chidx));
end


%compute spike properties - to save time, just copy them from the previous
%spike sort
spkSort.spkProps=spkSortGui.spkSort.spkProps;
 
%copy info from partial spkSort
spkSort.infoPartSort=spkSortGui.spkSort.info; %info is used for info from applySortFast


%save spkSort 
outname=[expname  '_p' num2str(probeID) '_spkSort'];
if ~isempty(tSuffix)
    outname=[outname '_' tSuffix];
end
save(fullfile(expFolder,animalID,expname,outname),'spkSort');


if copyToZ==1
    zbase='Z:\EphysNew\processedSpikes';
    save(fullfile(zbase,animalID,expname,outname),'spkSort');
end





