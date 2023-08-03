function applySortFast(expFolder,animalID,unitID,expID,probeID,name,copyToZ,jobID)
% apply sorting steps
% note: this will not respect any sorting steps that require units to be
% made invisible since it is applied to each jobs file independently
% generates one sort file per jobs file to be merged later (merging will
% handle
%
% input parameters:
% expFolder - base folder for experiments (string)
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% probeID - probe ID (number)
% name - name or initials of person running the script (for bookkeeping)
% copyToZ - copy id file to Z?
% jobID - job ID of raw spike file to process (number)

% output:
% one file
% structure spk with fields

%load history - assumption is this is saved as partSortHist
expname=[animalID '_u' unitID '_' expID];
histname=[expname '_partSortHist'];
load(fullfile(expFolder,animalID,expname,histname)); %generates sortHist

%get data - generates spk
load(fullfile(expFolder,animalID,expname,'SpikeFiles',[expname  '_j' num2str(jobID) '_p' num2str(probeID) '_spike']));

%start building output
spkSort.unitid=zeros(size(spk.detCh));
spkSort.spktimes=spk.spkTimesDet; %times of every spike
spkSort.detCh=spk.detCh; %channel of every spike
spkSort.detChSort=spk.detChSort; %channel as sorted according to probe

nrUnits=0;

%go through and apply the steps per the history
for i=1:length(sortHist)

    %apply view parameters (set everything, independent
    %of selected mode)
    viewProp.xProp=sortHist(i).viewX;
    viewProp.yProp=sortHist(i).viewY;
    viewProp.selectedCh=sortHist(i).selectedCh;

    viewProp.xPropTet=sortHist(i).viewXTet;
    viewProp.yPropTet=sortHist(i).viewYTet;
    viewProp.detChTet=sortHist(i).detChTet;
    viewProp.xChTet=sortHist(i).xChTet;
    viewProp.yChTet=sortHist(i).yChTet;

    viewProp.xPropCh=sortHist(i).viewXCh;
    viewProp.yPropCh=sortHist(i).viewYCh;
    viewProp.detChCh=sortHist(i).detChCh;
    viewProp.detChType=sortHist(i).detChType;

    viewProp.plotmode=sortHist(i).plotmode;

    %set unitid
    unitid=sortHist(i).unitid;

    switch sortHist(i).type
        case 'new'
            roi.Position=sortHist(i).roi{:};
            spkSort.unitid=assignUnit(spkSort.unitid,spkInfo,viewProp,unitid,roi,1);
            nrUnits=nrUnits+1;

        case 'addP'
            roi.Position=app.sortActions(i).roi{:};
            spkSort.unitid=assignUnit(spkSort.unitid,spkInfo,viewProp,unitid,roi,1);

        case 'addM'
            roi.Position=app.sortActions(i).roi{:};
            spkSort.unitid=assignUnit(spkSort.unitid,spkInfo,viewProp,unitid,roi,-1);

        case 'rem'
            spkSort.unitid(spkSort.unitid==unitid)=0;
            %renumber the remainder
            for u=unitid+1:nrUnits
                spkSort.unitid(spkSort.unitid==u)=u-1;
            end
            nrUnits=nrUnits-1;

        case 'merge'
            delUnit=max(unitid);
            keepUnit=min(unitid);
            spkSort.unitid(spkSort.unitid==delUnit)=keepUnit;
            %renumber the remainder
            for u=delUnit+1:nrUnits
                spkSort.unitid(spkSort.unitid==u)=u-1;
            end
            nrUnits=nrUnits-1;

    end
end

%save data
spkSort.info.probeid=ProbeID;
spkSort.info.name=name;
spkSort.info.date=date;
spkSort.info.job=jobID;
spkSort.info.expname=expname;
spkSort.info.historyFile=histname;

save(fullfile(expFolder,animalID,expname,'SpikeFiles',[expname  '_j' num2str(jobID) '_p' num2str(probeID) '_spkSort']),'spkSort');


end

function sortVec=assignUnit(sortVec,spkInfo,viewProp,unitid,roihandle,actionflag)
%edit unitid vector
%actionflag:
% 1: add a positive roi
% -1: add a negative roi
% this operates on the data as shown (i.e.,
% artefacts are not included)

%get data
%generate indices into full spk array
idxData=[1:length(sortVec)];

%remove any artefacts (labeled as -1)
%tf=app.spkInfo.unitid~=-1;

%remove all duplicates
tf = spkInfo.flagDuplicate==0;

%plotmode 1: probe, 2: tetrode, 3: channel
switch viewProp.plotmode
    case 1
        %reduce to selected channels or to correct probe only
        if ~isempty(viewProp.selectedCh) %all data from one probe
            tf = tf & ismember(spkInfo.detCh,viewProp.selectedCh);
        end
      
        %get data
        xData=spkInfo.(viewProp.xProp)(tf);
        yData=spkInfo.(viewProp.yProp)(tf);
        idxData=idxData(tf);

    case 2
        %we only need the data for events that are detected on the
        %detection channels
        tf=tf & (spkInfo.detCh==viewProp.detChTet);

        xData=spkInfo.(viewProp.xPropTet)(tf);
        yData=spkInfo.(viewProp.yPropTet)(tf);
        idxData=idxData(tf);

        if ~isempty(idxData)
            %turn into regular array and reshape
            nCh=size(xData{1},2);
            xData=cell2mat(xData);
            yData=cell2mat(yData);
            xData=reshape(xData,nCh,sum(tf));
            yData=reshape(yData,nCh,sum(tf));

            %now get the data for the channels selected for x and y
            chList=spkInfo.chListAll(tf);
            chList=chList{1};

            xidx=find(chList==viewProp.xChTet);
            yidx=find(chList==viewProp.yChTet);
            xData=xData(xidx,:);
            yData=yData(yidx,:);           
        end
    case 3
        %we only need the data for events that are detected on the
        %detection channel
        if viewProp.detChType==1 %use detChSort
            tf=tf & (spkInfo.detChSort==viewProp.detChCh);
        else
            tf=tf & (spkInfo.detCh==viewProp.detChCh);
        end
       
        xData=spkInfo.(viewProp.xPropCh)(tf);
        yData=spkInfo.(viewProp.yPropCh)(tf);
        idxData=idxData(tf);
end

%get data in ROI
selUnit=inpolygon(xData,yData,roihandle.Position(:,1),roihandle.Position(:,2)); %logical array

%edit unitid: add unit for all but negative rois
if actionflag~=-1
    sortVec(idxData(selUnit))=unitid;
else
    tf=selUnit & sortVec(idxData)==unitid;
    sortVec(idxData(tf))=0;
end

end

