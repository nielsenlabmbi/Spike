function convertOldSpikeFiles(FilePathIn,FileIn,FilePathOut)

%converts spikes or Spikes files generated with previous sorting pipeline
%handles conversion per directory (will convert multiple probe files if
%present)
%FilePathIn: z:\ephys_new\data\fexx0
%FileIn: fexx0_u000_000
%FilePathOut: as FilePathIn
%this will not generate spkSort files for files that only contain MU data

%initialize
minISI=0.0012;


%id file (for sample rate)
idname=fullfile(FilePathIn,FileIn,[FileIn '_id.mat']);

%spike and Spike file names
spikefile1=fullfile(FilePathIn,FileIn,[FileIn '*_spikes.mat']);
spikefile2=fullfile(FilePathIn,FileIn,'SpkFile',[FileIn '*_Spikes.mat']);

%check what is available
spikeList1=dir(spikefile1);
spikeList2=dir(spikefile2);

%construct one file list to work through
%the Spikes file has less potential issues than the spikes file, so take
%the first in cases where we have both
isDup=ismember(lower({spikeList1.name}),lower({spikeList2.name})); %these are the ones for which there is a spikes file
fileList={spikeList2.name}; %all the spikes files
fileType=ones(size(fileList)); %1=Spikes, 2=spikes

%add all of the Spikes files that don't exist as spikes files
for i=1:length(isDup)
    if isDup(i)==0
        fileList{end+1}=spikeList1(i).name;
        fileType(end+1)=2;
    end
end


for i=1:length(fileList)
    spkSort=struct;
    
    if fileType(i)==1 %Spikes File
        spikesIn=load(fullfile(FilePathIn,FileIn,'SpkFile',fileList{i}));
        load(idname);
        
        %figure out number of units
        uidx=unique(spikesIn.idk);
        uidx=uidx(uidx~=0); %0 are the unsorted spikes
        
        %we get unit id from idk 
        if ~isempty(uidx) %standard case
            spkSort.unitid=spikesIn.idk;
        else
            if isfield(spikesIn,'Spikes') && isfield(spikesIn.Spikes{1},'Unit')
                spkSort.unitid=spikesIn.Spikes{1}.Unit;
                uidx=unique(spkSort.unitid);
                uidx=uidx(uidx~=0); %0 are the unsorted spikes
            else
                %just add 1 to all of the units
                spkSort.unitid=spikesIn.idk+1;
                uidx=unique(spkSort.unitid);
                uidx=uidx(uidx~=0); %0 are the unsorted spikes
            end
        end
        
        %spike times are in properties
        spkSort.spktimes=spikesIn.Properties(15,:);
        spkSort.info.convNotes{1}='Time stamps from properties';
        
        %channel is also in properties
        spkSort.detChSort=spikesIn.Properties(16,:);
        spkSort.info.convNotes{2}='channel IDs based on old channel sort';
        
        %unit info - ignoring first (unsorted) unit
        minISIsample=round(minISI*id.sampleFreq);

        spkSort.unitinfo=cell(length(uidx),1);
        for u=1:length(uidx)
            %replicating Augusto's usual classification
            %based on ISI
            ts=spkSort.spktimes(spkSort.unitid==uidx(u)); 
            diffTs=diff(sort(ts));
            
            perISIv=sum(diffTs<minISIsample);
            if perISIv==0
                spkSort.unitinfo{u}='SU';
            else
                spkSort.unitinfo{u}='MU';
            end
            spkSort.unitisi(u)=perISIv/length(ts)*100;
            
            %footprint
            chidx=spkSort.detChSort(spkSort.unitid==u);
            spkSort.unitFP(u)=length(unique(chidx));
            
        end
        
        %determine whether this should be saved at all - will drop files
        %that only have MU
        saveFlag=1;
        if all(strcmp(spkSort.unitinfo,'MU'))
            saveFlag=0;
        end
            
        
        %create output name
        outname=replace(fileList{i},'Spikes','spkSort');
        
        
    else %spikes file
        spikesIn=load(fullfile(FilePathIn,FileIn,fileList{i}));
        load(idname);
        
        %spike times
        if isfield(spikesIn.Spikes{1},'TimeStampCorr') && length(spikesIn.Spikes{1}.TimeStampCorr)==length(spikesIn.Spikes{1}.Unit)
            spkSort.spktimes=spikesIn.Spikes{1}.TimeStampCorr;
            spkSort.info.convNotes{1}='TimeStampCorr from Spikes';
        else
            spkSort.spktimes=spikesIn.Spikes{1}.TimeStamp;
            spkSort.info.convNotes{1}='TimeStamp from Spikes';
        end
        
        %unit ids
        uidx=unique(spikesIn.Spikes{1}.Unit); %unit ids
        
        %there is a bug in the makeSpikeFiles that results in one extra
        %unit that should be part of the unsorted unit, set to NaN here and
        %remove later (NaN to work independent of unit numbers in the
        %original file)
        if length(uidx)>length(spikesIn.UnitType{1}) && min(uidx)==0
            spikesIn.Spikes{1}.Unit(spikesIn.Spikes{1}.Unit==0)=NaN;
        end
        

        %unit info - if first unit is mu there are 2 cases - either
        %first unit is unsorted spikes (then remove), or
        %first unit is one of many multi units (then keep)
        if spikesIn.UnitType{1}(1)>1 && any(spikesIn.UnitType{1}==1) %MU followed by SU

            %drop first unit as unsorted by setting it to 0
            if min(spikesIn.Spikes{1}.Unit)==1 
                spkSort.unitid=spikesIn.Spikes{1}.Unit-1;
            else
                spkSort.unitid=spikesIn.Spikes{1}.Unit;
            end
            
            Nunit=length(spikesIn.UnitType{1})-1;
            Ustart=2;
        else %first unit is SU, or first unit is MU with only MU after (we will drop that case later)
            spkSort.unitid=spikesIn.Spikes{1}.Unit;
            
            Nunit=length(spikesIn.UnitType{1});
            Ustart=1;
        end
        
        %set any NaNs to 0 in unitid
        spkSort.unitid(isnan(spkSort.unitid))=0;
        
        %sanity check
        if length(spkSort.unitid)~=length(spkSort.spktimes)
            disp('Mismatch between unitid and spiketime vectors')
        end
        
        spkSort.unitinfo=cell(Nunit,1);
        for u=Ustart:length(spikesIn.UnitType{1})
            %copy assignment
            if spikesIn.UnitType{1}(u)==1
                spkSort.unitinfo{u-Ustart+1}='SU';
            else
                spkSort.unitinfo{u-Ustart+1}='MU';
            end
        end
        
        %get channel if possible
        %can't reconstruct original channel because
        %experiment.mat gets overwritten when sorting more than
        %1 probe
        if isfield(spikesIn,'Properties')
            %need to organize to match stuff in Spikes structure
            [~,TimeOrder] = sort(spikesIn.Properties(15,:));
            spkSort.detChSort=spikesIn.Properties(16,TimeOrder);
            spkSort.info.convNotes{2}='channel IDs based on old channel sort';
              
            for u=1:length(spkSort.unitinfo)
                chidx=spkSort.detChSort(spkSort.unitid==u);
                spkSort.unitFP(u)=length(unique(chidx));
            end
        else
            spkSort.detChSort=nan;
        end
        
        %compute ISI - base this on what is copied into
        %spktimes to keep things consistent
        minISIsample=round(minISI*id.sampleFreq);
        for u=1:length(spkSort.unitinfo)
            ts=sort(spkSort.spktimes(spkSort.unitid==u));
            diffTs=diff(ts);
            perISIv=sum(diffTs<minISIsample)/length(ts);
            spkSort.unitisi(u)=perISIv*100;
        end
        
        %determine whether this should be saved at all
        saveFlag=1;
        if all(strcmp(spkSort.unitinfo,'MU'))
            saveFlag=0;
        end
        
        
        %create output name
        outname=replace(fileList{i},'spikes','spkSort');
        
    end %if fileType
    
    if saveFlag==1
        
        %need to figure out probe and add to file name if necessary
        idx=strfind(outname,'_P');
        if isempty(idx)
            probeId=1;
            outname=replace(outname,'spkSort','P1_spkSort');
        else
            probeId=str2num(outname(idx+2));
        end
        
        
        %add experiment name
        idx=strfind(outname,'_P');
        spkSort.info.expname=outname(1:idx-1);
        
        %add other general info
        spkSort.info.name='AL';
        spkSort.info.probeid=probeId;
        spkSort.info.convertedFile=1;
        
        spkSort.info.docu={'unitid - unit id per spike';...
            'spktimes - time of each spike (in samples)';...
            'detChSort - detection channel, renumbered from top to bottom per shank';...
            'unitinfo - category for each unit';...
            'unitisi - ISI for each unit';...
            'unitFP - unit footprint (detected on how many channels)'};
        
        %save
        outfile=fullfile(FilePathOut,FileIn,outname);
        
        save(outfile,'spkSort');
    end
    
end






