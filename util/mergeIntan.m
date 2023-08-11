function mergeIntan(expFolder,animalID,mergeStruct,parts,outID,name)
% merging of intat files
%
% input parameters:
% expFolder - base folder for experiments (string)
% animalID - animal ID (string)
% mergeStruct - cell array containing list of files to load (as
% {'u000_001','u000_002'})
% parts - number of blocks to divide each file into during reading
% operation
% outID - experiment ID of output file (unit will be set to uMMM per
% default)
% name - name or initials of person running the script (for bookkeeping)



%% mergeInfo structure
mergeInfo=struct;
mergeInfo.animal=animalID;
mergeInfo.outID=outID;
mergeInfo.files=mergeStruct;
mergeInfo.name=name;
mergeInfo.date=date;
mergeInfo.filesize=zeros(1,length(mergeStruct)); %necessary info for splitting files later
mergeInfo.bufferLength=1000;

%% make directory for output
outname=[animalID '_uMMM_' outID];
if ~exist(fullfile(expFolder,animalID,outname),'dir')
    mkdir(fullfile(expFolder,animalID),outname);
end


%% merge headers by making a minimal header for the new file 
%only need sampleFreq and stuff before it
%we could copy one of the headers, but this way it's clear that these are
%not raw intan files
%this also provides the number of channels (needed for buffer)

%read the first header - this only works if the sample freq is the same
%throughout anyways
headerFile=fullfile(expFolder,animalID,[animalID '_' mergeStruct{1}],[animalID '_' mergeStruct{1} '_info.rhd']);
headerIn=read_Intan_Header(headerFile);

%get number of channels
mergeInfo.nChannel=sum(headerIn.signal_group_num_amp_channels(headerIn.signal_group_enabled==1));

%generate output file
hOut=fullfile(fullfile(expFolder,animalID,outname,[outname '_info.rhd']));
fid = fopen(hOut, 'w');

%write magic number
magic_number= 0xC6912702;
fwrite(fid,magic_number,'uint32');

%fake file version
main_version=0;
second_version=0;
fwrite(fid,main_version,'int16');
fwrite(fid,second_version,'int16');

%sample rate
fwrite(fid,headerIn.sample_rate,'single');
fclose(fid);


%% merge amplifier files
%the first file could just be copied, but in the end we need to open it for
%reading anyways - just treat like ther est

%waitbar
totJobs=length(mergeStruct)*parts; 
h=waitbar(0/totJobs,'Merging files...');

%open output for write
outAmp=fullfile(expFolder,animalID,outname,[outname '_amplifier.dat']);

%check that this does not exist already
if exist(outAmp,'file')
    errordlg('Merge file already exists, aborting!','Error');
    return;
end

outFID = fopen(outAmp,'a+'); 

%go through files, read and append; broken into parts to make it
%manageable
for f=1:length(mergeStruct)
    %we need a buffer to mask the transition between the files (filtering
    %creates an artefact otherwise)
    if f>1
        fwrite(outFID,zeros(1,buffer),'int16');
    end
    
    mergeFile=fullfile(expFolder,animalID,[animalID '_' mergeStruct{f}],[animalID '_' mergeStruct{f} '_amplifier.dat']);
    fileinfo = dir(mergeFile);
    mergeInfo.filesize(f)=fileinfo.bytes;
    samplesPerJob = ceil(fileinfo.bytes/parts); % Number of samples in amplifier data file

    mergeFID = fopen(mergeFile,'r'); % Load amplifier data

    for i=0:parts-2 %need to treat the last one differently
        waitbar(((f-1)*parts+i+1)/totJobs,h)
        data = fread(mergeFID,samplesPerJob, 'int16');
        fwrite(outFID,data,'int16');
    end

    %last one - potentially less than full job left
    waitbar(f*parts/totJobs,h)
    samplesLeft=fileinfo.bytes-samplesPerJob*(parts-1);
    data = fread(mergeFID,samplesLeft, 'int16');
    fwrite(outFID,data,'int16');
 
    %use the data for transition to the next one

    fclose(mergeFID);
end

%close output file
fclose(outFID);
close(h);



%% save merge info file
save(fullfile(expFolder,animalID,outname,[outname '_mergeInfo']),'mergeInfo');