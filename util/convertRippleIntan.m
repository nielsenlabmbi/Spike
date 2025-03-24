%conversion from ripple files to intan
%goal: generate a header fil (_info.rhd) containing only the sample rate
%and an amplifier file (_amplifier.dat) containing the raw data 
%this will allow the rest of the pipeline to function normally
%Ripple uses the same nev/nsX file format as blackrock, so this uses the
%NPML library used to read blackrock files 
%rippleFile: name of ns5 file, with extension
%createBin: create binary file? if 0, then only write header 
%Nblocks: number of blocks to divide read/write operation into

function convertRippleIntan(rippleFile,createBinary,Nblocks)

% output file names
[outPath,outFile,~]=fileparts(rippleFile);


%% write header: only need sample rate
%this is written to be compatible with read_intan_header
%note even though this is shorter than the intan header, that function
%works (it returns empty fields for things that are not set)
%but the things before the sampleFreq need to exist

% read header info to get sample rate
header=openNSx(rippleFile,'noread');
sampleFreq=header.MetaTags.SamplingFreq;

hOut=fullfile(outPath,[outFile '_info.rhd']);
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
fwrite(fid,sampleFreq,'single');
fclose(fid);

%% write binary amplifier file
%intan format: int16 (intan scaling: samples to microvolt is x 0.195)

if createBinary==1

    %loop through ripple file
    h=waitbar(0,'Converting file');

    %full length of file:
    NSamples=header.MetaTags.DataPoints;
    NSamplesBlock=ceil(NSamples/Nblocks);

    hOut=fullfile(outPath,[outFile '_amplifier.dat']);
    fid = fopen(hOut, 'w');

    for i=1:Nblocks
        waitbar(i/Nblocks,h);
        %read chunk of data
        startRead=(i-1)*NSamplesBlock+1;
        stopRead=i*NSamplesBlock;
        if stopRead>NSamples
            stopRead=NSamples;
        end
        dat=openNSx(['t:' num2str(startRead) ':' num2str(stopRead)],'sample',rippleFile);
        %write
        fwrite(fid,dat.Data,'int16');
    end

    fclose(fid);

    %paused? report back
    if dat.RawData.PausedFile==1
        disp('File paused!')
    end
    close(h);
end

