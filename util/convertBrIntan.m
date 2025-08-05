%conversion from blackrock files to intan
%goal: generate a header fil (_info.rhd) containing only the sample rate
%and an amplifier file (_amplifier.dat) containing the raw data 
%this will allow the rest of the pipeline to function normally
%
%input: brFile - name of blackrock file, with .ns6 extension

function convertBrIntan(brFile)

% output file names
[outPath,outFile,~]=fileparts(brFile);


%% write header: only need sample rate
%this is written to be compatible with read_intan_header
%note even though this is shorter than the intan header, that function
%works (it returns empty fields for things that are not set)
%but the things before the sampleFreq need to exist

% read header info to get sample rate
header=openNSx(brFile,'noread');
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
%intan format: int16, to convert to voltage in microvolts multiply by 0.195
dat=openNSx('uV',brFile);

%paused? report back
if dat.RawData.PausedFile==1
    disp('File paused!')
end
disp(['Nr channels: ' num2str(dat.MetaTags.ChannelCount)])

hOut=fullfile(outPath,[outFile '_amplifier.dat']);
fid = fopen(hOut, 'w');
fwrite(fid,dat.Data./0.195,'int16');
fclose(fid);

