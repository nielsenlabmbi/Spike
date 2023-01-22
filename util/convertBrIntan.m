%conversion from blackrock files to intan
%function convertBrIntan(brFile)

brFile='D:\ephys\feai3\u000_000\feai3_u000_000.ns6';

%% read header info
header=openNSx(brFile,'noread');
sampleFreq=header.MetaTags.SamplingFreq;

samples=header.MetaTags.DataPoints;

%% output file names
[outPath,outFile,~]=fileparts(brFile);

%% write header: only need sample rate
%this is written to be compatible with read_intan_header
%note even though this is shorter than the intan header, that function
%works (it returns empty fields for things that are not set)
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

dat=openNSx('uV',filename);
dat.Data;