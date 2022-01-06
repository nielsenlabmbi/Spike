function extractTrials(physpath,anapath,animal,unit,exp)

%extracts trial information from the analyzer and digital file
%input:
%physpath: path to phys data (e.g., z:\ephysNew)
%anapath: path to analyzer data (e.g., z:\ephysNew\analyzer)
%animal: animal id
%unit: unit id (string)
%exp: exp id (string)
%
%output:
%structure trialInfo with
%dom: parameters changed in looper
%domval: parameter values for each condition
%blankid: condition id for the blanks
%triallist: condition id for each trial
%eventTimes: time stamp for each digital event (in samples)
%eventId: decimal number for event (0 corresponds to end of trial only)
%eventCh: same as eventId, but turned into on/off for 3 digital channels

%load analyzer and get all releveant info
load(fullfile(anapath,animal,[animal '_u' unit '_' exp '.analyzer']),'-mat');

%run helper functions on analyzer
[dom,domval,blankId] = getdomainvalue(Analyzer);
trialInfo.dom=dom;
trialInfo.domval=domval;
trialInfo.blankId=blankId;
trialInfo.triallist=getcondtrial(Analyzer);

%load digital file and extract event onsets
physname=fullfile(physpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp '_digitalin.dat']);
DigiFile = fopen(physname);
fileinfo = dir(physname);
digiData = fread(DigiFile, (fileinfo.bytes)/2, 'uint16');
dDigital = diff(digiData);

%non-zero entries in dDigital are events
trialInfo.eventTimes=find(dDigital~=0)+1;
trialInfo.eventId=digiData(trialInfo.eventTimes);

%in addition to raw event codes,translate into individual channels on/off
if max(trialInfo.eventId)<=7
    trialInfo.eventCh=de2bi(trialInfo.eventId,3);
else
    trialInfo.eventCh=de2bi(trialInfo.eventId,4);
end

%save
save(replace(physname,'digitalin.dat','trialInfo'),'trialInfo')