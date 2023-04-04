function extractTrialsBR(physpath,anapath,outpath,animal,unit,exp)

%extracts trial information from the analyzer and digital file
%input:
%physpath: path to phys data (e.g., z:\ephysNew\data)
%anapath: path to analyzer data (e.g., z:\ephysNew\analyzer)
%outpath: path to output data (e.g., z:\ephysNew\processedSpikes)
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
nrTrials=getnotrials(Analyzer);

[dom,domval,blankId] = getdomainvalue(Analyzer);
trialInfo.dom=dom;
trialInfo.domval=domval;
trialInfo.blankId=blankId;
trialInfo.triallist=getcondtrial(Analyzer);



%load digital file and extract event onsets
physname=fullfile(physpath,animal,['u' unit '_' exp],[animal '_u' unit '_' exp '.nev']);
NEV = openNEV(physname);
pulseTime = double(NEV.Data.SerialDigitalIO.TimeStamp);
pulseVal=double(NEV.Data.SerialDigitalIO.UnparsedData);
pulseVal=pulseVal';

timeOut=pulseTime;
valOut=pulseVal;

%first issue - get rid of noise spikes (pairs in quick succession)
%time is in samples, so anything shorter than a frame refresh (8ms)
%unreasonable - that's 240 samples at 30kHz; in reality the noise events
%have a difference time of 1 sample
diffTime=diff(pulseTime);
idx=find(diffTime<50);

valOut([idx idx+1])=[];
timeOut([idx idx+1])=[];


%translate events into individual channels
digOut=de2bi(valOut,5);
eventTimes=[];
eventId=[];

%easiest to clean up by just keeping the ones we need and creating new
%array
%column 3 is the trial start, stop marker
diffOne=diff([0;digOut(:,3)]); %otherwise the first one does not show up

tmpTrialStart=find(diffOne==1);
nTrialStart=length(tmpTrialStart);
eventTimes=timeOut(tmpTrialStart);
eventId=ones([1 nTrialStart]);
if nTrialStart~=nrTrials
    disp('mismatch between trial start events and nr trials!')
    disp(['nr trials: ' num2str(nrTrials)])
    disp(['nr start events: ' num2str(nTrialStart)])
end

tmpTrialStop=find(diffOne==-1);
nTrialStop=length(tmpTrialStop);
eventTimes=[eventTimes timeOut(tmpTrialStop)];
eventId=[eventId zeros([1 nTrialStop])];
if nTrialStop~=nrTrials
    disp('mismatch between trial stop events and nr trials!')
    disp(['nr trials: ' num2str(nrTrials)])
    disp(['nr stop events: ' num2str(nTrialStop)])
end

%as a safety check, display trial duration 
disp(['min trial time (s): ' num2str(min(timeOut(tmpTrialStop)-timeOut(tmpTrialStart))/30000)]);
disp(['max trial time (s): ' num2str(max(timeOut(tmpTrialStop)-timeOut(tmpTrialStart))/30000)]);


%column 4 is the stimulus start, stop marker 
diffTwo=diff([0;digOut(:,4)]); %adding an element makes the indexing easier 

tmpStimStart=find(diffTwo==1);
nStimStart=length(tmpStimStart);
eventTimes=[eventTimes timeOut(tmpStimStart)];
eventId=[eventId 3*ones([1 nStimStart])];
if nStimStart~=nrTrials
    disp('mismatch between stim start events and nr trials!')
    disp(['nr trials: ' num2str(nrTrials)])
    disp(['nr start events: ' num2str(nStimStart)])
end

tmpStimStop=find(diffTwo==-1);
nStimStop=length(tmpStimStop);
eventTimes=[eventTimes timeOut(tmpStimStop)];
eventId=[eventId ones([1 nStimStop])];
if nStimStop~=nrTrials
    disp('mismatch between stim stop events and nr trials!')
    disp(['nr trials: ' num2str(nrTrials)])
    disp(['nr stop events: ' num2str(nStimStop)])
end

%as a safety check, display stimulus duration 
disp(['min stim time (s): ' num2str(min(timeOut(tmpStimStop)-timeOut(tmpStimStart))/30000)]);
disp(['max stim time (s): ' num2str(max(timeOut(tmpStimStop)-timeOut(tmpStimStart))/30000)]);


%sort output according to time
[trialInfo.eventTimes,sortidx]=sort(eventTimes);
trialInfo.eventId=eventId(sortidx);
trialInfo.eventCh=de2bi(trialInfo.eventId,3);

%save
outname=fullfile(outpath,animal,[animal '_u' unit '_' exp],[animal '_u' unit '_' exp '_trialInfo.mat']);
save(outname,'trialInfo')

