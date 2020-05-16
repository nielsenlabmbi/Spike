function makeSpikesFile

% Load settings

load([fileparts(mfilename('fullpath')) '\Settings'])
experiment = inputdlg('Enter experiment name: ');
experiment = experiment{1};


% Get smaple rate from info file
sampleFrq = getSampleFreqFromInfoFile([expFolder '/' experiment '/' experiment '_info.rhd']);
%

% Get Events

DigiFile = fopen([expFolder '/' experiment '/' experiment '_digitalin.dat']);
fileinfo = dir([expFolder '/' experiment '/' experiment '_digitalin.dat']);
Digital = fread(DigiFile, (fileinfo.bytes)/2, 'uint16');
dDigital = diff(Digital);
One = find(dDigital == 1);
Two = find(dDigital == 2);
Three = find(dDigital == -2);
if length(Two) ~= length(Three)
Two = find(dDigital == 6);
end
if length(Two) ~= length(One)
Two = find(dDigital == 6 | dDigital == 2);
Three = find(dDigital == -6 | dDigital == -2);
end
Fourth = find(dDigital == -1);

% Get the ones that are not truncated

Events.Timestamp{1}=[];
Events.Timestamp{1}(:,1)=One(1:length(Fourth));
Events.Timestamp{1}(:,2)=Two(1:length(Fourth));
Events.Timestamp{1}(:,3)=Three(1:length(Fourth));
Events.Timestamp{1}(:,4)=Fourth(1:length(Fourth));
%

% Load spikes file
load([expFolder '/' experiment '/SpikeFiles/' experiment '_Spikes.mat']);

% Initialize spikes and set mapping
numCHs = max(Properties(16,:));
Spikes={};
Mapping(1,:) = 1:1:numCHs;
Mapping(2,:) = ones(1,numCHs);

% Add 1 to idk so that first unit number is 1 (is assumed by downstream
% code)
idk = idk+1;

% Populate Spikes structure and make unti position
Post = unique(idk);

%if missing spikes at end, add them to multiunit
if length(idk) ~= length(Properties(1,:))
    idk(end:length(Properties(1,:))) = 0;
end

for Unit = 1:length(Post)
UnitPosition(Unit) = mean(Properties(14,idk==Post(Unit)));
end

UnitType{1} = ones(length(Post),1);

[Spikes{1}.TimeStamp TimeOrder] = sort(Properties(15,:));
Spikes{1}.Unit = idk(TimeOrder);

Spikes{1}.Waveform=zeros(90,length(Mapping(1,:)),length(Post));

% Determine unit type by contamination below 1.2msec
samplesRefPeriod = round((sampleFrq/1000)*1.2);
for Unit = 1:length(Post)
TSt = Spikes{1}.TimeStamp(Spikes{1}.Unit==Post(Unit));
dTSt = diff(TSt);
PostCont = sum(dTSt<(samplesRefPeriod));
PostSpks = length(TSt);
if PostCont == 0
UnitType{1}(Unit) = 1;
else
if PostCont/PostSpks < 0.0005
UnitType{1}(Unit) = 2;
else
UnitType{1}(Unit) = 3;
end
end
end

% Save spikes file
save([spikesFolder '/' experiment '_spikes.mat'],'Spikes','UnitType','Events','Mapping','Properties','UnitPosition','sampleFrq')

% Load Analyzer file
load([analyzerFileFolder '/' experiment(1:5) '/' experiment analyzerFileEnding],'-mat');
% Analyze responses
Params = Analyzer.L.param{1,1};
Reps = length(Analyzer.loops.conds{1,1}.repeats);
BReps =length(Analyzer.loops.conds{1,end}.repeats);
Conds = length(Analyzer.loops.conds);
    
if isequal(Analyzer.loops.conds{1,end}.symbol{1,1}, 'blank')
    Conds = length(Analyzer.loops.conds)-1;
    for r = 1:BReps
        Trial = Analyzer.loops.conds{1,end}.repeats{1,r}.trialno;
        TrialInfo(Trial,:) = [repmat([(zeros(length(Analyzer.loops.conds{1,end-1}.val(:,:)),1)-1).' 0],[length(Trial) 1]) double(Events.Timestamp{1}(Trial,:))];
    end
end
for i = 1:Conds
    for r = 1:Reps
        Trial = Analyzer.loops.conds{1,i}.repeats{1,r}.trialno;
        if not(Trial > length(Fourth))
        TrialInfo(Trial,:) = [vertcat(Analyzer.loops.conds{1,i}.val{:,:}).' i double(Events.Timestamp{1}(Trial,:))];
        end
    end
    CondInfo(i,:)= vertcat(Analyzer.loops.conds{1,i}.val{:,:}).';
end
Parmass = Analyzer.L.param;
VarValDims = 1;
for i=1:length(Parmass(1,:))
    VarValDims = VarValDims*(length(eval(Parmass{1,i}{2})));
end

for site = 1:length(Spikes)
Units = unique(Spikes{site}.Unit);
Units = sort(Units);
if Units(1) == 0
    Units = Units(2:end);
end
Units = length(Units);
Data{site}.Spiking=cell(VarValDims, Units, Reps);
Data{site}.BSpiking=cell(Units, BReps);
Repetition = ones(VarValDims,1);
for T = 1:length(TrialInfo(:,1))
    if not(isequal(TrialInfo(T,1:length(Parmass(1,:))), (zeros(1,length(Parmass(1,:))) -1)))
    for Unit=1:Units
    VariableIndex = TrialInfo(T,length(Parmass(1,:))+1);
    Data{site}.Spiking{VariableIndex,Unit,Repetition(VariableIndex)} = double(Spikes{site}.TimeStamp(find(Spikes{site}.Unit == Unit & Spikes{site}.TimeStamp > TrialInfo(T,length(Parmass(1,:))+2) & Spikes{site}.TimeStamp < TrialInfo(T,length(Parmass(1,:))+5))))-TrialInfo(T,length(Parmass(1,:))+3);
    end
    Repetition(VariableIndex) = Repetition(VariableIndex)+1;
    end
end
Repetition = 1;
for T = 1:length(TrialInfo(:,1))
    if TrialInfo(T,1:length(Parmass(1,:))) == zeros(1,length(Parmass(1,:))) -1
    for Unit=1:Units
    Data{site}.BSpiking{Unit,Repetition} = double(Spikes{site}.TimeStamp(find(Spikes{site}.Unit == Unit & Spikes{site}.TimeStamp > TrialInfo(T,length(Parmass(1,:))+2) & Spikes{site}.TimeStamp < TrialInfo(T,length(Parmass(1,:))+5))))-TrialInfo(T,length(Parmass(1,:))+3);
    end
    Repetition = Repetition+1;
    end
end

Paramss = Analyzer.L.param;
FixIndx = zeros(length(Paramss(1,:)),1)-1;
for Unit = 1:length(UnitType{site})
%Calc Blank Responses
RepVal = [];
    for rep = 1:length(Data{site}.BSpiking(1,:))
        Spks = [Data{site}.BSpiking{Unit,rep}];
        RepVal(rep) =((sum(Spks>0 & Spks < sampleFrq*Analyzer.P.param{1,3}{3})/(Analyzer.P.param{1,3}{3}))-((sum(Spks<0))/Analyzer.P.param{1,1}{3}));
    end
Data{site}.BRespMean(Unit) = mean(RepVal);
Data{site}.BRespVar(Unit) = std(RepVal)/sqrt(length(RepVal));
Data{site}.AllBResp(:,Unit) = RepVal;

%Calc responses for all conditions
for i = 1:length(Data{site}.Spiking(:,1,1))
RepVal = [];
    for rep = 1:length(Data{site}.Spiking(1,1,:))
        Spks = [Data{site}.Spiking{i,Unit,rep}];
        RepVal(rep) =((sum(Spks>0 & Spks < sampleFrq*Analyzer.P.param{1,3}{3})/(Analyzer.P.param{1,3}{3}))-((sum(Spks<0))/Analyzer.P.param{1,1}{3}));
    end
Data{site}.RespVar(Unit,i) = std(RepVal)/sqrt(length(RepVal));
Data{site}.RespMean(Unit,i) = mean(RepVal);
Data{site}.AllResp(Unit,:,i) = RepVal;
end
end
end
RespFunc = [0 sampleFrq*Analyzer.P.param{1,3}{3} 0 -sampleFrq*Analyzer.P.param{1,1}{3}];
Values = cell(length(Paramss(1,:)),1);
for i = 1:length(Paramss(1,:))
    Params = Paramss{1,i};
    Variables{i} =Params{1,1};
    Values{i} = eval(Params{1,2});
end

% Save data file
save([dataFileFolder '/' experiment dataFileEnding],'Data', 'UnitType','Variables','Values','RespFunc','CondInfo','TrialInfo','sampleFrq');

% Write table if desiered
% Ask
choice = questdlg('Do you want to save this record in the summary file?', ...
        'Summary file save', ...
        'Yes','No','Yes');

switch choice
        case 'No'
            return;
end
% Write
[num,txt,raw] = xlsread(summaryFile,1);
Rows=[];
ExperimentName = experiment;

for i = 1:length(raw(:,1))
    if isequal(raw(i,3),{ExperimentName})
        Rows = [Rows i];
    end
end
if isempty(Rows)
    errordlg('Experiment has not been added to table yet')
    return
end
if length(Rows) > 1
    warndlg('Table showing previous units, deleting old units')
    raw(Rows(2:end),:) = [];
    Rows = Rows(1);
end
for site = 1:length(Data)
if site~=1
    raw = [raw(1:Rows(1),:);raw(Rows(1),:);raw(Rows(1)+1:end,:)];
    Rows(1)=Rows(1)+1;
end
raw = [raw(1:Rows(1),:);cell(length(UnitType{site})-1,length(raw(1,:)));raw(Rows(1)+1:end,:)];
raw(Rows(1):Rows(1)+length(UnitType{site})-1,1) = {lower(ExperimentName(1:5))};
raw(Rows(1):Rows(1)+length(UnitType{site})-1,3) = {ExperimentName};
raw(Rows(1):Rows(1)+length(UnitType{site})-1,11) = cell(length(UnitType{site}),1);
raw(Rows(1):Rows(1)+length(UnitType{site})-1,19) = raw(Rows(1),19);
raw(Rows(1):Rows(1)+length(UnitType{site})-1,2) = raw(Rows(1),2);
raw(Rows(1):Rows(1)+length(UnitType{site})-1,4) = raw(Rows(1),4);
raw(Rows(1):Rows(1)+length(UnitType{site})-1,5) = raw(Rows(1),5);
raw(Rows(1):Rows(1)+length(UnitType{site})-1,6) = raw(Rows(1),6);
raw(Rows(1):Rows(1)+length(UnitType{site})-1,7) = {site};
raw(Rows(1):Rows(1)+length(UnitType{site})-1,12:15) = repmat([{3} {date} {3} {date}],length(UnitType{site}),1);
A = min(min(Spikes{site}.Waveform,[],2),[],3);
if isempty(Spikes{site}.Waveform)
MxPos = zeros(1,length(Spikes{site}.TimeStamp));
else
for i = 1:length(A)
MxPos(i) = find(Spikes{site}.Waveform(i,:,:) == A(i),1);
end
MxPos = ceil(MxPos/length(Spikes{site}.Waveform(1,:,1)));
end
if length(MxPos) == length(Spikes{site}.TimeStamp)
for i = 1:length(UnitType{site})
    raw(Rows(1)+i-1,8) = {i};
    raw(Rows(1)+i-1,9) = {UnitType{site}(i)};
    Spks = Spikes{site}.Unit == i;
    TStmps = double(Spikes{site}.TimeStamp(Spks))-MxPos(Spks)';
    ISI = TStmps(2:end) - TStmps(1:end-1);
    Cont = sum(ISI<samplesRefPeriod)/length(ISI);
    raw(Rows(1)+i-1,10) = {Cont};
end
else
for i = 1:length(UnitType{site})
    raw(Rows(1)+i-1,8) = {i};
    raw(Rows(1)+i-1,9) = {UnitType{site}(i)};
    Spks = Spikes{site}.Unit == i;
    TStmps = double(Spikes{site}.TimeStamp(Spks));
    ISI = TStmps(2:end) - TStmps(1:end-1);
    Cont = sum(ISI<samplesRefPeriod)/length(ISI);
    raw(Rows(1)+i-1,10) = {Cont};
end
end
Rows(1)=Rows(1)+length(UnitType{site})-1;
end
xlswrite(summaryFile,raw,1);