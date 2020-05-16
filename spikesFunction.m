function spikesFunction(JobID)

% If no isBR in experiment file then set to 0
isBR = 0;
%

load('experiment.mat')

if ~isBR
    
	fileinfo = dir([experiment '/' experiment '_amplifier.dat']);
	samples = fileinfo.bytes/(2*length(CHs));
	samplesPerJob = ceil(samples/200);
	firstSample = samplesPerJob*JobID - 60;

	if firstSample<0
	    firstSample=0;
	end

	File = [experiment '/' experiment '_amplifier.dat'];
	DataFile = fopen(File,'r');
	fseek(DataFile,2*length(CHs)*firstSample,'bof');

	if JobID == 199
	samplesLeft = samples - samplesPerJob*199 + 60;
	Data = fread(DataFile, [length(CHs) samplesLeft], 'int16');
	else
	Data = fread(DataFile, [length(CHs) samplesPerJob], 'int16');
	end

else
	addpath(genpath('NPMK'));

	Data = openNSx([experiment '/' experiment '.ns6'],'t:1:2');
	samples = Data.MetaTags.DataPoints;
	samplesPerJob = ceil(samples/200);
	firstSample = samplesPerJob*JobID - 60;

	if firstSample<1
	    firstSample=1;
	end

	if JobID == 199
	Data = openNSx([experiment '/' experiment '.ns6'],['t:' num2str(firstSample) ':' num2str(samples)]);
	else
	Data = openNSx([experiment '/' experiment '.ns6'],['t:' num2str(firstSample) ':' num2str(firstSample+samplesPerJob-1)]);
	end
	Data = Data.Data;

end

% Filter and normalize to threshold
for i = 1:length(CHs)
Data(i,:) = filter(b1, a1, double(Data(i,:)));
end
Data = Data(CHs,:);
Data = Data(~BadCh(CHs),:);
Data=double(Data);
numChs = sum(~BadCh);

for i = 1:numChs
Data(i,:) = Data(i,:)./abs(Th(i));
end

TimeStamps = [];
Amps = [];
Width = [];
Amplitude = [];
Energy = [];
pos = [];
with = [];
Ch = [];
AboveTh = Data>-1;

% After calculating AboveTh un-normalize data
for i = 1:numChs
Data(i,:) = Data(i,:).*abs(Th(i));
end

minData = [zeros(4,length(Data(1,:)));Data;zeros(4,length(Data(1,:)))];

for i = [1 2 4 8 5]
minData = min(circshift(minData,[0 i]),minData);
end
minData = circshift(minData,[0 -10]);
for i = [1 2 4 1]
minData = min(circshift(minData,[i 0]),minData);
end
minData = circshift(minData,[-4 0]);
minData = minData(5:end-4,:);
CrossTh = AboveTh & circshift(~AboveTh,[0 -1]);
for i = [1 2 4 2]
CrossTh = CrossTh | circshift(CrossTh,[0 i]);
end

% Make sure there is no repeted value of max during artificial refractory
% period (necessary for BR)
RepetedMax = zeros(length(Data(:,1)),length(Data(1,:)));

% Is the first max in 10 smaples across neighbor channels
for ch = -4:4
for sm = -10:-1
RepetedMax = RepetedMax | Data == circshift(Data,[ch sm]);
end
end

% Is the first max at that sample across neighbor channels
for ch = 1:4
RepetedMax = RepetedMax | Data == circshift(Data,[ch 0]);
end

Spikes = CrossTh & minData==Data & ~RepetedMax;
Spikes(:,1:30)=0;
Spikes(:,end-30:end)=0;

for i = 1:numChs
Wvfrms=[];
Times = find(Spikes(i,:)>0);
TimeStamps =[TimeStamps Times+firstSample];
trd=1;
for dt = -15:25
Wvfrms(:,:,trd) = Data(max(i-4,1):min(numChs,i+4),Times+dt);
trd=trd+1;
end
if not(isempty(Wvfrms))
Wvfrms = Wvfrms - repmat(mean(Wvfrms(:,:,1:5),3),[1 1 length(Wvfrms(1,1,:))]);
En = sum(Wvfrms.^2,3);
Mn = squeeze(Wvfrms(:,:,16));
else
Wvfrms=[];
En=[];
Mn=[];
end


% Get max and max pos as first peak after min
Mx=[];
if not(isempty(Wvfrms))
for TheCh = 1:length(Wvfrms(:,1,1))
if ~isempty(Wvfrms)
redWvrm = squeeze(Wvfrms(TheCh,:,16:end));
    if length(Wvfrms(1,:,1))==1
        redWvrm = redWvrm';
    end
dWvrm = redWvrm - circshift(redWvrm,[0 1]);
dWvrm(:,end) = [];
dWvrm(dWvrm==0) = dWvrm(circshift(dWvrm,[0 -1])==0);
dWvrm = dWvrm./abs(dWvrm);
dWvrm = dWvrm - circshift(dWvrm,[0 1]);
dWvrm(:,end) = [];
dWvrm = circshift(dWvrm,[0 -1]);
dWvrm(:,1) = 0;
[tmp MxPos] = min(dWvrm,[],2);
Mx(TheCh,:) = redWvrm((MxPos-1).*length(dWvrm(:,1)) + [1:1:length(dWvrm(:,1))]')';
else
Mx=[];
MxPos = [];
end
end
%

if i<5
Mn = [zeros(5-i,length(Mn(1,:)));Mn];
Mx = [zeros(5-i,length(Mx(1,:)));Mx];
end
if i>(numChs-4)
Mn = [Mn;zeros(i-(numChs-4),length(Mn(1,:)))];
Mx = [Mx;zeros(i-(numChs-4),length(Mx(1,:)))];
end
Pks = Mx-Mn;
Amps =[Amps Pks];
Amplitude = [Amplitude Mn(5,:)];
Width =[Width MxPos'-1];
Energy =[Energy En(5,:)];
Ch(end+1:end+length(Times)) = i;
Mn(Mn>0)=0;
pos = [pos i+(sum(Mn.*(repmat([-4:1:4]',1,length(Mn(1,:)))),1)./sum(Mn,1))];
with = [with (Mn(5,:)./sum(Mn,1))];
end
end

Properties = [Amps;Amplitude;Energy;Width;with;pos;TimeStamps;Ch];
save([experiment '/SpikeFiles/Spikes_' num2str(JobID) '.mat'],'Properties')
fclose all;
