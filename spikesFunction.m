function spikesFunction(experiment, parts, JobID)
% SPIKESFUNCTION computes spike waveforms and waveform properties.
load Settings.mat
load([expFolderLocal '\' experiment '\experiment.mat'])
Properties=[];
partsOverlapSamples = floor(2/(1000/sampleFrq)); % get 2msec overlap between samples
fileinfo = dir([expFolderLocal '\' experiment '\' experiment '_amplifier.dat']);
samples = fileinfo.bytes/(2*length(CHs)); % Number of samples in amplifier data file
samplesPerJob = ceil(samples/parts); % Number of samples to allocate to each of the 200 jobs
firstSample = samplesPerJob*JobID - partsOverlapSamples; % Sets first sample to process; each job has 1st 2msec overlap with previous job and last 2msec overlap with next job

if firstSample<0
    firstSample=0;
end
File = [expFolderLocal '/' experiment '/' experiment '_amplifier.dat'];
DataFile = fopen(File,'r'); % Load amplifier data
fseek(DataFile,2*length(CHs)*firstSample,'bof'); % Offset from beginning of file

if JobID == parts-1 % The last job - first JobID is 0
    samplesLeft = samples - samplesPerJob*(parts-1) + partsOverlapSamples; % samplesLeft=TotalSamples-SamplesDone+Overhang
    Data = fread(DataFile, [length(CHs) samplesLeft], 'int16'); % If JobID is the last job, read all the samples left
else
    Data = fread(DataFile, [length(CHs) samplesPerJob], 'int16'); % If JobID isn't the last job, read samplesPerJob samples past the file position set by fseek
end

% Filter and normalize to threshold
for i = 1:length(CHs) % For every channel
    Data(i,:) = filter(b1, a1, double(Data(i,:))); % Filter data for this job using filter coefficients from experiment.mat
end
Data = Data(CHs,:); % Order data by channel depth
Data = Data(~BadCh(CHs),:); % Only include data from good channels
Data=double(Data); % Convert data from int16 to double

Th=Th(CHs,:);
Th=Th(~BadCh(CHs));

numChs = sum(~BadCh); % Number of good channels

for i = 1:numChs
    Data(i,:) = Data(i,:)./abs(Th(i)); % Normalize each channel to threshold from experiment.mat to determine which data points are above and below threshold
end
TimeStamps = [];
Amps = [];
Width = [];
NegPeak = [];
Energy = [];
pos = [];
with = [];
Ch = [];
x_cm = [];
y_cm = [];
x_detected = [];
y_detected = [];
shank=[];
% make sure things that refer to channels only include good channels
yPosition=yPosition(~BadCh(CHs));
xPosition=xPosition(~BadCh(CHs));
shankNum=shankNum(~BadCh(CHs));
AboveTh = Data>-1; % If Data<-1 after normalization, voltage doesn't pass threshold to be counted as a spike

% After calculating AboveTh un-normalize data
for i = 1:numChs
    Data(i,:) = Data(i,:).*abs(Th(i));
end

%% implement artificial threshold across .8 msec or a maximum of 30 samples of time (after spike detection) and 4 channels of space up-left and down-right
% this is implemented by detecting minimum peaks of voltage across time and
% scape (channels). For a spike to be valid it has to by the minimum value
% of voltage acros 21 samples (including itself) and +/- 4 channels (as
% long as those channels exist). This condition is in addition to this
% minimum value accourring after th crossing (see below).

minData = [zeros(4,length(Data(1,:)));Data;zeros(4,length(Data(1,:)))]; % Initialize minData with buffers of zeros above and below Data so top and bottom channels aren't compared to each other
circs = 0;
for i = [1 2 4 8 16 30] % Circshift sizes to use: Shift 1 away from original, then 3, 7, 15, 31, 61 as long as comparison window dows not exceed refrTime (half refrTime before and after)
    if i*(1000/sampleFrq) < refrTime
        minData = min(circshift(minData,[0 i]),minData); 
        prevTimeDif = i*(1000/sampleFrq);
        circs = circs + i;
    else
        remTime = refrTime-prevTimeDif;
        CircSamples = floor(remTime/((1000/sampleFrq)));
        if rem(CircSamples,2) ~= 0
            CircSamples=CircSamples-1;
        end
        minData = min(circshift(minData,[0 CircSamples]),minData); 
        circs = circs + CircSamples;
    end
end

minData = circshift(minData,[0 -((circs-1)/2)]); % Centers the timepoint window on the middle of the comparison window

%% Apply artificial threshold across channels. Faster but does not consider distance

% for i = [1 2 4 1] % Shift 1, 3, 7, 8
%     minData = min(circshift(minData,[i 0]),minData); % Find minimum voltage in a 9 channel window
% end
% 
% minData = circshift(minData,[-4 0]); % Centers the 9 channel window on the channel where the minimum originally occurred

%% Mod to do +/- 8chs but only if ch distance is less than refrSpace

chPos(:,1,1) = [zeros(4,1);xPosition;zeros(4,1)];
chPos(:,1,2) = [zeros(4,1);yPosition;zeros(4,1)];

preMinData = minData;

for i = -8:8
    chDistance = sqrt(sum((chPos-circshift(chPos,[i 0 0])).^2,3));
    chDistance = chDistance > refrSpace;
    chDistance = double(chDistance);
    chDistance(chDistance==1)=100;
    minData = min(circshift(preMinData,[i 0])+repmat(chDistance,[1 length(preMinData(1,:))]),minData); % Find minimum voltage in a 16 channel window but ignoring channels more than 100um away
end

minData = minData(5:end-4,:); % Removes the buffer of zeros

%% Detect threshold crossings within .3msec

CrossTh = AboveTh & circshift(~AboveTh,[0 -1]); % Matrix of threshold crossings, same size as Data
for i = [1 2 4 8] % Shift 1, 3, 7, 15
    if i*(1000/sampleFrq) < thrTime
        CrossTh = CrossTh | circshift(CrossTh,[0 i]); % Find all samples that are threshold crossings or within thrTime after a threshold crossing
        prevTimeDif = i*(1000/sampleFrq);
    else
        remTime = thrTime-prevTimeDif;
        CircSamples = floor(remTime/((1000/sampleFrq)));
        if rem(CircSamples,2) ~= 0
            CircSamples=CircSamples-1;
        end
        CrossTh = CrossTh | circshift(CrossTh,[0 CircSamples]); % Find all samples that are threshold crossings or within thrTime after a threshold crossing
    end
end

%% Make sure there is no repeted value of max during artificial refractory period (necessary for raw data with low bit depth)

RepetedMax = zeros(length(Data(:,1)),length(Data(1,:)));
% Is the first max in 15 samples across neighbor channels
for ch = -8:8 % channel window
    for sm = -30:0 % sample window
        if ~(sm==0 && ch==0)
        RepetedMax = RepetedMax | Data == circshift(Data,[ch sm]);
        end
    end
end

%% Final spikes detection

Spikes = CrossTh & minData==Data & ~RepetedMax; % Boolean array; spikes are detected when they are within .3ms of a threshold crossing and are a minimum value within .6msec and 100um and are the first instance of that minmum value.
Spikes(:,1:floor(partsOverlapSamples/2))=0; % Removes spikes detected (with minimum) in the first 1msec overlap at the beginning of each job. This is important as some of these may go beyond recording to get waveform.
Spikes(:,end-floor(partsOverlapSamples/2):end)=0; % Removes from the last msec of overlap at the end as well. This is important as some of these may go beyond recording to get waveform.

%% Find the channel(s) where the shank number changes
shankShift=circshift(shankNum,1);
shankStart=find(shankNum(1:end)~=shankShift(1:end)); % Index of first channel after each shank change
shankStart(end+1)=numChs+1;
shankStart = [1;shankStart];

%%
for i = 1:numChs % For every channel
    shankInd=find(shankStart<=i,1,'last'); % Find which shank channel i is on
    Wvfrms=[]; % Initializes waveform
    Times = find(Spikes(i,:)>0); % Find coordinates where spikes occurred
    TimeStamps =[TimeStamps Times+firstSample]; % Get timestamps of these spikes
    firstId = -floor((wvfTime/2)/(1000/sampleFrq));
    lastId = floor(wvfTime/(1000/sampleFrq));
    peakSample = -firstId+1;
    for dt = firstId:lastId % At timeshifts of dt from spike timestamps get waveform that is minus half of wvfTime from spike to wvfTime after spike
        Wvfrms(:,:,dt-firstId+1) = Data(max(i-4,shankStart(shankInd)):min(shankStart(shankInd+1)-1,i+4),Times+dt);
    end
    if not(isempty(Wvfrms))
        Wvfrms = Wvfrms - repmat(mean(Wvfrms(:,:,1:floor(-firstId/2)),3),[1 1 length(Wvfrms(1,1,:))]); % Normalizes waveform by subtracting baseline; takes mean across first half of time window of waveform before peak, subtracts from every waveform
        En = sum(Wvfrms.^2,3); % Calculates waveform energy for the waveform on each channel of each spike; spikes are dimension 2, channels dimension 1, times dimension 3
        Mn = squeeze(Wvfrms(:,:,peakSample)); % baseline subtracted negative peak of spike
    else
        Wvfrms=[];
        En=[];
        Mn=[];
    end

    % Get max and max pos as first peak of waveform after min
    Mx=[];
    if ~isempty(Wvfrms)
        for TheCh = 1:length(Wvfrms(:,1,1)) % For each of the channels spanned by Wvfrms
            if ~isempty(Wvfrms)
                redWvrm = squeeze(Wvfrms(TheCh,:,peakSample:end)); % Only need iterations of Wvfrm after the spike's negative peak; also squeeze to 2 dimensions, #samples by #timeshifts dt
                if length(Wvfrms(1,:,1))==1 % If there's only one sample, squeeze will make the third dimension of Wvfrms the first instead of the second dimension of redWvrm
                    redWvrm = redWvrm'; % Correct redWvrm's dimensions
                end
                dWvrm = redWvrm - circshift(redWvrm,[0 1]); % From each element of redWvrm, subtract the element before it in its row
                dWvrm(:,end) = []; % The end value is the difference between samples that aren't adjacent in time - remove it
                dWvrm(dWvrm==0) = dWvrm(circshift(dWvrm,[0 -1])==0); % For any sample with no difference between it (i) and the previous (i-1), substitute the difference between (i-1) and and (i-2)
                dWvrm = dWvrm./abs(dWvrm); % dWvrm only contains sign of change, not magnitude
                dWvrm = dWvrm - circshift(dWvrm,[0 1]); % dWvrm is change in sign of change
                dWvrm(:,end) = [];
                dWvrm = circshift(dWvrm,[0 -1]);
                dWvrm(:,1) = 0;
                [~, MxPos] = min(dWvrm,[],2); % MxPos contains times of spike positive peaks
                Mx(TheCh,:) = redWvrm((MxPos-1).*length(dWvrm(:,1)) + [1:1:length(dWvrm(:,1))]')'; % Mx contains the values at those peaks
            else
                Mx=[];
                MxPos = [];
            end
        end

        Mn = [zeros(max(0,shankStart(shankInd)+4-i),length(Mn(1,:)));Mn; zeros(max(0,i-(shankStart(shankInd+1)-5)),length(Mn(1,:)))];
        Mx = [zeros(max(0,shankStart(shankInd)+4-i),length(Mx(1,:)));Mx; zeros(max(0,i-(shankStart(shankInd+1)-5)),length(Mx(1,:)))];

        Pks = Mx-Mn; % Spike amplitude peak to peak
        Amps =[Amps Pks]; % Amp(-4 to +4)
        NegPeak = [NegPeak Mn(5,:)]; % Negative peak at the channel of detection
        Width =[Width MxPos'-1]; % Waveform width
        Energy =[Energy En(5,:)]; % Waveform energy at the channel of detection
        Ch(end+1:end+length(Times)) = i; % CH detected
        Mn(Mn>0)=0;
        pos = [pos i+(sum(Mn.*(repmat([-4:1:4]',1,length(Mn(1,:)))),1)./sum(Mn,1))]; % CH Pos, center of mass
        with = [with (Mn(5,:)./sum(Mn,1))]; % CHs width
        
        x_surround=xPosition(max(i-4,shankStart(shankInd)):min(shankStart(shankInd+1)-1,i+4)); % x positions of the surrounding channels
        y_surround=yPosition(max(i-4,shankStart(shankInd)):min(shankStart(shankInd+1)-1,i+4)); % y positions of the surrounding channels

        x_detected(end+1:end+length(Times))=xPosition(i);
        y_detected(end+1:end+length(Times))=yPosition(i);

        x_cm = [x_cm sum(En.*repmat(x_surround,1,size(En,2)),1)./sum(En,1)]; % calculate x center of mass
        y_cm = [y_cm sum(En.*repmat(y_surround,1,size(En,2)),1)./sum(En,1)]; % calculate y center of mass

        shank(end+1:end+length(Times))=shankNum(i);

    end
end

Properties = [Amps;NegPeak;Energy;Width;with;pos;TimeStamps;Ch;x_detected;y_detected;x_cm;y_cm;shank];
PropTitles{1} = 'Amp (-4)';
PropTitles{2} = 'Amp (-3)';
PropTitles{3} = 'Amp (-2)';
PropTitles{4} = 'Amp (-1)';
PropTitles{5} = 'Amp (0)';
PropTitles{6} = 'Amp (1)';
PropTitles{7} = 'Amp (2)';
PropTitles{8} = 'Amp (3)';
PropTitles{9} = 'Amp (4)';
PropTitles{10} = 'Pk2Pk Amp';
PropTitles{11} = 'Energy';
PropTitles{12} = 'Wvf width';
PropTitles{13} = 'CHs width';
PropTitles{14} = 'CHs pos';
PropTitles{15} = 'Time(samples)';
PropTitles{16} = 'Ch detected';
PropTitles{17} = 'X_Detected';
PropTitles{18} = 'Y_Detected';
PropTitles{19} = 'X_Pos';
PropTitles{20} = 'Y_Pos';
PropTitles{21} = 'Shank';
if ~exist([expFolderLocal '\' experiment '\SpikeFiles'], 'dir')
mkdir([expFolderLocal '\' experiment '\SpikeFiles'])
end
save([expFolderLocal '\' experiment '\SpikeFiles\Spikes_' num2str(JobID) '.mat'],'Properties','parts','PropTitles')
fclose all;