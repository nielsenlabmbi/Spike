function makeExperimentFile

load([fileparts(mfilename('fullpath')) '\Settings'])
experiment = inputdlg('Enter experiment name: ');
experiment = experiment{1};
probeType = inputdlg('Enter probr type (D,F,etc..): ');
probeType = probeType{1};
isBR = 0; % BR=Black Rock
switch probeType
    case 'M'
        CHs = [56,10,55,11,54,12,53,13,52,14,51,15,50,16,49,17,48,18,47,19,46,20,45,21,44,22,43,23,42,24,41,25,40,26,39,27,38,28,37,29,36,30,35,31,34,32,33,9,57,8,58,7,59,6,60,5,61,4,62,3,63,2,64,1];
        BadCh = zeros(1,64);
    case 'D'
        CHs = [17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54];
        BadCh = zeros(1,64);
    case 'F'
        CHs = [13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49] ;
        BadCh = zeros(1,64);
    case 'FF'
        CHs = [[13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49] [13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,64)];
            case '2'
            BadCh = [ones(1,64) zeros(1,64)];
        end
    case 'FD'
        CHs = [[13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49] [17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,64)];
            case '2'
            BadCh = [ones(1,64) zeros(1,64)];
        end

    case 'DF'
        CHs = [[17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54] [13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,64)];
            case '2'
            BadCh = [ones(1,64) zeros(1,64)];
        end
    case 'DD'
        CHs = [[17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54] [17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,64)];
            case '2'
            BadCh = [ones(1,64) zeros(1,64)];
        end
    case 'DM'
        CHs = [[17,47,1,18,46,64,19,45,2,20,44,63,21,43,3,22,42,62,23,41,4,24,40,61,25,39,5,26,38,60,27,37,6,28,36,59,29,35,7,30,34,58,31,33,8,32,48,57,16,49,9,15,50,56,14,51,10,13,52,55,12,53,11,54] [56,10,55,11,54,12,53,13,52,14,51,15,50,16,49,17,48,18,47,19,46,20,45,21,44,22,43,23,42,24,41,25,40,26,39,27,38,28,37,29,36,30,35,31,34,32,33,9,57,8,58,7,59,6,60,5,61,4,62,3,63,2,64,1]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,64)];
            case '2'
            BadCh = [ones(1,64) zeros(1,64)];
        end
    case 'FM'
        CHs = [[13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49] [56,10,55,11,54,12,53,13,52,14,51,15,50,16,49,17,48,18,47,19,46,20,45,21,44,22,43,23,42,24,41,25,40,26,39,27,38,28,37,29,36,30,35,31,34,32,33,9,57,8,58,7,59,6,60,5,61,4,62,3,63,2,64,1]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,64)];
            case '2'
            BadCh = [ones(1,64) zeros(1,64)];
        end
    case 'FT'
        CHs = [[13,20,12,21,11,22,10,23,9,24,8,25,7,26,6,27,5,28,4,29,3,30,2,31,1,32,14,19,15,18,16,17,45,52,44,53,43,54,42,55,41,56,40,57,39,58,38,59,37,60,36,61,35,62,34,63,33,64,46,51,47,50,48,49] [1,2,3,4]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,4)];
            case '2'
            BadCh = [ones(1,64) zeros(1,4)];
        end
    case 'DBR'
        CHs = [31 62 2 29 60 33 27 58 4 25 56 35 23 54 6 21 52 37 19 50 8 17 48 39 15 46 10 13 44 41 11 42 12 9 40 43 7 38 14 5 36 45 3 34 16 1 64 47 32 63 18 30 61 49 28 59 20 26 57 51 24 55 22 53];
        isBR=1;
        BadCh = zeros(1,64);
    case 'FBR'
        CHs = [26:-1:1 28 27 30 29 32 31 58:-1:33 60 59 62 61 64 63];
        isBR=1;
        BadCh = zeros(1,64);
    case 'DD1BR'
        CHs = [[31 62 2 29 60 33 27 58 4 25 56 35 23 54 6 21 52 37 19 50 8 17 48 39 15 46 10 13 44 41 11 42 12 9 40 43 7 38 14 5 36 45 3 34 16 1 64 47 32 63 18 30 61 49 28 59 20 26 57 51 24 55 22 53] [31 2 29 27 4 25 23 6 21 19 8 17 15 10 13 11 12 9 7 14 5 3 16 1 32 18 30 28 20 26 24 22]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,32)];
            case '2'
            BadCh = [ones(1,64) zeros(1,32)];
        end
        isBR=1;
    case 'FF1BR'
        CHs = [[26:-1:1 28 27 30 29 32 31 58:-1:33 60 59 62 61 64 63] [26:-1:1 28 27 30 29 32 31]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,32)];
            case '2'
            BadCh = [ones(1,64) zeros(1,32)];
        end
        isBR=1;
    case 'FF2BR'
        CHs = [[26:-1:1 28 27 30 29 32 31 58:-1:33 60 59 62 61 64 63] [58:-1:33 60 59 62 61 64 63]+32];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,32)];
            case '2'
            BadCh = [ones(1,64) zeros(1,32)];
        end
        isBR=1;
    case 'DD2BR'
        CHs = [[31 62 2 29 60 33 27 58 4 25 56 35 23 54 6 21 52 37 19 50 8 17 48 39 15 46 10 13 44 41 11 42 12 9 40 43 7 38 14 5 36 45 3 34 16 1 64 47 32 63 18 30 61 49 28 59 20 26 57 51 24 55 22 53] [62 60 33 58 56 35 54 52 37 50 48 39 46 44 41 42 40 43 38 36 45 34 64 47 63 61 49 59 57 51 55 53]+32];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,32)];
            case '2'
            BadCh = [ones(1,64) zeros(1,32)];
        end
        isBR=1;
    case 'DF1BR'
        CHs = [[31 62 2 29 60 33 27 58 4 25 56 35 23 54 6 21 52 37 19 50 8 17 48 39 15 46 10 13 44 41 11 42 12 9 40 43 7 38 14 5 36 45 3 34 16 1 64 47 32 63 18 30 61 49 28 59 20 26 57 51 24 55 22 53] [26:-1:1 28 27 30 29 32 31]+64];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,32)];
            case '2'
            BadCh = [ones(1,64) zeros(1,32)];
        end
        isBR=1;
    case 'DF2BR'
        CHs = [[31 62 2 29 60 33 27 58 4 25 56 35 23 54 6 21 52 37 19 50 8 17 48 39 15 46 10 13 44 41 11 42 12 9 40 43 7 38 14 5 36 45 3 34 16 1 64 47 32 63 18 30 61 49 28 59 20 26 57 51 24 55 22 53] [58:-1:33 60 59 62 61 64 63]+32];
        probeToSort = inputdlg('Select probe to sort: ');
        switch probeToSort{1}
            case '1'
            BadCh = [zeros(1,64) ones(1,32)];
            case '2'
            BadCh = [ones(1,64) zeros(1,32)];
        end
        isBR=1;
end

mkdir([expFolder '/' experiment '/'],'SpikeFiles')
% If is not BR get smaple rate from info file and load .dat
if ~(isBR)
sampleFrq = getSampleFreqFromInfoFile([expFolder '/' experiment '/' experiment '_info.rhd']);
File = [expFolder '/' experiment '/' experiment '_amplifier.dat']; % Path to voltage data
DataFile = fopen(File,'r'); % Opens voltage data file for reading
Data = fread(DataFile, [length(BadCh) sampleFrq*30], 'int16'); % reads first 30 seconds of voltage data for all channels

% If it is BR then open NSx
else
sampleFrq = 30000;
disp([expFolder '/' experiment '/' experiment '.ns6'])
Data = openNSx([expFolder '/' experiment '/' experiment '.ns6'],['t:1:' num2str(sampleFrq*30)]); % loads first 30 seconds of voltage data for all channels    
Data = Data.Data;
end
%


[b1, a1] = butter(3, [highPass/sampleFrq,lowPass/sampleFrq]*2, 'bandpass'); % bandpass filters voltage data

dataWindow = figure('pos',[60 450 750 500]);
Th = ones(length(BadCh),1).*(-200); % default voltage threshold

for i = 1:length(CHs)
Done = BadCh(i);
while not(Done)
Data(i,:) = filter(b1, a1, double(Data(i,:)));
figure(dataWindow)
plot(Data(i,1:sampleFrq*2))
set(gca,'YLim',[Th(i)*2 -Th(i)*2])
hold on; 
hold on; plot([0 sampleFrq*2],[Th(i) Th(i)],'r');hold off;
option = questdlg('',['Ch: ' num2str(i)],'Good','Bad','Change Th','Good');
switch option
    case 'Good'
        Done = 1;
    case 'Bad'
        BadCh(i) = 1;
        Done = 1;
    case 'Change Th'
        thInp = inputdlg('Enter new threshold: ');
        Th(i) = eval(thInp{1});
end
Th(i+1)=Th(i);
end
end
Th=Th(~BadCh);
save([expFolder '/' experiment '/experiment.mat'],'experiment','Th','BadCh','CHs','b1','a1','isBR');