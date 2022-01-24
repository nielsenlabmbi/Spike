function [psth,binVec]=computePSTH(data,baseTime,stimTime,binWidth)
%compute psth
%input:
%data: data for all trials for 1 unit/channel
%format: cell array as provided for the unit or channel by SUTrialData or
%MUThreshTrialData (1 entry per trial, time of an event is set to 0)
%baseTime: duration of baseline interval before stimulus onset (in s)
%stimTime: duration of stimulus interval after stimulus onset (in s)
%binWidth: width of PSTH bin in ms
%
%output:
%psth: psth for the unit/channel, in spikes/s
%binVec: vector with bin edges
%you can plot the results using
%histogram('BinEdges',binVec,'BinCounts',psth)


%construct bin vector; goal: 0 should be one of the bin edges
startBin=ceil(-baseTime*1000/binWidth)*binWidth; %need multiple of binWidth to make 0 an edge
stopBin=floor(stimTime*1000/binWidth)*binWidth;
binVec=[startBin:binWidth:stopBin];

%compute histogram
N=zeros(length(data),length(binVec)-1);
for i=1:length(data)
    N(i,:)=histcounts(data{i},binVec);
end
N=N./(binWidth/1000);

psth=sum(N,1);