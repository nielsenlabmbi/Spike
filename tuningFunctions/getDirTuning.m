function dirData=getDirTuning(FRate,dirdom,plotCurve)

%compute orientation tuning parameters for a unit
%input:
%FRate - firing rate matrix, rep x cond (ori conditions only)
%dirdom - list with orientations per condition 
%plotCurve - 1: plot curves
%
%output - structure oriData with fields
%prefDir - preferred dir
%nullDir - null dir
%prefResp - preferred response 
%nullResp - null response 
%DSI - (pref-null)/pref
%Ldir - length vector (1-DirCircVar)
%Adir - angle of vector
%BW - bandwidth
%BWsmooth - smooth bandwidth, but after smoothing tuning curve first


dirData=struct;

%% average rate and other shortcuts
avgRate=mean(FRate,1,'omitnan');
if size(dirdom,1)>1
    dirdom=dirdom';
end

%% preferred direction - accounts for the fact that the maximum rate might
%%occur multiple times
[prefResp,prefCond]=max(avgRate);

if length(find(avgRate==prefResp))==1
    %only one max - can use it directly
    dirData.prefDir=dirdom(prefCond);
else
    %use convolution to pick the maximum that is surrounded by
    %higher values; using square convolution   
    tc=[avgRate(end) avgRate avgRate(1)];  %wrap around tuning curve first (need 1 on either side)
    idx=find(avgRate==prefResp)+1; %indices for peaks
    rateConv=conv(tc,[1 1 1],'same'); %convolution
    mxConv=rateConv(idx); %these are the convolution values for the peaks, pick the bigger one    
    mx=idx(mxConv==max(mxConv)); 
    mx=mx(1);%if there are multiple, we'll just take the first
    dirData.prefDir=dirdom(mx-1);
end

%% preferred response
dirData.prefResp=avgRate(dirdom==dirData.prefDir);

%% null direction and response
dirData.nullDir=mod(dirData.prefDir+180,360);
dirData.nullResp=avgRate(dirdom==dirData.nullDir);

%% OSI
dirData.DSI=(dirData.prefResp-dirData.nullResp)/dirData.prefResp;

%% vector length - better with rectified response
avgRateR=avgRate;
avgRateR(avgRateR<0)=0;
dirVec=sum(avgRateR.*exp(1i*dirdom*pi/180)); 
dirData.Ldir=abs(dirVec)/sum(avgRateR);
dirData.Adir=mod(angle(dirVec)*180/pi/2,180);

%% bandwidth - as in ringach 2002, with and without smoothing
wHan=hanning(3);
wHan=wHan/sum(wHan);
avgRateS=conv(avgRate,wHan,'same');


%criterion response: peak resp/sqrt(2)
critResp=dirData.prefResp/sqrt(2);
critRespSmooth=max(avgRateS)/sqrt(2);

%wrap around to avoid edge arfecats (need to do that before interpolation
tcWrap=[avgRate avgRate avgRate];
tcWrapS=[avgRateS avgRateS avgRateS];

%interpolate tuning function and orientation vector
nInterp=3000; %3*1000 because of wrap around
domInter=linspace(dirdom(1)-360,dirdom(end)+360,nInterp); %equivalent to 3x range

tcIP=griddedInterpolant(tcWrap); %generates object to be evaluated
tcInter=tcIP(linspace(1,length(tcWrap),nInterp));

tcIPS=griddedInterpolant(tcWrapS); %smooth version
tcInterS=tcIPS(linspace(1,length(tcWrapS),nInterp));

%find crossings with criterion response
[~,~,cIdx]=zerocrossrate(tcInter,'Level',critResp);
cIdx(1)=0;
if isempty(cIdx) || sum(cIdx)<4 %less than 4: only 1 data point of the tc falls below the level
    dirData.BW=NaN;
else
    %find the ones closest to the peak
    crossOri=domInter(cIdx);
    cross1=max(crossOri(crossOri<dirData.prefDir));
    cross2=min(crossOri(crossOri>dirData.prefDir));
    dirData.BW=(cross2-cross1)/2;
end

[~,~,cIdxS]=zerocrossrate(tcInterS,'Level',critRespSmooth);
cIdxS(1)=0;
if isempty(cIdxS) || sum(cIdxS)<4
    dirData.BWSmooth=NaN;
else
    %find the ones closest to the peak
    crossOri=domInter(cIdxS);
    cross1S=max(crossOri(crossOri<dirData.prefDir));
    cross2S=min(crossOri(crossOri>dirData.prefDir));
    dirData.BWSmooth=(cross2S-cross1S)/2;
end



%% plot (if selected)
if plotCurve==1
    figure
    subplot(1,2,1)
    plot([dirdom dirdom(1)+360],[avgRate avgRate(1)],'o-')
    hold on
    plot([cross1 cross2],[critResp critResp],'o--')
    plot([dirData.prefDir dirData.prefDir],[0 dirData.prefResp],'g--')
    
    xlabel('Direction (deg)')
    ylabel('Avg firing rate (Hz)')
    title({'Raw tuning curve';['DSI: ' num2str(dirData.DSI) ...
        ' Ldir: ' num2str(dirData.Ldir) ...
        ' BW: ' num2str(dirData.BW)]});

    subplot(1,2,2)
    plot([dirdom dirdom(1)+360],[avgRateS avgRateS(1)],'o-')
    hold on
    plot([cross1S cross2S],[critRespSmooth critRespSmooth],'o--')
    plot([dirData.prefDir dirData.prefDir],[0 dirData.prefResp],'g--')
    
    xlabel('Direction (deg)')
    ylabel('Avg firing rate (Hz)')
    title({'Smoothed tuning curve';['BWsmooth: ' num2str(dirData.BWSmooth)]})

end
