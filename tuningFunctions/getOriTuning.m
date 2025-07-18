function oriData=getOriTuning(FRate,oridom,plotCurve)

%compute orientation tuning parameters for a unit
%input:
%FRate - firing rate matrix, rep x cond (ori conditions only, collapsed across directions)
%oridom - list with orientations per condition 
%plotCurve - 1: plot curves
%
%output - structure oriData with fields
%prefOri - preferred ori (0 - 180 deg)
%orthOri - orthogonal ori (0 - 180 deg)
%prefResp - preferred response (avg across dir)
%orthResp - orthogonal response (avg across dir)
%OSI - (pref-orth)/pref
%Lori - length vector (1-CircVar)
%Aori - angle of vector
%BW - bandwidth
%BWsmooth - smooth bandwidth, but after smoothing tuning curve first
%sigOri - result of Hotellings T-squared test


oriData=struct;

%% average rate and other shortcuts
avgRate=mean(FRate,1,'omitnan');
nRep=size(FRate,1);
if size(oridom,1)>1
    oridom=oridom';
end

%% preferred orientation - accounts for the fact that the maximum rate might
%%occur multiple times
[prefResp,prefCond]=max(avgRate);

if length(find(avgRate==prefResp))==1
    %only one max - can use it directly
    oriData.prefOri=oridom(prefCond);
else
    %use convolution to pick the maximum that is surrounded by
    %higher values; using square convolution
    tc=[avgRate(end) avgRate avgRate(1)];  %wrap around tuning curve first (need 1 on either side)
    peakIdx=find(avgRate==prefResp)+1; %indices for peaks (shift by 1 because of wrap around)
    rateConv=conv(tc,[1 1 1],'same'); %convolution
    mxConv=rateConv(peakIdx); %these are the convolution values for the peaks, pick the bigger one
    mx=peakIdx(mxConv==max(mxConv)); 
    mx=mx(1);%if there are multiple, we'll just take the first
    oriData.prefOri=oridom(mx-1);
end

%% preferred response
oriData.prefResp=avgRate(oridom==oriData.prefOri);

%% orthogonal orientation and response
oriData.orthOri=mod(oriData.prefOri+90,180);
oriData.orthResp=avgRate(oridom==oriData.orthOri);

%% OSI
oriData.OSI=(oriData.prefResp-oriData.orthResp)/oriData.prefResp;

%% vector length - better with rectified response
avgRateR=avgRate;
avgRateR(avgRateR<0)=0;
oriVec=sum(avgRateR.*exp(2*1i*oridom*pi/180)); 
oriData.Lori=abs(oriVec)/sum(avgRateR);
oriData.Aori=mod(angle(oriVec)*180/pi/2,180);

%% bandwidth - as in ringach 2002, with and without smoothing
wHan=hanning(3);
wHan=wHan/sum(wHan);
avgRateS=conv(avgRate,wHan,'same');

%criterion response: peak resp/sqrt(2)
critResp=oriData.prefResp/sqrt(2);
critRespSmooth=max(avgRateS)/sqrt(2);

%wrap around to avoid edge arfecats (need to do that before interpolation
tcWrap=[avgRate avgRate avgRate];
tcWrapS=[avgRateS avgRateS avgRateS];

%interpolate tuning function and orientation vector
nInterp=1000*3;
domInter=linspace(oridom(1)-180,oridom(end)+180,nInterp); %equivalent to 3x range

tcIP=griddedInterpolant(tcWrap); %generates object to be evaluated
tcInter=tcIP(linspace(1,length(tcWrap),nInterp));

tcIPS=griddedInterpolant(tcWrapS); %generates object to be evaluated
tcInterS=tcIPS(linspace(1,length(tcWrapS),nInterp));

%find crossings with criterion response
[~,~,cIdx]=zerocrossrate(tcInter,'Level',critResp);
cIdx(1)=0;
if isempty(cIdx) || sum(cIdx)<4
    oriData.BW=NaN;
else
    %find the ones closest to the peak
    crossOri=domInter(cIdx);
    cross1=max(crossOri(crossOri<oriData.prefOri));
    cross2=min(crossOri(crossOri>oriData.prefOri));
    oriData.BW=(cross2-cross1)/2;
end

[~,~,cIdxS]=zerocrossrate(tcInterS,'Level',critRespSmooth);
cIdxS(1)=0;
if isempty(cIdxS) || sum(cIdxS)<4
    oriData.BWSmooth=NaN;
else
    %find the ones closest to the peak
    crossOri=domInter(cIdxS);
    cross1S=max(crossOri(crossOri<oriData.prefOri));
    cross2S=min(crossOri(crossOri>oriData.prefOri));
    oriData.BWSmooth=(cross2S-cross1S)/2;
end


%% significance test
%based on Mazurek 2014
%compute length of ori vec for every repeat
oriRep=zeros(1,nRep);
FRateR=FRate;
FRateR(FRateR<0)=0;
FRateR(isnan(FRateR))=0;
for r=1:nRep
    %get all ori for one repeat
    ftmp=squeeze(FRateR(r,:));
    %vector sum on that repeat
    if sum(ftmp)==0
        oriRep(r)=0;
    else
        oriRep(r)=sum(ftmp.*exp(2*1i*oridom*pi/180))/sum(ftmp);
    end
end

%get x and y coordinates for each vector
x=real(oriRep);
y=-imag(oriRep);

%hotellings tsquare test against 0
[~,p]=hotellingt2test([x' y'],[0 0]);
oriData.sigOri=p;

%% plot (if selected)
if plotCurve==1
    figure
    subplot(1,2,1)
    plot([oridom oridom(1)+180],[avgRate avgRate(1)],'o-')
    hold on
    if ~isnan(oriData.BW)
       plot([cross1 cross2],[critResp critResp],'o--')
    end
    plot([oriData.prefOri oriData.prefOri],[0 oriData.prefResp],'g--')
    
    xlabel('Orientation (deg)')
    ylabel('Avg firing rate (Hz)')
    title({'Raw tuning curve';['OSI: ' num2str(oriData.OSI) ...
        ' Lori: ' num2str(oriData.Lori) ...
        ' BW: ' num2str(oriData.BW)]});

    subplot(1,2,2)
    plot([oridom oridom(1)+180],[avgRateS avgRateS(1)],'o-')
    hold on
    if ~isnan(oriData.BWSmooth)
        plot([cross1S cross2S],[critRespSmooth critRespSmooth],'o--')
    end
    plot([oriData.prefOri oriData.prefOri],[0 oriData.prefResp],'g--')
    
    xlabel('Orientation (deg)')
    ylabel('Avg firing rate (Hz)')
    title({'Smoothed tuning curve';...
        ['BWsmooth: ' num2str(oriData.BWSmooth)]});

end
