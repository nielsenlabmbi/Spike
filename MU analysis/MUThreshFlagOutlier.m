function trialExclude=MUThreshFlagOutlier(MUThresh,MUThreshInfo,plotOpt)

%flag trials as outliers; outliers based on probe-wide metric

%shortcuts
nTrial=length(MUThreshInfo.triallist);
nCh=length(MUThresh);


%compute summary metric for every trial
spkMat=[MUThresh.baseNspk]; %all sites for all trials
spkMat=reshape(spkMat,nTrial,nCh);
trialSum=mean(spkMat,2);


%find outliers
nScl = 5;
scale = -1/(sqrt(2)*erfcinv(3/2));
scaledMad = scale*median(abs(trialSum-median(trialSum)));
% outProbe = isoutlier(trialSum);
% outProbe = trialSum > mean(trialSum)+(nScl*std(trialSum));
outProbe = trialSum > median(trialSum)+(nScl*scaledMad);
trialExclude=outProbe(:); %outliers = 1

%plot if selected
if plotOpt==1
    idxOut=find(trialExclude==1);

    figure
    t=tiledlayout(1,2,'TileSpacing','None');
    nexttile
    imagesc(spkMat);
    colormap('jet')
    xlabel('Ch')
    ylabel('Trials')
    if ~isempty(idxOut)
        hold on
        for i=1:length(idxOut)
            plot([1 nCh],[idxOut(i) idxOut(i)],'r-')
        end
    end
    
    nexttile
    plot(trialSum(:),[1:nTrial],'k-o')
    set(gca,'YLim',[0.5 nTrial+0.5],'YTickLabel',[],'YDir','reverse')
    xlabel('Sum Ch')

    xline(mean(trialSum),'c-','LineWidth',2)
    for s = 1:nScl
        xline(mean(trialSum)+(s*std(trialSum)),'c--')
    end
    xline(median(trialSum),'m-','LineWidth',2)
    for s = 1:nScl
        xline(median(trialSum)+(s*scaledMad),'m--')
    end

    if ~isempty(idxOut)
        hold on
        xL=get(gca,'XLim');
        for i=1:length(idxOut)
            plot([xL(1) xL(2)],[idxOut(i) idxOut(i)],'r-')
        end
    end

end
