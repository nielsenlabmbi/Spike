function splitThreshold(oldFilename)
%splits threshold files into 2 probes

load(oldFilename); %generates structure thresholding


if length(thresholding.thresholds)<=64 %one probe only
    newFilename=replace(oldFilename,'_thres','_p1_thres');
    save(newFilename,'thresholding');
else
    tmp=thresholding;
    for p=1:2
        thresholding.thresholds=tmp.thresholds((p-1)*64+1:p*64);
        thresholding.badChannels=tmp.badChannels((p-1)*64+1:p*64);
        thresholding.chkChannels=tmp.chkChannels((p-1)*64+1:p*64);
         newFilename=replace(oldFilename,'_thres',['_p' num2str(p) '_thres']);
         save(newFilename,'thresholding');
    end
end
    
