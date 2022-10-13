function compile_spkProps(spkpath,sortpath,animal,expname,probeID,jobStart,jobStop)
% input parameters:
% spkpath (character) - base path to folder containing SpikeFiles subfolder
% sortpath (character) - base path to folder containing spkSort and id files
% animal (character) - animal ID 
% expname (character) - experiment name 
% probeID (number) - probe number to process 
% jobStart (number) - start job ID 
% jobStop (number) - stop job ID 

%% compile spike waveform properties from spkinfo files
%extract sampling frequency
load(fullfile(sortpath,animal,expname,[expname '_id.mat'])); 

%build property vectors
wb = waitbar(0,'Compiling spike properties');
for job = jobStart:jobStop
    waitbar((job-jobStart)/(jobStop-jobStart),wb);
   
    %load spkinfo file 
    load(fullfile(spkpath,animal,expname,'SpikeFiles',[expname  '_j' num2str(job) '_p' num2str(probeID) '_spkinfo']),'spk');   
    
    if job==jobStart      
        AmpMin = spk.AmpMinDet;
        AmpMaxBefore = spk.AmpMaxBeforeDet;
        AmpMax = spk.AmpMaxDet;               
        PksBefore = spk.PksBeforeDet;
        Pks = spk.PksDet;       
        WidthI = (spk.WidthIDet/id.sampleFreq)*1000; %convert to ms
        En = spk.EnDet;       
    else               
        AmpMin = [AmpMin, spk.AmpMinDet];
        AmpMaxBefore = [AmpMaxBefore, spk.AmpMaxBeforeDet];
        AmpMax = [AmpMax, spk.AmpMaxDet];      
        PksBefore = [PksBefore, spk.PksBeforeDet];
        Pks = [Pks, spk.PksDet];
        WidthI = [WidthI, spk.WidthIDet];
        En = [En, spk.EnDet];
    end         
    clear spk;  
    
end
WidthI = (WidthI/id.sampleFreq)*1000; %convert from samples to ms
MaxMinRatio = abs(AmpMax./AmpMin); %ratio of peak to trough
close(wb);

%% create spkProps struct containing the mean and standard deviation of unit waveform properties 
%extract unit info
load(fullfile(sortpath,animal,expname,[expname '_p' num2str(probeID) '_spkSort.mat'])); 

%build spkProps
for unit = 1:length(spkSort.unitinfo)   
    spkProps.avgEn(unit) = mean(En(spkSort.unitid==unit)); 
    spkProps.avgAmpMin(unit) = mean(AmpMin(spkSort.unitid==unit)); 
    spkProps.avgAmpMax(unit) = mean(AmpMax(spkSort.unitid==unit)); 
    spkProps.avgAmpMaxBefore(unit) = mean(AmpMaxBefore(spkSort.unitid==unit)); 
    spkProps.avgMaxMinRatio(unit) = mean(MaxMinRatio(spkSort.unitid==unit));
    spkProps.avgPks(unit) = mean(Pks(spkSort.unitid==unit)); 
    spkProps.avgPksBefore(unit) = mean(PksBefore(spkSort.unitid==unit)); 
    spkProps.avgWidthI(unit) = mean(WidthI(spkSort.unitid==unit));     
    
    spkProps.stdEn(unit) = std(En(spkSort.unitid==unit)); 
    spkProps.stdAmpMin(unit) = std(AmpMin(spkSort.unitid==unit)); 
    spkProps.stdAmpMax(unit) = std(AmpMax(spkSort.unitid==unit));
    spkProps.stdAmpMaxBefore(unit) = std(AmpMaxBefore(spkSort.unitid==unit));
    spkProps.stdMaxMinRatio(unit) = std(MaxMinRatio(spkSort.unitid==unit));
    spkProps.stdPks(unit) = std(Pks(spkSort.unitid==unit)); 
    spkProps.stdPksBefore(unit) = std(PksBefore(spkSort.unitid==unit)); 
    spkProps.stdWidthI(unit) = std(WidthI(spkSort.unitid==unit)); 
end
 
%% nest spkProps within spkSort and save 
spkSort.spkProps = spkProps;
save(fullfile(sortpath,animal,expname,[expname '_p' num2str(probeID) '_spkSort.mat']),'spkSort');





