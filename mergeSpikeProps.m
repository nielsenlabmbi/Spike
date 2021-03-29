function mergeSpikeProps(expFolder,animalID,unitID,expID,jobIDstart,jobIDstop)

% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% expFolder - base folder for experiments (string)
% jobIDstart - job ID of first file to merge
% jobIDstop - job ID of last file to merge


expname=[animalID '_u' unitID '_' expID];

f=waitbar(0,'Merging...');
for i=jobIDstart:jobIDstop
    waitbar((i-jobIDstart)/(jobIDstop-jobIDstart),f);
    
    filename=fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(i) '_spkinfo']);

    load(filename); %generates spk
    dataVar=fieldnames(spk);
    
    for n=1:length(dataVar)
       if i==jobIDstart 
           spkInfo.(dataVar{n})=spk.(dataVar{n});
       else
           if ~strcmp(dataVar{n},'expname') && ~strcmp(dataVar{n},'dateProps')
               spkInfo.(dataVar{n})=[spkInfo.(dataVar{n}) spk.(dataVar{n})];
           end
       end
      
    end
    
end

outname=fullfile(expFolder,animalID,expname,[expname  '_spk']);
save(outname,'spkInfo');
close(f);