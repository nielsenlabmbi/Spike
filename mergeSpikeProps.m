function mergeSpikeProps(expFolder,animalID,unitID,expID,probeId,jobIDstart,jobIDstop,name)

%input parameters:
% animalID - animal ID (string)
% unitID - unit ID (string)
% expID - experiment ID (string)
% expFolder - base folder for experiments (string)
% jobIDstart - job ID of first file to merge
% jobIDstop - job ID of last file to merge
% name - name of person executing the merge
%
% output: spk file
% also updates id file


expname=[animalID '_u' unitID '_' expID];



f=waitbar(0,'Merging...');
for i=jobIDstart:jobIDstop
    waitbar((i-jobIDstart)/(jobIDstop-jobIDstart),f);
    
    filename=fullfile(expFolder,animalID,expname,'SpikeFiles',[expname '_j' num2str(i) '_p' num2str(probeId) '_spkinfo']);

    load(filename); %generates spk
    dataVar=fieldnames(spk);
    
    for n=1:length(dataVar)
       if i==jobIDstart 
           spkInfo.(dataVar{n})=spk.(dataVar{n});
       else
           if ~strcmp(dataVar{n},'expname') & ~strcmp(dataVar{n},'probeId')
               spkInfo.(dataVar{n})=[spkInfo.(dataVar{n}) spk.(dataVar{n})];
           end
       end
      
    end
    
end

%add which jobs were merged (just in case)
spkInfo.jobRange=[jobIDstart jobIDstop];

%add info about probes

outname=fullfile(expFolder,animalID,expname,[expname '_p' num2str(probeId) '_spk']);
save(outname,'spkInfo','-v7.3');


%load id file, update id and save
outname=fullfile(expFolder,animalID,expname,[expname  '_id']);
load(outname);
id.mergeSpikeProps(probeId).name=name;
id.mergeSpikeProps(probeId).date=date;
save(outname,'id');


close(f);