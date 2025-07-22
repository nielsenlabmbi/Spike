function compile_spkPos(spkpath,animal,unit,exp,probeID)
% computes position information for all units in a spksort file
% (independent of classification)
%
% input parameters:
% spkpath (character) - base path to folder containing SpikeFiles subfolder
% animal (character) - animal ID 
% unit (character) - unit
% exp (character) - experiment
% probeID (number) - probe number to process 


%% compile position info from spkInfo files

%load spksort and id file
expname=[animal '_u' unit '_' exp];
load(fullfile(spkpath,animal,expname,[expname '_p' num2str(probeID) '_spkSort.mat'])); 
load(fullfile(spkpath,animal,expname,[expname '_id.mat'])); 


%number of units
nrUnit=length(spkSort.unitinfo);

%figure out job file for each spike
jobBin=discretize(spkSort.spktimes,spkSort.info.jobEdges);
jobBin=jobBin-1; %jobs start with 0
uniqueJobs=unique(jobBin); %may not need to read all jobs


%accumulators
comXMin=zeros(1,nrUnit);
comZMin=zeros(1,nrUnit);
comXEn=zeros(1,nrUnit);
comZEn=zeros(1,nrUnit);
Nspk=zeros(1,nrUnit);

%loop through jobs
wb = waitbar(0,'Compiling spike properties');
for job = uniqueJobs
    waitbar(job/max(uniqueJobs),wb);
   
    %load spkinfo file 
    load(fullfile(spkpath,animal,expname,'SpikeFiles',[expname  '_j' num2str(job) '_p' num2str(probeID) '_spkinfo']));   
    for u=1:nrUnit
        idxIn=find(spkSort.unitid==u & jobBin==job); %find the spikes for this unit in this job
        if ~isempty(idxIn)
            spkPar=[spkSort.spktimes(idxIn);spkSort.detCh(idxIn)]'; %need to be able to identify the correct spikes in the jobs file
            idxOut=ismember([spk.spkTimesDet;spk.detCh]',spkPar,'rows'); %find the spikes
            comXMin(u)=comXMin(u)+sum(spk.comXMinDet(idxOut)); %accumulate
            comZMin(u)=comZMin(u)+sum(spk.comZMinDet(idxOut));
            comXEn(u)=comXEn(u)+sum(spk.comXEnDet(idxOut));
            comZEn(u)=comZEn(u)+sum(spk.comZEnDet(idxOut));
            Nspk(u)=Nspk(u)+sum(idxOut);
        end
    end
    
    clear spk;  
end
close(wb);

%transfer
spkSort.spkProps.comXMin=comXMin./Nspk;
spkSort.spkProps.comZMin=comZMin./Nspk;
spkSort.spkProps.comXEn=comXEn./Nspk;
spkSort.spkProps.comZEn=comZEn./Nspk;

%detCh is saved in spkSort, so no jobs loop needed, but we need probe info

for u=1:nrUnit
    chidx=spkSort.detCh(spkSort.unitid==u);
    chdist=accumarray(chidx',1,[id.probes(probeID).nChannels 1]);
    [~,mxCh]=max(chdist);

    spkSort.spkProps.detChMain(u)=mxCh;
    spkSort.spkProps.detChMainX(u)=id.probes(probeID).x(mxCh);
    spkSort.spkProps.detChMainZ(u)=id.probes(probeID).z(mxCh);
end

 
%% save 
save(fullfile(spkpath,animal,expname,[expname '_p' num2str(probeID) '_spkSort.mat']),'spkSort');





