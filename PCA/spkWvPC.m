%Update spike info files with scores (weights) of waveforms on 1st 3 PCs
%from PCA (wavePCA.m => FEXX_000_000_pX_pca.mat)

function spkWvPC(animalID,unitID,expID,probeID,job,dataFold)
% INPUTS:
% animalID = animal ID (string, ex: ‘febe0’)
% unitID = unit ID (string, ex: ‘000’)
% expID = experiment ID (string, ex: ‘000’)
% probeID = which probe in this experiment (number)
% job = the index in the parfor loop ('j') or whatever job you are running
%       this on
% dataFold = path to folder storing data to be processed (string, ex:
%            'C:\Users\brand\Documents\data')

disp(['job:' num2str(job)])
fileBase = [animalID '_u' unitID '_' expID];

load(fullfile(dataFold,animalID,fileBase,[fileBase '_p' num2str(probeID) '_pca.mat']))
coeff = PC.coeff;
% score = PC.score;
% latent = PC.latent;
% explained = PC.explained;
mu = PC.mu;
% wvs = PC.wvfrms;

load(fullfile(dataFold,animalID,fileBase,'spikeFiles',[fileBase '_j' num2str(job) '_p' num2str(probeID) '_spkinfo.mat']))
load(fullfile(dataFold,animalID,fileBase,'spikeFiles',[fileBase '_j' num2str(job) '_p' num2str(probeID) '_spike.mat']))

nSpksJob = length(spk.spkTimesDet);
jobWvfrms = zeros(nSpksJob,length(mu));
% nbrWvfrms{nSpksJob,1} = [];
% nbrWvfrmsA{nSpksJob,1} = [];
% PC1All = cell(size(nbrWvfrms));
% PC2All = cell(size(nbrWvfrms));
% PC3All = cell(size(nbrWvfrms));

s=0;
for ch = 1:length(spikeData) %loop through all channels for this job and collect waveforms in variable 'jobWvfrms'

    if isnan(spikeData(ch).spikeTimes) %if there aren't any spikes on this channel (in this job) move onto the next
        continue
    end
    nSpkCh = length(spikeData(ch).spikeTimes); %number of spikes that occur on this channel
    
    jobWvfrms(1+s:nSpkCh+s,:) = spikeData(ch).Wvfrms(:,:,1); %fill in the next 'nSpkCh' rows of 'jobWvfrms' with the spikes from this channel
    
    s = s+nSpkCh; %increase counter by the number of spikes you've filled in

%     for spike = 1:nSpkCh
%         s=s+1;
%         nbrWvfrms{s} = squeeze(spikeData(ch).Wvfrms(spike,:,:))';
%         nbrWvfrmsA{s} = alignWvs(nbrWvfrms{s});
%         
%         spk.PC1All{s} = (nbrWvfrmsA{s}-mu)*coeff(:,1);
%         spk.PC2All{s} = (nbrWvfrmsA{s}-mu)*coeff(:,2);
%         spk.PC3All{s} = (nbrWvfrmsA{s}-mu)*coeff(:,3);
%     end
    
end

if isempty(jobWvfrms)
    disp(['job ' num2str(job) ' has no spikes'])
else
    score = (jobWvfrms-mu)*coeff;

    spk.PC1Det = score(:,1)'; %store projections of spikes in first 3 PCs as fields in structure 'spk'
    spk.PC2Det = score(:,2)';
    spk.PC3Det = score(:,3)';

    cd([dataFold '/' animalID '/' fileBase '/SpikeFiles'])
    save([fileBase '_j' num2str(job) '_p' num2str(probeID) '_spkinfo.mat'],'spk')
end

end