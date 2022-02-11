%PCA on spike waveforms

function wavePCA(animalID,unitID,expID,probeID,nJobs,propSpks,dataFold)
% INPUTS:
% animalID = animal ID (string, ex: ‘febe0’)
% unitID = unit ID (string, ex: ‘000’)
% expID = experiment ID (string, ex: ‘000’)
% probeID = which probe in this experiment (number)
% nJobs = number of jobs SpikeFiles are broken up into (number)
% propSpks = proportion of spikes used to calculate PCs (number)
% dataFold = path to folder storing data to be processed (string ex:
%            'C:\Users\brand\OneDrive - Johns Hopkins\Documents\data')
int = 0;
align = 0;

fileBase = [animalID '_u' unitID '_' expID];

%% Make data matrix
%Collect waveforms
waves = []; %nxp matrix where n is number of data points (nSpks total) and p is number of waveform samples
jobId = [];
chId = [];
flagJob = [];
nSpksJob = zeros(1,nJobs);
nSpksSam = zeros(1,nJobs);
for job = 0:nJobs-1 %go through each job and collect all spike waveforms in 'waves'

    disp(['job' num2str(job)])
    load(fullfile(dataFold,animalID,fileBase,'SpikeFiles',[fileBase '_j' num2str(job) '_p' num2str(probeID) '_spike.mat']))
    load(fullfile(dataFold,animalID,fileBase,'SpikeFiles',[fileBase '_j' num2str(job) '_p' num2str(probeID) '_spkInfo.mat']))
    
    if ~isfield(spk,'spkTimesDet') %if there is no spike data for this job, move onto the next
            flagJob = [flagJob job];
            disp(['     job ' num2str(job) ' has no spikes'])
            continue
    end
    
    nSpksJob(job+1) = length(spk.spkTimesDet); %number of spikes (across channels) that occured in this job 
    
    for chan = 1:length(spikeData)
        
        chWvfrms = spikeData(chan).Wvfrms(:,:,1);
        if isnan(chWvfrms) %if there are no spikes of this channel (in this job), move onto next channel
            continue
        end
        [nSpksCh,nSam] = size(chWvfrms);
        samID = sort(randsample(nSpksCh,round(nSpksCh*propSpks)));
        nSpks = length(samID);
        nSpksSam(job+1) = nSpksSam(job+1)+nSpks; %number of spikes (across channels) that we are sampling from this job
        samWvfrms = chWvfrms(samID,:);
        
        chId = [chId repmat(chan,1,nSpks)];
        waves = [waves;samWvfrms];
        
    end
    
    jobId = [jobId repmat(job,1,nSpksSam(job+1))];

end

%Interpolate waveforms 
if int == 1
    intMethod = 'pchip'; dInterp=0.1;
    
    % wvsTmp=griddedInterpolant(waves,intMethod);
    % wvsInt=wvsTmp({[1:size(waves,1)],[1:dInterp:nSam]});
    
    vq = [1:dInterp:nSam];
    wvsInt = interp1(waves',vq,intMethod); wvsInt = wvsInt';
    
    waves = wvsInt;
end 

%Align waveforms
if align == 1
    M = min(waves,[],2); %find the minimum for each waveform
    for w = 1:length(M)
        mID(w) = find(waves(w,:)==M(w)); %find which sample ID the minimum occurs at for each waveform
    end

    if length(unique(mID))>1 %if not all waveforms are aligned at minimum
        disp('aligning waveforms to minimum at center')
        wavesA = alignWvs(waves); %center all waveforms at their minimum
        waves = wavesA;
    end
end

%% PCA
data = waves;

% [coeff,score,latent,tsquared,explained,mu] = pca(data);

mu = mean(data);Xo = data-mu;
[coeff,latent] = eigs(cov(Xo),size(data,2));
latent = diag(latent);explained = (latent./sum(latent))*100;
score = Xo*coeff;

PC = struct('fileName',[fileBase '_p' num2str(probeID) '_pca'],'coeff',coeff,'latent',latent,...
    'explained',explained,'score',score,'wvfrms',data,'mu',mu,'propSpks',propSpks);
if int == 1
    PC.dInterp = dInterp;
end

%Plot PCA data
figure;
subplot(2,2,1)
plot3(PC.score(:,1),PC.score(:,2),PC.score(:,3),'.')
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

subplot(2,2,2)
plot(PC.coeff(:,1),'LineWidth',2);hold on
plot(PC.coeff(:,2),'LineWidth',2)
plot(PC.coeff(:,3),'LineWidth',2)
plot(PC.coeff(:,4),'LineWidth',2)
legend(['PC1;explained var:' num2str(PC.explained(1)) '%'],['PC2;explained var:' num2str(PC.explained(2)) '%'],...
    ['PC3;explained var:' num2str(PC.explained(3)) '%'],['PC4;explained var:' num2str(PC.explained(4)) '%'])
axis tight

subplot(2,2,3)
plot(cumsum(PC.explained),'ro');hold on
plot([0,length(PC.explained)],[90,90],'--','Color','k')
minPC = find(cumsum(PC.explained)>90,1);
plot([minPC minPC],[0 90],'--','Color','k')
xlabel('PC')
ylabel('cumulative % explained var')
axis tight

subplot(2,2,4);hold on
id = randsample(size(PC.wvfrms,1),5000);
plot(PC.wvfrms(id,:)','LineWidth',1)
plot(PC.mu,'Color','k','LineWidth',2)
axis tight

% figure;
% z1 = 1; z2 = 2;
% subplot(2,2,1)
% plot(PC.score(:,z1),PC.score(:,z2),'.')
% xlabel(['PC' num2str(z1)]);ylabel(['PC' num2str(z2)])
% subplot(2,2,2)
% histogram(PC.score(:,z2));set(gca, 'XDir','reverse');camroll(-90)
% title(['PC' num2str(z2)])
% ylim([0 3000])
% linkaxes([subplot(2,2,1),subplot(2,2,2)],'y')
% subplot(2,2,3)
% histogram(PC.score(:,z1));set(gca, 'YDir','reverse')
% title(['PC' num2str(z1)])
% linkaxes([subplot(2,2,1),subplot(2,2,3)],'x')

%% SAVE
%Save pca inputs/outputs in 'fileBase'_pca.mat and scores for each job in
%respective 'fileBase'_jX_pX_spkInfo.mat file

cd(fullfile(dataFold,animalID,fileBase))
save([fileBase '_p' num2str(probeID) '_pca'],'PC')

end



















%Interpolate and align waveforms at new minimums

% sfreq = id.sampleFreq; %sample frequency 
% sampfactor = 5; %upsample factor 
% minshift = 16 + 15*(sampfactor-1); %new index for old minimum after upsampling
% minsearch = [minshift-(sampfactor-1):minshift+(sampfactor-1)]; %search space for new minimum after upsampling
% x = [1:31]; %old waveforms = 31 samples, with minimum at 16
% xx = [1:1/sampfactor:31]; %new sample points w/ interpolation 
% 
% wavesInt = spline(x,waves,xx); %interpolate
% [minamp minidx] = min(wavesInt(:,minsearch),[],2);
% minidx = minsearch(minidx); 
% minoffset = minidx - minshift; 
% for w = 1:length(wavesInt(:,1))
%     wvfshift = circshift(wavesInt(w,:),-1*minoffset(w),2);
%     wavesNew(w,:) = downsample(wvfshift,sampfactor);
% end
