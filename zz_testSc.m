clc; clear; close all;

load('/Users/ramanujan/Downloads/meaa1_u000_011_Spikes.mat')

cc = [];
for ii=1:16
    for jj=(ii+1):16
        cc = [cc; ii jj];
    end
end

unit1 = 0;
unit2 = 1;

clf;
set(gcf,'pos',[182,1,804,954],'color','w');
ha = tight_subplot(12,10,0.005,0.005,0.005);

for ii=1:120
    u1 = Properties(cc(ii,:),idk==unit1);
    u2 = Properties(cc(ii,:),idk==unit2);
    
    plot(ha(ii),u1(2,:),u1(1,:),'k.'); hold(ha(ii),'on');
    plot(ha(ii),u2(2,:),u2(1,:),'r.')
    grid(ha(ii),'on');
    axis(ha(ii),'square');
%     axis(ha(ii),'equal');
    set(ha(ii),'xtick',[],'ytick',[]);
end
    