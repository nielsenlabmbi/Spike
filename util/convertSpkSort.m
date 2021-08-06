function convertSpkSort(spkFile)

load(spkFile); %geenrates spkSort

spkSort.info.probeid=spkSort.probeid;
spkSort.info.jobStart=spkSort.jobStart;
spkSort.info.jobStop=spkSort.jobStop;
spkSort.info.name=spkSort.name;
spkSort.info.artRej=spkSort.artRej;
spkSort.info.date=spkSort.date;
spkSort.info.expname=spkSort.expname;

spkSort=rmfield(spkSort,{'probeid','jobStart','jobStop','name','artRej','date','expname'});

if isfield(spkSort,'artNCh')
    spkSort.info.artNCh=spkSort.artNCh;
    spkSort.info.artThreshold=spkSort.artThreshold;
    spkSort=rmfield(spkSort,{'artNCh','artThreshold'});
end

save(spkFile,'spkSort')