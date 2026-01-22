function projectOut=getProjectFiles(project,includeOnly,varargin)

%input:
%project: string defining the project in the database
%includeOnly: all files or only those with includeFlag=1
%additional arguments: parameters to retrieve from database about every
%project file

%output
%projectOut: table with predefined plus added columns
%predefined:
%experimentId: FEXXX
%unitNr: unit
%experimentNr: experiment
%probeId: probe id
%datatype: SU or MU
%fileBase: full filename without probe
%filePhys: full filename with probe

%open database
connSQL=database('ephysDatabase','','');
if ~isempty(connSQL.Message)
    %needs password, most likely - assumption is this is stored as a secret
    %called 'dbPassword'
    passw=getSecret('dbPassword');
    connSQL=database('ephysDatabase','nielsenlab',passw);

end


%get data for project
selquery=['SELECT * FROM tblproject WHERE projName="' project '"'];
data=select(connSQL,selquery);

if isempty(data)
    errordlg('Project does not exist!','Read error');
    projectOut=table;
    close(connSQL);
    return;
end

%return the relevant information
projectOut=table;

if includeOnly==1
    includeIdx=(data.includeFlag==1);
else
    includeIdx=ones(height(data),1);
end

projectOut.experimentId=data.experimentId(includeIdx);
projectOut.unitNr=data.unitNr(includeIdx);
projectOut.experimentNr=data.experimentNr(includeIdx);
projectOut.probeId=data.probeId(includeIdx);
projectOut.datatype=data.datatype(includeIdx);
projectOut.sortSuffix=data.sortSuffix(includeIdx);

%add filenames for easier use
projectOut.fileBase=strcat(projectOut.experimentId,'_u',projectOut.unitNr,...
    '_',projectOut.experimentNr);
projectOut.filePhys=strcat(projectOut.experimentId,'_u',projectOut.unitNr,...
    '_',projectOut.experimentNr,'_p',num2str(projectOut.probeId));

%now add extra if required
if nargin>2
    %get column names so that we can figure out which table to search in
    sqlColumns=sqlfind(connSQL,'');
    tblA=sqlColumns.Columns{strcmp(sqlColumns.Table,'tblanimal')};
    tblE=sqlColumns.Columns{strcmp(sqlColumns.Table,'tblexperiment')};
    tblR=sqlColumns.Columns{strcmp(sqlColumns.Table,'tblrecording')};
    tblM=sqlColumns.Columns{strcmp(sqlColumns.Table,'tblmanipulation')};
    columnNames=[tblA tblE tblR tblM];

    %vector of table names to make the loop below easier
    tblNames(1:length(tblA))={'tblanimal'};
    tblNames(end+1:end+length(tblE))={'tblexperiment'};
    tblNames(end+1:end+length(tblR))={'tblrecording'};
    tblNames(end+1:end+length(tblM))={'tblmanipulation'};

    %search criteria are different for every table
    searchCrit(1:length(tblA))={'experimentId'};
    searchCrit(end+1:end+length(tblE))={'experimentKey'};
    searchCrit(end+1:end+length(tblR))={'recordingKey'};
    searchCrit(end+1:end+length(tblM))={'experimentKey'};

    %need to make entries with keys into cells
    data=convertvars(data,{'experimentKey','recordingKey'},'string');
    data=convertvars(data,{'experimentKey','recordingKey'},'cellstr');

    for i=1:length(varargin)
        idx=find(strcmp(columnNames,varargin{i}));
       
        count=1;
        for j=1:height(data) %need data to get the relevant search keys
            if includeIdx(j)==1
                selquery=['SELECT ' varargin{i} ' FROM ' tblNames{idx} ...
                    ' WHERE ' searchCrit{idx} '="' data.(searchCrit{idx}){j} '"'];

                dTmp=select(connSQL,selquery);
                if isempty(dTmp)
                    disp('Warning; could not locate information in DB');
                else
                    projectOut.(varargin{i})(count)=dTmp.(varargin{i});
                end

                count=count+1;

            end %if
        end %for j
    end
end %nargin


%close database
close(connSQL);