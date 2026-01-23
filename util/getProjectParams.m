function pOut=getProjectParams(project,includeOnly,addAnalyzer,varargin)

%input:
%project: string defining the project in the database
%includeOnly: all files or only those with includeFlag=1
%addAnalyzer: use analyzer files to fill in missing data?
%additional arguments: parameters to retrieve from database 

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

%project='DSdev';
%includeOnly=1;


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
    pOut=table;
    close(connSQL);
    return;
end

%return the relevant information
pOut=table;

if includeOnly==1
    data=data(data.includeFlag==1,:);
end

pOut.experimentId=data.experimentId;
pOut.unitNr=data.unitNr;
pOut.experimentNr=data.experimentNr;

%reduce to unique
[pOut,idx]=unique(pOut,'rows','stable');
pOut.experimentKey=data.experimentKey(idx);

%loop through all entries
idxMissing=zeros(height(pOut),length(varargin));
for i=1:height(pOut) 
    %find all entries in parameter table
    selquery=['SELECT * FROM tblparams WHERE ' ...
        'experimentKey="' num2str(pOut.experimentKey(i)) '"'];
    dTmp=select(connSQL,selquery);

    %dtmp is a table, find relevant entries
    for j=1:length(varargin)
        idx=find(strcmp(dTmp.pname,varargin{j}));

        if ~isempty(idx)
            pOut.(varargin{j})(i)=dTmp.pval(idx);
        else
            idxMissing(i,j)=1;
        end

    end
end 


%close database
close(connSQL);

%if selected, add missing information from analyzer files
if addAnalyzer==1
    for i=1:height(pOut) 
        if any(idxMissing(i,:))
            %read analyzer file
            fname=[pOut.experimentId{i} '_u' pOut.unitNr{i} '_' pOut.experimentNr{i} '.analyzer'];
            fnameFull=fullfile('z:\ephysnew\analyzer',pOut.experimentId{i},fname);
            load(fnameFull,'-mat'); %generates Analyzer
            
            for j=1:length(varargin)
                pval=getparam(varargin{j},Analyzer);
                pOut.(varargin{j})(i)={num2str(pval)};
            end

        end
    end
end