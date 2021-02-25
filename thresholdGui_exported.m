classdef thresholdGui_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ThresholdGUIUIFigure    matlab.ui.Figure
        GridLayout              matlab.ui.container.GridLayout
        LeftPanel               matlab.ui.container.Panel
        panel_location          matlab.ui.container.Panel
        button_getExpt          matlab.ui.control.Button
        button_save             matlab.ui.control.Button
        button_load             matlab.ui.control.Button
        ProbeAreaLabel_2        matlab.ui.control.Label
        dropdown_probe2         matlab.ui.control.DropDown
        ProbeAreaLabel          matlab.ui.control.Label
        dropdown_probe1         matlab.ui.control.DropDown
        dropdown_area1          matlab.ui.control.DropDown
        dropdown_area2          matlab.ui.control.DropDown
        panel_fileInfo          matlab.ui.container.Panel
        table_fileInfo          matlab.ui.control.Table
        RightPanel              matlab.ui.container.Panel
        panel_signal            matlab.ui.container.Panel
        axis_snippet            matlab.ui.control.UIAxes
        ThresholdSliderLabel    matlab.ui.control.Label
        slider_threshold        matlab.ui.control.Slider
        button_randomSnippet    matlab.ui.control.Button
        button_nextChan         matlab.ui.control.Button
        button_prevChan         matlab.ui.control.Button
        button_rejectChan       matlab.ui.control.Button
        axis_spikePreview       matlab.ui.control.UIAxes
        panel_extract           matlab.ui.container.Panel
        radio_compSelect        matlab.ui.container.ButtonGroup
        radio_local             matlab.ui.control.RadioButton
        radio_remote            matlab.ui.control.RadioButton
        remoteipEditFieldLabel  matlab.ui.control.Label
        text_ip                 matlab.ui.control.EditField
        button_checkConn        matlab.ui.control.Button
        button_copyData         matlab.ui.control.Button
        button_extract          matlab.ui.control.Button
        button_firstSnippet     matlab.ui.control.Button
        button_lastSnippet      matlab.ui.control.Button
        text_channelNum         matlab.ui.control.NumericEditField
        panel_probe             matlab.ui.container.Panel
        axis_probe              matlab.ui.control.UIAxes
        label_selectedProbe     matlab.ui.control.Label
        button_nextProbe        matlab.ui.control.Button
        button_prevProbe        matlab.ui.control.Button
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    properties (Access = private)
        settings
        exptName
        probe
        sampleFreq
        nChannels
        thresholds
        badChannels
        nTotalSamples
        snippet
        butterParams
        
        currentProbeNum
        currentChannel
        currentTc
        
        thresholdLineHandle
        signalLineHandle
        channelMarkerHandle
    end
    
    methods (Access = private)
        function loadSettings(app)
            if ~exist('Settings.mat','file')
                a = inputdlg('Enter your name:','Name',1,{'Nielsen Lab'});
                name = a{1};
                a = inputdlg('Enter your email:','Email',1,{'nielsenlabmbi@gmail.com'});
                email = a{1};
                remoteIP = '172.30.6.202';
                doRemote = 1;
                defaultFilter.highPass = 250;
                defaultFilter.lowPass = 5000;
                expFolder = '/Users/ramanujan/Downloads';
                spikeFileEnding = '_spikes.mat';
                spikesFolder = 'C:\Users\nielsenlab\Documents\MATLAB\ephysAnalysis\out\spikeStage';
                summaryFile = 'C:\Users\nielsenlab\Documents\MATLAB\ephysAnalysis\out\summary\Ephys Data.xlsx';
                dataFileFolder = 'C:\Users\nielsenlab\Documents\MATLAB\ephysAnalysis\out\dataStage';
                dataFileEnding = '_data.mat';
                analyzerFileFolder = 'Z:\EphysNew\analyzer';
                analyzerFileEnding = '.analyzer';
                save('Settings.mat','name','email','remoteIP','doRemote','defaultFilter','expFolder','spikesFolder','spikeFileEnding','analyzerFileEnding','analyzerFileFolder','dataFileEnding','dataFileFolder','summaryFile','summaryFile');
            end
            loadedFile = load('Settings.mat');
            app.settings = loadedFile;
        end
        
        function refreshFileInfo(app)
            % display Settings
            % display FileProps
            
            data = {
                    'expt name' '';...
                    'exp loc' app.settings.expFolder;...
                    'n bad chans' '0';...
                    'probes' 'XXX,YYY';...
                    
                    'name' app.settings.name;...
                    'email' app.settings.email;...
                    
                    'high pass' num2str(app.settings.defaultFilter.highPass);...
                    'low pass' num2str(app.settings.defaultFilter.lowPass);...
                    'sample freq' 'XXX';...
                    
                    'analyzer loc' app.settings.analyzerFileFolder;...
                    'spike loc' app.settings.spikesFolder;...
                    'data loc' app.settings.dataFileFolder;...
                   };
                    
            app.table_fileInfo.Data = cell2table(data);
        end
        
        function refreshRemoteInfo(app)
            if ~isfield(app.settings,'doRemote')
                app.settings.doRemote = false;
                app.settings.remoteIP = '';
            end
            
            if app.settings.doRemote
                app.radio_compSelect.SelectedObject = app.radio_remote;
                app.text_ip.Value = app.settings.remoteIP;
                app.text_ip.Enable = true;
                app.button_checkConn.Enable = true;
                app.button_copyData.Enable = true;
            else
                app.radio_compSelect.SelectedObject = app.radio_local;
                app.text_ip.Value = app.settings.remoteIP;
                app.text_ip.Enable = false;
                app.button_checkConn.Enable = false;
                app.button_copyData.Enable = false;
            end
            
            doRemote = app.settings.doRemote;
            remoteIP = app.settings.remoteIP;
            
            save('Settings.mat','doRemote','remoteIP','-append');
        end
        
        function refreshSnippet(app,getNew,pos)
            if getNew
                app.getNewSnippet(pos);
            end
            
            cla(app.axis_snippet);
            app.currentTc = filter(app.butterParams.b1, app.butterParams.a1, app.snippet(app.currentChannel,:));
            app.signalLineHandle = plot(app.axis_snippet,linspace(0,1000,app.sampleFreq),app.currentTc,'k');
            hold(app.axis_snippet,'on')
            
            if app.badChannels(app.currentChannel)
                app.button_rejectChan.Text = 'Accept Channel';
                set(app.axis_snippet,'color',[0.6 0.3 0.3])
            else
                app.button_rejectChan.Text = 'Reject Channel';
                set(app.axis_snippet,'color','w')
            end
            app.plotThreshold()
            app.highlightChannel()
            app.previewSpikes()
        end
        
        function getNewSnippet(app,pos)
            dataFileName = [app.settings.expFolder '/' app.exptName '/' app.exptName '_amplifier.dat'];
            DataFile = fopen(dataFileName,'r');
            
            frewind(DataFile);
            if strcmp(pos,'first')
                startSample = 50;
            elseif strcmp(pos,'last')
                startSample = app.nTotalSamples - 120*app.sampleFreq;
            else % random
                startSample = randi(app.nTotalSamples - 120*app.sampleFreq);
            end
            fseek(DataFile,2*startSample*app.nChannels,'bof');

            app.snippet = fread(DataFile, [app.nChannels app.sampleFreq], 'int16');
        end
        
        function plotThreshold(app)
            if ~isempty(app.thresholdLineHandle)
                if isvalid(app.thresholdLineHandle)
                    delete(app.thresholdLineHandle);
                end
            end
            app.thresholdLineHandle = line(app.axis_snippet,[0 1000],[app.thresholds(app.currentChannel) app.thresholds(app.currentChannel)],'color',[0.7 0.3 0.2],'linewidth',2);
            app.previewSpikes();
            set(app.axis_snippet,'YLim',sort(app.thresholds(app.currentChannel)*[2 -2]));
        end
        
        function refreshProbe(app)
            cla(app.axis_probe);
            if ~isempty(app.probe)
                app.label_selectedProbe.Text = num2str(app.currentProbeNum);
                probeTmp = app.probe(app.currentProbeNum);
                plot(app.axis_probe,probeTmp.x,probeTmp.z,'sqr', 'MarkerSize',11);
                hold(app.axis_probe,"on");
                for i=1:probeTmp.nChannels
                    text(probeTmp.x(i)-5,probeTmp.z(i),...
                        num2str(probeTmp.channels(i)+1),'FontSize',9,'parent',app.axis_probe);
                end
                axis(app.axis_probe,[min(probeTmp.x)-50 max(probeTmp.x)+50 min(probeTmp.z)-50 max(probeTmp.z)+50]);
                axis(app.axis_probe,'equal');
                set(app.axis_probe,'FontSize',10,'TickDir','out');
                app.highlightChannel();
            end
        end
        
        function highlightChannel(app)
            if ~isempty(app.channelMarkerHandle)
                delete(app.channelMarkerHandle);
            end
            probeTmp = app.probe(app.currentProbeNum);
            cc = app.currentChannel;
            if app.currentProbeNum == 2 && (cc > app.probe(1).nChannels)
                cc = cc - app.probe(1).nChannels;
                app.channelMarkerHandle = plot(app.axis_probe,probeTmp.x(cc),probeTmp.z(cc),'bo', 'MarkerSize',11,'LineWidth',2);    
            elseif app.currentProbeNum == 1 && (cc <= app.probe(1).nChannels)
                app.channelMarkerHandle = plot(app.axis_probe,probeTmp.x(cc),probeTmp.z(cc),'bo', 'MarkerSize',11,'LineWidth',2);    
            end
        end
        
        function previewSpikes(app)
            sig = app.currentTc;
            th = app.thresholds(app.currentChannel);
            if th<0
                crossings = find(sig<th);
            else
                crossings = find(sig>th);
            end
            crossings(crossings<61) = [];
            crossings(crossings+60>length(sig)) = [];
            tt = linspace(-60*1000/app.sampleFreq,60*1000/app.sampleFreq,121);
            
            cla(app.axis_spikePreview);
            if ~isempty(crossings)
                if length(crossings) > 100
                    crossings = crossings(randperm(length(crossings),100));
                end
                clip = nan(length(crossings),121);
                for ii=1:length(crossings)
                    clip(ii,:) = sig(crossings(ii)-60:crossings(ii)+60);
                end
                plot(app.axis_spikePreview,tt,clip','k');
                grid(app.axis_spikePreview,'on');
            end
        end
        
        function [sample_freq,v1,v2] = getSampleFreq(app)
            filename = [app.settings.expFolder '/' app.exptName '/' app.exptName '_info.rhd'];
            fid2 = fopen(filename, 'r');
            
            % Check 'magic number' at beginning of file to make sure this is an Intan
            % Technologies RHD2000 data file.
            magic_number = fread(fid2, 1, 'uint32');
            if magic_number ~= hex2dec('c6912702')
                error('Unrecognized file type.');
            end
            
            % Read version number.
            v1 = fread(fid2, 1, 'int16');
            v2 = fread(fid2, 1, 'int16');
            
            % Read information of sampling rate and amplifier frequency settings.
            sample_freq = fread(fid2, 1, 'single');
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            addpath('microprobe_wiring/')
            addpath('microprobe_wiring/headstage and pcb wiring/')
            app.loadSettings();
            app.refreshFileInfo();
            app.refreshRemoteInfo();
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.ThresholdGUIUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {433, 433};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {282, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Button pushed function: button_getExpt
        function button_getExptPushed(app, event)
            f = uigetdir(app.settings.expFolder,'Please select the experiment folder');
            if isunix
                app.exptName = f(find(f == '/',1,'last')+1 : end);
                expFolder = f(1:find(f == '/',1,'last')-1);
            else
                app.exptName = f(find(f == '\',1,'last')+1 : end);
                expFolder = f(1:find(f == '\',1,'last')-1);
            end
            save('Settings.mat','expFolder','-append');
            
            app.settings.expFolder = expFolder;
            app.table_fileInfo.Data{1,2}{1} = app.exptName;
            app.table_fileInfo.Data{2,2}{1} = expFolder;
                        
            app.sampleFreq = app.getSampleFreq();
            app.table_fileInfo.Data{9,2}{1} = num2str(app.sampleFreq);
            
            idFileExists = exist([app.settings.expFolder '/' app.exptName '/' app.exptName '_id.mat'],'file');
            if idFileExists
                load([app.settings.expFolder '/' app.exptName '/' app.exptName '_id.mat'],'id')

                if ~strcmp(app.settings.name,id.name)
                    ButtonName = questdlg(['This file was processed by ' id.name '. Change to your name?'],...
                        'Name Mismatch', 'Yes','No','Yes');
                    if strcmp(ButtonName,'No')
                        app.settings.name = id.name;
                        app.settings.email = id.email;
                        
                        app.table_fileInfo.Data{5,2}{1} = app.settings.name;
                        app.table_fileInfo.Data{6,2}{1} = app.settings.email;
                    end
                end
                
                app.probe = id.probes;
                app.dropdown_probe1.Value = app.probe(1).type;
                app.dropdown_area1.Value = app.probe(1).area;
                if length(app.probe) == 2
                    app.dropdown_probe2.Value = app.probe(2).type;
                    app.dropdown_area2.Value = app.probe(2).area;
                    app.dropdown_area2.Enable = true;
                    app.table_fileInfo.Data{4,2}{1} = [app.probe(1).type ', ' app.probe(2).type];
                    app.nChannels = app.probe(1).nChannels + app.probe(2).nChannels;
                else
                    app.dropdown_probe2.Value = 'none';
                    app.table_fileInfo.Data{4,2}{1} = app.probe(1).type;
                    app.nChannels = app.probe(1).nChannels;
                end
            end
        end

        % Value changed function: dropdown_probe1
        function dropdown_probe1ValueChanged(app, event)
            type = app.dropdown_probe1.Value;
            if strcmp(type,'single')
                app.probe(1).channels = 0;
                app.probe(1).x = 0;
                app.probe(1).y = 0;
                app.probe(1).z = 0;
                app.probe(1).shaft = 1;
                app.probe(1).tipelectrode = 1;
                app.probe(1).type = type;
                app.probe(1).wiring = [0 0 0 0 1];
                app.probe(1).nShanks = 1;
                app.probe(1).nChannels = 1;
                app.probe(1).config = 1;
            elseif strcmp(type,'tetrode')
                app.probe(1).channels = [0 1 2 3]';
                app.probe(1).x = 20*[-1 -1 1 1]';
                app.probe(1).y = [0 0 0 0]';
                app.probe(1).z = 20*[-1 1 -1 1]';
                app.probe(1).shaft = [1 1 1 1]';
                app.probe(1).tipelectrode = 1;
                app.probe(1).type = type;
                app.probe(1).wiring = [app.probe(1).channels app.probe(1).x app.probe(1).y app.probe(1).z app.probe(1).shaft];
                app.probe(1).nShanks = 1;
                app.probe(1).nChannels = 4;
                app.probe(1).config = [1 2 3 4];
            elseif strcmp(type,'select')
                app.probe(1) = [];
            else
                eval(['probe_' type]); close;
            
                app.probe(1).channels = s.channels;
                app.probe(1).x = s.x;
                app.probe(1).y = s.y;
                app.probe(1).z = s.z;
                app.probe(1).shaft = s.shaft;
                app.probe(1).tipelectrode = s.tipelectrode;
                app.probe(1).type = type;
                app.probe(1).wiring = probewiring;
                app.probe(1).nShanks = length(unique(app.probe(1).shaft));
                app.probe(1).nChannels = numel(app.probe(1).channels);
                
                a = app.probe(1).wiring(app.probe(1).wiring(:,5) == 1,[1 4]);
                [~,idx] = sort(a(:,2),'descend');
                app.probe(1).config = a(idx,1);
                
                a = app.probe(1).wiring(app.probe(1).wiring(:,5) == 2,[1 4]);
                [~,idx] = sort(a(:,2),'descend');
                app.probe(1).config = [app.probe(1).config; a(idx,1)]';
            end
            
            if ~isempty(app.probe)
                app.probe(1).area = app.dropdown_area1.Value;
                if length(app.probe) == 1
                    app.nChannels = app.probe(1).nChannels;
                    app.table_fileInfo.Data{4,2}{1} = app.probe(1).type;
                else
                    app.nChannels = app.probe(1).nChannels + app.probe(2).nChannels;
                    app.table_fileInfo.Data{4,2}{1} = [app.probe(1).type ', ' app.probe(2).type];
                end
                app.badChannels = zeros(1,app.nChannels);
                app.thresholds = ones(1,app.nChannels).*(-200);
                app.currentChannel = 1;
                app.currentProbeNum = 1;
                
                app.refreshProbe();
            else
                cla(app.axis_probe);
            end
        end

        % Value changed function: dropdown_probe2
        function dropdown_probe2ValueChanged(app, event)
            type = app.dropdown_probe2.Value;
            
            if strcmp(type,'none')
                app.currentProbeNum = 1;
                app.currentChannel = 1;
                app.nChannels = app.probe(1).nChannels;
                app.table_fileInfo.Data{4,2}{1} = app.probe(1).type;
                app.dropdown_area2.Enable = false;
                app.probe(2) = [];
            else
                app.probe(2).area = app.dropdown_area2.Value;
                if strcmp(type,'single')
                    app.probe(2).channels = 0;
                    app.probe(2).x = 0;
                    app.probe(2).y = 0;
                    app.probe(2).z = 0;
                    app.probe(2).shaft = 1;
                    app.probe(2).tipelectrode = 1;
                    app.probe(2).type = type;
                    app.probe(2).wiring = [0 0 0 0 1];
                    app.probe(2).nShanks = 1;
                    app.probe(2).nChannels = 1;
                    app.probe(2).config = 1;
                elseif strcmp(type,'tetrode')
                    app.probe(2).channels = [0 1 2 3]';
                    app.probe(2).x = 20*[-1 -1 1 1]';
                    app.probe(2).y = [0 0 0 0]';
                    app.probe(2).z = 20*[-1 1 -1 1]';
                    app.probe(2).shaft = [1 1 1 1]';
                    app.probe(2).tipelectrode = 1;
                    app.probe(2).type = type;
                    app.probe(2).wiring = [app.probe(1).channels app.probe(1).x app.probe(1).y app.probe(1).z app.probe(1).shaft];
                    app.probe(2).nShanks = 1;
                    app.probe(2).nChannels = 4;
                    app.probe(2).config = [1 2 3 4];
                else
                    eval(['probe_' type]); close;
                    
                    app.probe(2).channels = s.channels;
                    app.probe(2).x = s.x;
                    app.probe(2).y = s.y;
                    app.probe(2).z = s.z;
                    app.probe(2).shaft = s.shaft;
                    app.probe(2).tipelectrode = s.tipelectrode;
                    app.probe(2).type = type;
                    app.probe(2).wiring = probewiring;
                    app.probe(2).nShanks = length(unique(app.probe(2).shaft));
                    app.probe(2).nChannels = numel(app.probe(2).channels);
                    
                    a = app.probe(2).wiring(app.probe(2).wiring(:,5) == 1,[1 4]);
                    [~,idx] = sort(a(:,2),'descend');
                    app.probe(2).config = a(idx,1);
                    
                    a = app.probe(2).wiring(app.probe(2).wiring(:,5) == 2,[1 4]);
                    [~,idx] = sort(a(:,2),'descend');
                    app.probe(2).config = [app.probe(2).config; a(idx,1)]';
                end
                
                app.nChannels = app.probe(1).nChannels + app.probe(2).nChannels;
                app.currentChannel = app.probe(1).nChannels + 1;
                app.table_fileInfo.Data{4,2}{1} = [app.probe(1).type ', ' app.probe(2).type];
                app.currentProbeNum = 2;
                
                app.dropdown_area2.Enable = true;
            end
            
            app.badChannels = zeros(1,app.nChannels);
            app.thresholds = ones(1,app.nChannels).*(-200);

            app.refreshProbe();
        end

        % Value changed function: dropdown_area1
        function dropdown_area1ValueChanged(app, event)
            app.probe(1).area = app.dropdown_area1.Value;
        end

        % Value changed function: dropdown_area2
        function dropdown_area2ValueChanged(app, event)
            app.probe(2).area = app.dropdown_area2.Value;
        end

        % Button pushed function: button_nextProbe
        function button_nextProbeButtonPushed(app, event)
            if length(app.probe) == 1
                app.currentProbeNum = 1;
            else
                app.currentProbeNum = 3-app.currentProbeNum;
            end
            app.label_selectedProbe.Text = num2str(app.currentProbeNum);
            app.refreshProbe();
        end

        % Button pushed function: button_prevProbe
        function button_prevProbeButtonPushed(app, event)
            if length(app.probe) == 1
                app.currentProbeNum = 1;
            else
                app.currentProbeNum = 3-app.currentProbeNum;
            end
            app.label_selectedProbe.Text = num2str(app.currentProbeNum);
            app.refreshProbe();
        end

        % Button pushed function: button_load
        function button_loadButtonPushed(app, event)
            if isempty(app.probe) || isempty(app.exptName)
                errordlg('Please select experiment and probe','Incomplete','non-modal')
                return;
            end
            
            fileinfo = dir([app.settings.expFolder '/' app.exptName '/' app.exptName '_amplifier.dat']);
            app.nTotalSamples = fileinfo.bytes/(2*app.nChannels);
            
            oldExptFileExists = exist([app.settings.expFolder '/' app.exptName '/experiment.mat'],'file');
            idFileExists = exist([app.settings.expFolder '/' app.exptName '/' app.exptName '_id.mat'],'file');
            threshFileExists = exist([app.settings.expFolder '/' app.exptName '/' app.exptName '_threshold.mat'],'file');
            
            if oldExptFileExists && ~idFileExists && ~threshFileExists
                load([app.settings.expFolder '/' app.exptName '/experiment.mat'],'b1','a1','BadCh','Th','probes')
                app.butterParams.b1 = b1;
                app.butterParams.a1 = a1;
                app.badChannels = BadCh;
                app.thresholds(BadCh == 0) = Th;
                app.table_fileInfo.Data{3,2}{1} = num2str(sum(app.badChannels));
                
                app.dropdown_probe1.Value = probes{1};
                if length(probes) == 2
                    app.dropdown_probe2.Value = probes{2};
                    app.table_fileInfo.Data{4,2}{1} = [probes{1} ', ' probes{2}];
                else
                    app.dropdown_probe2.Value = 'none';
                    app.table_fileInfo.Data{4,2}{1} = probes{1};
                end
            end
            
            if threshFileExists
                load([app.settings.expFolder '/' app.exptName '/' app.exptName '_threshold.mat'],'thresholding')
                app.thresholds = thresholding.thresholds;
                app.badChannels = thresholding.badChannels;
                app.butterParams = thresholding.butter;
                app.table_fileInfo.Data{3,2}{1} = num2str(sum(app.badChannels));
            end
            
            if ~isfield(app,'butterParams')
                [app.butterParams.b1, app.butterParams.a1] = butter(3, [app.settings.defaultFilter.highPass/app.sampleFreq,app.settings.defaultFilter.lowPass/app.sampleFreq]*2, 'bandpass');
                app.butterParams.highPass = app.settings.defaultFilter.highPass;
                app.butterParams.lowPass = app.settings.defaultFilter.lowPass;
            end
            
            app.currentChannel = 1;
            app.currentProbeNum = 1;
            app.refreshProbe();
            
            app.refreshSnippet(true,'random');
            app.slider_threshold.Value = app.thresholds(1);
        end

        % Button pushed function: button_save
        function button_saveButtonPushed(app, event)
            if ~exist([app.settings.expFolder '/' app.exptName '/SpikeFiles'],'dir')
                mkdir([app.settings.expFolder '/' app.exptName '/SpikeFiles'])
            end
            
            % for backward compatibility, for now
            experiment = app.exptName;
            Th = app.thresholds(~app.badChannels);
            BadCh = app.badChannels;
            b1 = app.butterParams.b1;
            a1 = app.butterParams.a1;
            isBR = 0;
            probes = {app.probe.type};
            if length(app.probe) == 1
                CHs = app.probe(1).config;
            else
                CHs = [app.probe(1).config app.probe(2).config];
            end
            save([app.settings.expFolder '/' app.exptName '/experiment.mat'],'experiment','Th','CHs','BadCh','b1','a1','isBR','probes');
            
            % the new way
            thresholding.thresholds = app.thresholds;
            thresholding.badChannels = app.badChannels;
            thresholding.butter = app.butterParams;
            
            id.exptId = app.exptName;
            id.probes = app.probe;
            id.name = app.settings.name;
            id.email = app.settings.email;
            id.isBR = false;
            id.sampleFreq = app.sampleFreq;
            id.filter = []
            
            save([app.settings.expFolder '/' app.exptName '/' app.exptName '_threshold.mat'],'thresholding');
            save([app.settings.expFolder '/' app.exptName '/' app.exptName '_id.mat'],'id');
            
            msgbox('Thresholds and id saved.','Saved','non-modal')
        end

        % Button pushed function: button_lastSnippet
        function button_lastSnippetButtonPushed(app, event)
            app.refreshSnippet(true,'last');
        end

        % Value changed function: slider_threshold
        function slider_thresholdChanged(app, event)
            app.thresholds(app.currentChannel) = event.Value;
            app.plotThreshold();
        end

        % Button pushed function: button_nextChan
        function button_nextChanButtonPushed(app, event)
            if app.currentChannel < app.nChannels
                app.currentChannel = app.currentChannel + 1;
                app.text_channelNum.Value = app.currentChannel;
                app.slider_threshold.Value = app.thresholds(app.currentChannel);
                
                if app.currentChannel > app.probe(1).nChannels && app.currentProbeNum == 1
                    app.currentProbeNum = 2;
                    app.refreshProbe();
                end
                app.refreshSnippet(false);
            end
        end

        % Button pushed function: button_prevChan
        function button_prevChanButtonPushed(app, event)
            if app.currentChannel > 1
                app.currentChannel = app.currentChannel - 1;
                app.text_channelNum.Value = app.currentChannel;
                app.slider_threshold.Value = app.thresholds(app.currentChannel);
                if app.currentChannel <= app.probe(1).nChannels && app.currentProbeNum == 2
                    app.currentProbeNum = 1;
                    app.refreshProbe();
                end
                app.refreshSnippet(false);
            end
        end

        % Value changed function: text_channelNum
        function text_channelNumValueChanged(app, event)
            app.currentChannel = app.text_channelNum.Value;
            if app.currentChannel > app.probe(1).nChannels && length(app.probe) > 1
                app.currentProbeNum = 2;
            elseif app.currentChannel <= app.probe(1).nChannels
                app.currentProbeNum = 1;
            else
                app.currentChannel = event.PreviousValue;
                app.text_channelNum.Value = event.PreviousValue;
                app.currentProbeNum = 1;
            end
            app.refreshProbe();
            app.refreshSnippet(false);
        end

        % Button pushed function: button_firstSnippet
        function button_firstSnippetButtonPushed(app, event)
            app.refreshSnippet(true,'first');
        end

        % Button pushed function: button_randomSnippet
        function button_randomSnippetPushed(app, event)
            app.refreshSnippet(true,'random');
        end

        % Button pushed function: button_rejectChan
        function button_rejectChanPushed(app, event)
            if app.badChannels(app.currentChannel)
                app.badChannels(app.currentChannel) = 0;
                app.button_rejectChan.Text = 'Reject Channel';
            else
                app.badChannels(app.currentChannel) = 1;
                app.button_rejectChan.Text = 'Accept Channel';
            end
            app.table_fileInfo.Data{3,2}{1} = num2str(sum(app.badChannels));
            app.refreshSnippet(false)
        end

        % Selection changed function: radio_compSelect
        function radio_compSelectSelectionChanged(app, event)
            selectedButton = app.radio_compSelect.SelectedObject;
            
            if strcmp(selectedButton.Text,'remote')
                app.settings.doRemote = true;     
            else
                app.settings.doRemote = false;
            end
            
            app.refreshRemoteInfo();
        end

        % Value changed function: text_ip
        function text_ipValueChanged(app, event)
            app.settings.remoteIP = app.text_ip.Value;
            app.refreshRemoteInfo();
        end

        % Button pushed function: button_checkConn
        function button_checkConnButtonPushed(app, event)
            % [] = system(['ssh nielsenlab@' app.settings.remoteIP]);
        end

        % Button pushed function: button_copyData
        function button_copyDataButtonPushed(app, event)
            
        end

        % Button pushed function: button_extract
        function button_extractButtonPushed(app, event)
            if app.settings.doRemote
                
            else
                
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ThresholdGUIUIFigure and hide until all components are created
            app.ThresholdGUIUIFigure = uifigure('Visible', 'off');
            app.ThresholdGUIUIFigure.AutoResizeChildren = 'off';
            app.ThresholdGUIUIFigure.Position = [100 100 1358 433];
            app.ThresholdGUIUIFigure.Name = 'Threshold GUI';
            app.ThresholdGUIUIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.ThresholdGUIUIFigure);
            app.GridLayout.ColumnWidth = {282, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Scrollable = 'on';

            % Create panel_location
            app.panel_location = uipanel(app.LeftPanel);
            app.panel_location.AutoResizeChildren = 'off';
            app.panel_location.Title = 'File';
            app.panel_location.Position = [8 265 269 151];

            % Create button_getExpt
            app.button_getExpt = uibutton(app.panel_location, 'push');
            app.button_getExpt.ButtonPushedFcn = createCallbackFcn(app, @button_getExptPushed, true);
            app.button_getExpt.Position = [10 98 100 22];
            app.button_getExpt.Text = 'Get Experiment';

            % Create button_save
            app.button_save = uibutton(app.panel_location, 'push');
            app.button_save.ButtonPushedFcn = createCallbackFcn(app, @button_saveButtonPushed, true);
            app.button_save.Position = [144 8 116 22];
            app.button_save.Text = 'Save';

            % Create button_load
            app.button_load = uibutton(app.panel_location, 'push');
            app.button_load.ButtonPushedFcn = createCallbackFcn(app, @button_loadButtonPushed, true);
            app.button_load.Position = [10 8 116 22];
            app.button_load.Text = 'Load';

            % Create ProbeAreaLabel_2
            app.ProbeAreaLabel_2 = uilabel(app.panel_location);
            app.ProbeAreaLabel_2.HorizontalAlignment = 'right';
            app.ProbeAreaLabel_2.Position = [8 39 66 22];
            app.ProbeAreaLabel_2.Text = 'Probe/Area';

            % Create dropdown_probe2
            app.dropdown_probe2 = uidropdown(app.panel_location);
            app.dropdown_probe2.Items = {'none', '64D', '64E', '64F', '64G', '64H', '64M', '128P', 'single', 'tetrode', ''};
            app.dropdown_probe2.ValueChangedFcn = createCallbackFcn(app, @dropdown_probe2ValueChanged, true);
            app.dropdown_probe2.Position = [89 39 100 22];
            app.dropdown_probe2.Value = 'none';

            % Create ProbeAreaLabel
            app.ProbeAreaLabel = uilabel(app.panel_location);
            app.ProbeAreaLabel.HorizontalAlignment = 'right';
            app.ProbeAreaLabel.Position = [8 69 66 22];
            app.ProbeAreaLabel.Text = 'Probe/Area';

            % Create dropdown_probe1
            app.dropdown_probe1 = uidropdown(app.panel_location);
            app.dropdown_probe1.Items = {'select', '64D', '64E', '64F', '64G', '64H', '64M', '128P', 'single', 'tetrode'};
            app.dropdown_probe1.ValueChangedFcn = createCallbackFcn(app, @dropdown_probe1ValueChanged, true);
            app.dropdown_probe1.Position = [89 69 100 22];
            app.dropdown_probe1.Value = 'select';

            % Create dropdown_area1
            app.dropdown_area1 = uidropdown(app.panel_location);
            app.dropdown_area1.Items = {'select', 'V1', 'PSS', 'LGN', 'S1'};
            app.dropdown_area1.ValueChangedFcn = createCallbackFcn(app, @dropdown_area1ValueChanged, true);
            app.dropdown_area1.Position = [197 69 63 22];
            app.dropdown_area1.Value = 'select';

            % Create dropdown_area2
            app.dropdown_area2 = uidropdown(app.panel_location);
            app.dropdown_area2.Items = {'none', 'V1', 'PSS', 'LGN', 'S1'};
            app.dropdown_area2.ValueChangedFcn = createCallbackFcn(app, @dropdown_area2ValueChanged, true);
            app.dropdown_area2.Enable = 'off';
            app.dropdown_area2.Position = [197 39 63 22];
            app.dropdown_area2.Value = 'none';

            % Create panel_fileInfo
            app.panel_fileInfo = uipanel(app.LeftPanel);
            app.panel_fileInfo.AutoResizeChildren = 'off';
            app.panel_fileInfo.Title = 'File Details';
            app.panel_fileInfo.Position = [8 12 269 244];

            % Create table_fileInfo
            app.table_fileInfo = uitable(app.panel_fileInfo);
            app.table_fileInfo.ColumnName = {'Name'; 'Prop'};
            app.table_fileInfo.RowName = {''; ''};
            app.table_fileInfo.Position = [10 6 250 213];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Scrollable = 'on';

            % Create panel_signal
            app.panel_signal = uipanel(app.RightPanel);
            app.panel_signal.Position = [294 12 770 404];

            % Create axis_snippet
            app.axis_snippet = uiaxes(app.panel_signal);
            xlabel(app.axis_snippet, 'Time (ms)')
            ylabel(app.axis_snippet, 'Signal (mV)')
            app.axis_snippet.PlotBoxAspectRatio = [1.31331168831169 1 1];
            app.axis_snippet.GridAlpha = 0.15;
            app.axis_snippet.MinorGridAlpha = 0.25;
            app.axis_snippet.Box = 'on';
            app.axis_snippet.Position = [8 9 455 351];

            % Create ThresholdSliderLabel
            app.ThresholdSliderLabel = uilabel(app.panel_signal);
            app.ThresholdSliderLabel.HorizontalAlignment = 'center';
            app.ThresholdSliderLabel.VerticalAlignment = 'top';
            app.ThresholdSliderLabel.Position = [474 359 76 22];
            app.ThresholdSliderLabel.Text = 'Threshold';

            % Create slider_threshold
            app.slider_threshold = uislider(app.panel_signal);
            app.slider_threshold.Limits = [-500 500];
            app.slider_threshold.MajorTicks = [-500 -400 -300 -200 -100 0 100 200 300 400 500];
            app.slider_threshold.Orientation = 'vertical';
            app.slider_threshold.ValueChangedFcn = createCallbackFcn(app, @slider_thresholdChanged, true);
            app.slider_threshold.MinorTicks = [-500 -490 -480 -470 -460 -450 -440 -430 -420 -410 -400 -390 -380 -370 -360 -350 -340 -330 -320 -310 -300 -290 -280 -270 -260 -250 -240 -230 -220 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500];
            app.slider_threshold.Position = [479 44 3 308];
            app.slider_threshold.Value = -200;

            % Create button_randomSnippet
            app.button_randomSnippet = uibutton(app.panel_signal, 'push');
            app.button_randomSnippet.ButtonPushedFcn = createCallbackFcn(app, @button_randomSnippetPushed, true);
            app.button_randomSnippet.Position = [77 369 61 22];
            app.button_randomSnippet.Text = 'Random';

            % Create button_nextChan
            app.button_nextChan = uibutton(app.panel_signal, 'push');
            app.button_nextChan.ButtonPushedFcn = createCallbackFcn(app, @button_nextChanButtonPushed, true);
            app.button_nextChan.Position = [297 369 23 23];
            app.button_nextChan.Text = '>';

            % Create button_prevChan
            app.button_prevChan = uibutton(app.panel_signal, 'push');
            app.button_prevChan.ButtonPushedFcn = createCallbackFcn(app, @button_prevChanButtonPushed, true);
            app.button_prevChan.Position = [229 369 25 23];
            app.button_prevChan.Text = '<';

            % Create button_rejectChan
            app.button_rejectChan = uibutton(app.panel_signal, 'push');
            app.button_rejectChan.ButtonPushedFcn = createCallbackFcn(app, @button_rejectChanPushed, true);
            app.button_rejectChan.Position = [330 369 133 22];
            app.button_rejectChan.Text = 'Reject Channel';

            % Create axis_spikePreview
            app.axis_spikePreview = uiaxes(app.panel_signal);
            xlabel(app.axis_spikePreview, 'Time (ms)')
            ylabel(app.axis_spikePreview, 'Signal (mV)')
            app.axis_spikePreview.PlotBoxAspectRatio = [1.31331168831169 1 1];
            app.axis_spikePreview.GridAlpha = 0.15;
            app.axis_spikePreview.MinorGridAlpha = 0.25;
            app.axis_spikePreview.Box = 'on';
            app.axis_spikePreview.Position = [548 183 219 188];

            % Create panel_extract
            app.panel_extract = uipanel(app.panel_signal);
            app.panel_extract.Tooltip = {'Work in Progress'};
            app.panel_extract.Position = [549 9 210 173];

            % Create radio_compSelect
            app.radio_compSelect = uibuttongroup(app.panel_extract);
            app.radio_compSelect.SelectionChangedFcn = createCallbackFcn(app, @radio_compSelectSelectionChanged, true);
            app.radio_compSelect.Position = [9 112 100 52];

            % Create radio_local
            app.radio_local = uiradiobutton(app.radio_compSelect);
            app.radio_local.Enable = 'off';
            app.radio_local.Text = 'local';
            app.radio_local.Position = [11 25 47 22];
            app.radio_local.Value = true;

            % Create radio_remote
            app.radio_remote = uiradiobutton(app.radio_compSelect);
            app.radio_remote.Enable = 'off';
            app.radio_remote.Text = 'remote';
            app.radio_remote.Position = [11 3 60 22];

            % Create remoteipEditFieldLabel
            app.remoteipEditFieldLabel = uilabel(app.panel_extract);
            app.remoteipEditFieldLabel.HorizontalAlignment = 'right';
            app.remoteipEditFieldLabel.Enable = 'off';
            app.remoteipEditFieldLabel.Position = [8 81 56 22];
            app.remoteipEditFieldLabel.Text = 'remote ip';

            % Create text_ip
            app.text_ip = uieditfield(app.panel_extract, 'text');
            app.text_ip.ValueChangedFcn = createCallbackFcn(app, @text_ipValueChanged, true);
            app.text_ip.Enable = 'off';
            app.text_ip.Position = [79 81 120 22];

            % Create button_checkConn
            app.button_checkConn = uibutton(app.panel_extract, 'push');
            app.button_checkConn.ButtonPushedFcn = createCallbackFcn(app, @button_checkConnButtonPushed, true);
            app.button_checkConn.Enable = 'off';
            app.button_checkConn.Position = [9 52 88 23];
            app.button_checkConn.Text = 'check conn';

            % Create button_copyData
            app.button_copyData = uibutton(app.panel_extract, 'push');
            app.button_copyData.ButtonPushedFcn = createCallbackFcn(app, @button_copyDataButtonPushed, true);
            app.button_copyData.Enable = 'off';
            app.button_copyData.Position = [111 52 88 23];
            app.button_copyData.Text = 'copy data';

            % Create button_extract
            app.button_extract = uibutton(app.panel_extract, 'push');
            app.button_extract.ButtonPushedFcn = createCallbackFcn(app, @button_extractButtonPushed, true);
            app.button_extract.Enable = 'off';
            app.button_extract.Position = [111 18 88 23];
            app.button_extract.Text = 'extract';

            % Create button_firstSnippet
            app.button_firstSnippet = uibutton(app.panel_signal, 'push');
            app.button_firstSnippet.ButtonPushedFcn = createCallbackFcn(app, @button_firstSnippetButtonPushed, true);
            app.button_firstSnippet.Position = [8 369 61 22];
            app.button_firstSnippet.Text = 'First';

            % Create button_lastSnippet
            app.button_lastSnippet = uibutton(app.panel_signal, 'push');
            app.button_lastSnippet.ButtonPushedFcn = createCallbackFcn(app, @button_lastSnippetButtonPushed, true);
            app.button_lastSnippet.Position = [147 369 61 22];
            app.button_lastSnippet.Text = 'Last';

            % Create text_channelNum
            app.text_channelNum = uieditfield(app.panel_signal, 'numeric');
            app.text_channelNum.ValueChangedFcn = createCallbackFcn(app, @text_channelNumValueChanged, true);
            app.text_channelNum.HorizontalAlignment = 'center';
            app.text_channelNum.Position = [257 369 37 22];
            app.text_channelNum.Value = 1;

            % Create panel_probe
            app.panel_probe = uipanel(app.RightPanel);
            app.panel_probe.AutoResizeChildren = 'off';
            app.panel_probe.Position = [8 12 269 404];

            % Create axis_probe
            app.axis_probe = uiaxes(app.panel_probe);
            title(app.axis_probe, '')
            xlabel(app.axis_probe, 'X')
            ylabel(app.axis_probe, 'Y')
            app.axis_probe.PlotBoxAspectRatio = [1 1.53 1];
            app.axis_probe.Position = [10 9 249 350];

            % Create label_selectedProbe
            app.label_selectedProbe = uilabel(app.panel_probe);
            app.label_selectedProbe.HorizontalAlignment = 'center';
            app.label_selectedProbe.Position = [121 368 27 24];
            app.label_selectedProbe.Text = '1';

            % Create button_nextProbe
            app.button_nextProbe = uibutton(app.panel_probe, 'push');
            app.button_nextProbe.ButtonPushedFcn = createCallbackFcn(app, @button_nextProbeButtonPushed, true);
            app.button_nextProbe.Position = [156 369 23 23];
            app.button_nextProbe.Text = '>';

            % Create button_prevProbe
            app.button_prevProbe = uibutton(app.panel_probe, 'push');
            app.button_prevProbe.ButtonPushedFcn = createCallbackFcn(app, @button_prevProbeButtonPushed, true);
            app.button_prevProbe.Position = [89 369 25 23];
            app.button_prevProbe.Text = '<';

            % Show the figure after all components are created
            app.ThresholdGUIUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = thresholdGui_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ThresholdGUIUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ThresholdGUIUIFigure)
        end
    end
end