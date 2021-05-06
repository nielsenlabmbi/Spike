classdef thresholdGui_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ThresholdGUIUIFigure  matlab.ui.Figure
        GridLayout            matlab.ui.container.GridLayout
        LeftPanel             matlab.ui.container.Panel
        panel_location        matlab.ui.container.Panel
        button_getExpt        matlab.ui.control.Button
        button_save           matlab.ui.control.Button
        ProbeAreaLabel_2      matlab.ui.control.Label
        dropdown_probe2       matlab.ui.control.DropDown
        ProbeAreaLabel        matlab.ui.control.Label
        dropdown_probe1       matlab.ui.control.DropDown
        dropdown_area1        matlab.ui.control.DropDown
        dropdown_area2        matlab.ui.control.DropDown
        SettingsButton        matlab.ui.control.Button
        panel_fileInfo        matlab.ui.container.Panel
        table_fileInfo        matlab.ui.control.Table
        RightPanel            matlab.ui.container.Panel
        panel_signal          matlab.ui.container.Panel
        axis_snippet          matlab.ui.control.UIAxes
        ThresholdSliderLabel  matlab.ui.control.Label
        slider_threshold      matlab.ui.control.Slider
        button_nextChan       matlab.ui.control.Button
        button_prevChan       matlab.ui.control.Button
        button_rejectChan     matlab.ui.control.Button
        axis_spikePreview     matlab.ui.control.UIAxes
        text_channelNum       matlab.ui.control.NumericEditField
        ChannelLabel          matlab.ui.control.Label
        SegmentstartsLabel    matlab.ui.control.Label
        slider_time           matlab.ui.control.Slider
        text_minThreshold     matlab.ui.control.NumericEditField
        text_maxThreshold     matlab.ui.control.NumericEditField
        PresetthresholdPanel  matlab.ui.container.Panel
        RMSSpinnerLabel       matlab.ui.control.Label
        spinner_RMS           matlab.ui.control.Spinner
        button_applyRMS       matlab.ui.control.Button
        lamp_applyRMS         matlab.ui.control.Lamp
        axis_timecourse       matlab.ui.control.UIAxes
        GotosEditFieldLabel   matlab.ui.control.Label
        text_gotoTime         matlab.ui.control.NumericEditField
        check_channel         matlab.ui.control.CheckBox
        panel_probe           matlab.ui.container.Panel
        axis_probe            matlab.ui.control.UIAxes
        label_selectedProbe   matlab.ui.control.Label
        button_nextProbe      matlab.ui.control.Button
        button_prevProbe      matlab.ui.control.Button
        ProbeLabel            matlab.ui.control.Label
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
        nProbes
        thresholds
        badChannels
        chkChannels %has a channel been looked at
        nTotalSamples
        snippet %holds unfiltered waveform x channel
        butterParams
        snippetRMS %holds the RMS of the different channels
        timeCourse %holds long time course for current channel
        
        currentProbeNum
        currentChannel
        currentDisplayChannel %to make it easier for multiple probes
        currentTc %filtered waveform for current channel
        currentStart
        
        thresholdLineHandle
        signalLineHandle
        channelMarkerHandle
        timeLineHandle
    end
    
    methods (Access = private)
        
        function refreshFileInfo(app)
            % display Settings
            % display FileProps
            
            data = {
                'expt name' '';...
                'exp loc' app.settings.expFolder;...
                'n probes' 'XXX';...
                'n channels' 'XXX';...
                'sample freq' 'XXX';...
                'n bad chans' '0';...
                
                'high pass' num2str(app.settings.filterHighPass);...
                'low pass' num2str(app.settings.filterLowPass);...
                'nr snippets rms' num2str(app.settings.nrSnippetsRMS);...
                'use GPU' num2str(app.settings.useGPU);...
                'length timecourse' num2str(app.settings.lengthTC)
                
                'name' app.settings.name;...
                };
            
            app.table_fileInfo.Data = cell2table(data);
        end
        
        
        function refreshSnippet(app,getNew)
            if getNew
                app.getNewSnippet();
            end
            
            cla(app.axis_snippet);
            app.currentTc = filter(app.butterParams.b1, app.butterParams.a1, app.snippet(app.currentChannel,:));
            app.signalLineHandle = plot(app.axis_snippet,linspace(0,1,app.sampleFreq)+app.slider_time.Value,app.currentTc,'k');
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
        
        function getNewSnippet(app)
            dataFileName = [app.settings.expFolder '/' app.exptName '/' app.exptName '_amplifier.dat'];
            DataFile = fopen(dataFileName,'r');
            frewind(DataFile);
            fseek(DataFile,2*floor(app.slider_time.Value*app.sampleFreq)*app.nChannels,'bof');
            
            app.snippet = fread(DataFile, [app.nChannels app.sampleFreq], 'int16');
        end
        
        
        function plotThreshold(app)
            if ~isempty(app.thresholdLineHandle)
                if isvalid(app.thresholdLineHandle)
                    delete(app.thresholdLineHandle);
                end
            end
            app.thresholdLineHandle = line(app.axis_snippet,[0 1000]+app.slider_time.Value,...
                [app.thresholds(app.currentChannel) app.thresholds(app.currentChannel)],'color',[0.7 0.3 0.2],'linewidth',2);
            app.previewSpikes();
            set(app.axis_snippet,'XLim',[0 1]+app.slider_time.Value);
            set(app.axis_snippet,'YLim',sort(app.thresholds(app.currentChannel)*[2 -2]));
        end
        
        function refreshProbe(app)
            cla(app.axis_probe);
            if isfield(app.probe,'type')
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
        
        function refreshChannel(app)
            
            %check to make sure threshold is not out of range for
            %slider
            if app.thresholds(app.currentChannel)<app.slider_threshold.Limits(1)
                app.slider_threshold.Limits(1)=floor(app.thresholds(app.currentChannel)/100)*100;
                app.text_minThreshold.Value=app.slider_threshold.Limits(1);
            end
            if app.thresholds(app.currentChannel)>app.slider_threshold.Limits(2)
                app.slider_threshold.Limits(2)=ceil(app.thresholds(app.currentChannel)/100)*100;
                app.text_maxThreshold.Value=app.slider_threshold.Limits(2);
            end
            
            app.slider_threshold.Value = app.thresholds(app.currentChannel);
            app.check_channel.Value = app.chkChannels(app.currentChannel);
            app.text_channelNum.Value = app.currentChannel;
        end
        
        function highlightChannel(app)
            if isfield(app.probe,'type')
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
        end
        
        function previewSpikes(app)
            sig = app.currentTc;
            th = app.thresholds(app.currentChannel);
            
            %find crossings - we want the transitions
            if th<0
                crossings = find(sig(2:end)<th & sig(1:end-1)>=th);
            else
                crossings = find(sig(2:end)>th & sig(1:end-1)<=th);
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
        
        
        function getNewTimeCourse(app)
            %load only data for the current channel
            dataFileName = [app.settings.expFolder '/' app.exptName '/' app.exptName '_amplifier.dat'];
            DataFile = fopen(dataFileName,'r');
            frewind(DataFile);
            fseek(DataFile,2*floor(app.currentStart*app.sampleFreq)*app.nChannels+2*(app.currentChannel-1),'bof');
            
            tc = fread(DataFile, app.settings.lengthTC*app.sampleFreq, 'int16',2*(app.nChannels-1));
            app.timeCourse=filter(app.butterParams.b1, app.butterParams.a1, tc);
            clear tc;
        end
        
        function refreshTimeCourse(app)
            cla(app.axis_timecourse);
            minTc=max(-2000,min(app.timeCourse)); %to keep within reasonable limits
            maxTc=min(2000,max(app.timeCourse));
            
            plot(app.axis_timecourse,linspace(app.currentStart,app.currentStart+app.settings.lengthTC,app.settings.lengthTC*app.sampleFreq),app.timeCourse,'k');
            set(app.axis_timecourse,'YLim',[minTc maxTc],'XLim',[app.currentStart app.currentStart+app.settings.lengthTC]); %fix to something reasonable for now
            
            app.refreshTimeMarker();
        end
        
        function refreshTimeMarker(app)
            startTimeSnippet=app.slider_time.Value;
            minTc=max(-2000,min(app.timeCourse)); %to keep within reasonable limits
            maxTc=min(2000,max(app.timeCourse));
            
            if ~isempty(app.timeLineHandle)
                if isvalid(app.timeLineHandle)
                    delete(app.timeLineHandle);
                end
            end
            hold(app.axis_timecourse,'on')
            if startTimeSnippet>=app.currentStart && startTimeSnippet<=app.currentStart+app.settings.lengthTC
                app.timeLineHandle=plot(app.axis_timecourse,[startTimeSnippet startTimeSnippet],[minTc maxTc],'r-');
            end
            
        end
        
        
    end
    
    methods (Access = public)
        
        function getSettings(app,settings)
            app.settings=settings;
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            addpath('microprobe_wiring/')
            addpath('microprobe_wiring/headstage and pcb wiring/')
            
            %initialize settings
            app.settings.name='Nielsen Lab';
            app.settings.filterHighPass = 250; %Hz
            app.settings.filterLowPass = 5000; %Hz
            app.settings.expFolder = 'z:\ephys_new';
            app.settings.nrSnippetsRMS = 100;
            app.settings.useGPU=1;
            app.settings.lengthTC=20; %sec
            app.refreshFileInfo();
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.ThresholdGUIUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {578, 578};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {287, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Button pushed function: SettingsButton
        function button_settingsButtonPushed(app, event)
            %get settings
            settingsHandle=thresholdGuiSettings(app);
            waitfor(settingsHandle)
            
            %execute startup functions
            app.refreshFileInfo();
            %app.refreshRemoteInfo();
        end

        % Button pushed function: button_getExpt
        function button_getExptPushed(app, event)
            %get experiment folder and name
            f = uigetdir(app.settings.expFolder,'Please select the experiment folder');
            if f==0
                return;
            end
            figure(app.ThresholdGUIUIFigure)
            [expFolder,app.exptName]=fileparts(f);
            app.settings.expFolder = expFolder;
            app.table_fileInfo.Data{1,2}{1} = app.exptName;
            app.table_fileInfo.Data{2,2}{1} = expFolder;
            
            %read intan header
            header=read_Intan_Header(fullfile(app.settings.expFolder,app.exptName,[app.exptName '_info.rhd']));
            
            %get sampling frequency and update table
            app.sampleFreq = header.sample_rate;
            app.table_fileInfo.Data{5,2}{1} = num2str(app.sampleFreq);
            
            %check whether ID file exists
            idFileExists = exist(fullfile(app.settings.expFolder,app.exptName,[app.exptName '_id.mat']),'file');
            
            if ~idFileExists %file not opened before, make guesses based on intan header
                %if there are 2 probes, port A and B should be active
                %look for port A first - this will usually be the
                %headstage present if only one is used
                portIdx=find(strcmp(header.signal_group_name,'Port A')==1);
                if ~isempty(portIdx)
                    portEnable=header.signal_group_enabled(portIdx);
                    if portEnable
                        app.probe(1).nChannels=header.signal_group_num_amp_enabled(portIdx); %this just prevents errors of mismatch between data and loaded probe
                        app.nChannels=app.probe(1).nChannels;
                        app.nProbes=1;
                    end
                end
                
                %look for port B
                portIdx=find(strcmp(header.signal_group_name,'Port B')==1);
                if ~isempty(portIdx)
                    portEnable=header.signal_group_enabled(portIdx);
                    if portEnable
                        app.dropdown_area2.Enable = true;
                        app.probe(2).nChannels=header.signal_group_num_amp_enabled(portIdx); %this just prevents errors of mismatch between data and loaded probe
                        app.nChannels=app.nChannels+app.probe(2).nChannels;
                        app.nProbes=2;
                    end
                end
                
                %update table, dropdowns etc
                app.table_fileInfo.Data{3,2}{1} = num2str(app.nProbes);
                app.dropdown_area1.Enable=1;
                app.dropdown_probe1.Enable=1;
                if app.nProbes==2
                    app.dropdown_area2.Enable = true;
                    app.dropdown_probe2.Enable = 1;
                    app.table_fileInfo.Data{4,2}{1} = [num2str(app.probe(1).nChannels) ', ' num2str(app.probe(2).nChannels)];
                else
                    app.dropdown_probe2.Value = 'select';
                    app.dropdown_probe2.Enable=0;
                    app.dropdown_area2.Value='select';
                    app.dropdown_area2.Enable=0;
                    app.table_fileInfo.Data{4,2}{1} = num2str(app.probe(1).nChannels);
                end
                
            else  %file already opened before, so use the previous settings
                load(fullfile(app.settings.expFolder,app.exptName,[app.exptName '_id.mat']),'id')
                
                if ~strcmp(app.settings.name,id.name)
                    ButtonName = questdlg(['This file was processed by ' id.name '. Change to your name?'],...
                        'Name Mismatch', 'Yes','No','Yes');
                    if strcmp(ButtonName,'No')
                        app.settings.name = id.name;
                        
                        app.table_fileInfo.Data{12,2}{1} = app.settings.name;
                    end
                end
                
                %update probe
                app.probe = id.probes;
                app.nProbes = length(app.probe);
                app.table_fileInfo.Data{3,2}{1} = num2str(app.nProbes);
                
                app.dropdown_probe1.Value = app.probe(1).type;
                app.dropdown_area1.Value = app.probe(1).area;
                app.dropdown_probe1.Enable=1;
                app.dropdown_area1.Enable=1;
                
                if app.nProbes == 2
                    app.dropdown_probe2.Value = app.probe(2).type;
                    app.dropdown_area2.Value = app.probe(2).area;
                    app.dropdown_area2.Enable = true;
                    app.dropdown_probe2.Enable = 1;
                    app.nChannels = app.probe(1).nChannels + app.probe(2).nChannels;
                    app.table_fileInfo.Data{4,2}{1} = [num2str(app.probe(1).nChannels) ', ' num2str(app.probe(2).nChannels)];
                else
                    app.dropdown_probe2.Value = 'select';
                    app.dropdown_probe2.Enable=0;
                    app.dropdown_area2.Value='select';
                    app.dropdown_area2.Enable=0;
                    app.nChannels = app.probe(1).nChannels;
                    app.table_fileInfo.Data{4,2}{1} = num2str(app.probe(1).nChannels);
                end
            end
            
            %get duration data
            fileinfo = dir(fullfile(app.settings.expFolder,app.exptName,[app.exptName '_amplifier.dat']));
            app.nTotalSamples = fileinfo.bytes/(2*app.nChannels);
            maxTime=floor((app.nTotalSamples-app.sampleFreq)/app.sampleFreq);
            app.slider_time.Limits=[0 maxTime];
            app.slider_time.MajorTicks=[0:100:maxTime];
            app.slider_time.Value=0;
            app.text_gotoTime.Limits=[0 maxTime];
            app.currentStart=0; %start for long time course [in sec]
            
            %initialize arrays
            app.badChannels = zeros(1,app.nChannels);
            app.thresholds = ones(1,app.nChannels).*(-200);
            app.chkChannels = zeros(1,app.nChannels);
            app.snippetRMS = [];
            
            %check whether file has already been thresholded and load those
            %settings
            oldExptFileExists = exist(fullfile(app.settings.expFolder,app.exptName,'experiment.mat'),'file');
            threshFileExists = exist(fullfile(app.settings.expFolder,app.exptName,[app.exptName '_threshold.mat']),'file');
            
            if oldExptFileExists && ~threshFileExists %for backwards compatibility
                load(fullfile(app.settings.expFolder,app.exptName,'experiment.mat'),'b1','a1','BadCh','Th')
                app.butterParams.b1 = b1;
                app.butterParams.a1 = a1;
                app.badChannels = BadCh;
                app.thresholds(BadCh == 0) = Th;
                app.table_fileInfo.Data{6,2}{1} = num2str(sum(app.badChannels));
            end
            
            if threshFileExists
                load(fullfile(app.settings.expFolder,app.exptName,[app.exptName '_threshold.mat']),'thresholding')
                app.thresholds = thresholding.thresholds;
                app.badChannels = thresholding.badChannels;
                app.butterParams = thresholding.butter;
                app.chkChannels = thresholding.chkChannels;
                app.settings.filterHighPass=app.butterParams.highPass;
                app.settings.filterLowPass=app.butterParams.lowPass;
                
                app.table_fileInfo.Data{6,2}{1} = num2str(sum(app.badChannels));
                app.table_fileInfo.Data{7,2}{1} = num2str(app.settings.filterHighPass);
                app.table_fileInfo.Data{8,2}{1} = num2str(app.settings.filterLowPass);
            end
            
            if ~isfield(app,'butterParams')
                [app.butterParams.b1, app.butterParams.a1] = butter(3, [app.settings.filterHighPass/app.sampleFreq,app.settings.filterLowPass/app.sampleFreq]*2, 'bandpass');
                app.butterParams.highPass = app.settings.filterHighPass;
                app.butterParams.lowPass = app.settings.filterLowPass;
            end
            
            %plot probe (if possible)
            app.currentChannel = 1;
            app.currentDisplayChannel = 1;
            app.currentProbeNum = 1;
            app.refreshProbe();
            if app.nProbes==1
                app.button_nextProbe.Enable=0;
                app.button_prevProbe.Enable=0;
            else
                app.button_nextProbe.Enable=1;
                app.button_prevProbe.Enable=1;
            end
            
            %plot snippets
            app.refreshSnippet(true);
            app.getNewTimeCourse();
            app.refreshTimeCourse();
            
            %enable snippet selection etc
            app.slider_threshold.Value = app.thresholds(1);
            app.check_channel.Value = app.chkChannels(1);
            app.slider_threshold.Enable=1;
            app.text_minThreshold.Enable=1;
            app.text_maxThreshold.Enable=1;
            app.slider_time.Enable=1;
            app.button_nextChan.Enable=1;
            app.button_prevChan.Enable=1;
            app.button_rejectChan.Enable=1;
            app.text_channelNum.Enable=1;
            app.button_save.Enable=1;
            app.button_applyRMS.Enable=1;
            app.spinner_RMS.Enable=1;
            
            
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
                app.refreshProbe();
            else
                cla(app.axis_probe);
            end
        end

        % Value changed function: dropdown_probe2
        function dropdown_probe2ValueChanged(app, event)
            type = app.dropdown_probe2.Value;
            
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
                app.probe(2).config = [1 2 3 4];
            elseif strcmp(type,'select')
                app.probe(2) = [];
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
                
                a = app.probe(2).wiring(app.probe(2).wiring(:,5) == 1,[1 4]);
                [~,idx] = sort(a(:,2),'descend');
                app.probe(2).config = a(idx,1);
                
                a = app.probe(2).wiring(app.probe(2).wiring(:,5) == 2,[1 4]);
                [~,idx] = sort(a(:,2),'descend');
                app.probe(2).config = [app.probe(2).config; a(idx,1)]';
            end
            
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
            %resets probe and channel number to avoid confusion
            app.currentProbeNum = 3-app.currentProbeNum;
            if app.currentProbeNum==1
                app.currentChannel=1;
            else
                app.currentChannel=app.probe(1).nChannels+1;
            end

            app.label_selectedProbe.Text = num2str(app.currentProbeNum);
            app.refreshChannel();
            app.refreshProbe();
        end

        % Button pushed function: button_prevProbe
        function button_prevProbeButtonPushed(app, event)
            %resets probe and channel number to avoid confusion
            app.currentProbeNum = 3-app.currentProbeNum;
            if app.currentProbeNum==1
                app.currentChannel=1;
            else
                app.currentChannel=app.probe(1).nChannels+1;
            end
            
            app.label_selectedProbe.Text = num2str(app.currentProbeNum);
            app.refreshChannel();
            app.refreshProbe();
        end

        % Value changed function: slider_threshold
        function slider_thresholdChanged(app, event)
            app.thresholds(app.currentChannel) = event.Value;
            app.plotThreshold();
        end

        % Value changed function: text_minThreshold
        function text_minThresholdValueChanged(app, event)
            currentLimits=app.slider_threshold.Limits;
            currentLimits(1)=app.text_minThreshold.Value;
            tSteps=[round(currentLimits(1),-2):100:round(currentLimits(2),-2)]; %round with -2 rounds to 100
            app.slider_threshold.Limits=currentLimits;
            app.slider_threshold.MajorTicks=tSteps;
        end

        % Value changed function: text_maxThreshold
        function text_maxThresholdValueChanged(app, event)
            currentLimits=app.slider_threshold.Limits;
            currentLimits(2)=app.text_maxThreshold.Value;
            tSteps=[round(currentLimits(1),-2):100:round(currentLimits(2),-2)]; %round with -2 rounds to 100
            app.slider_threshold.Limits=currentLimits;
            app.slider_threshold.MajorTicks=tSteps;
        end

        % Value changed function: slider_time
        function slider_timeValueChanged(app, event)
            %update snippet
            app.refreshSnippet(true);
            
            %update text in 'go to' window
            app.text_gotoTime.Value=round(app.slider_time.Value);
            
            %figure out whether we need to renew the large time course
            if app.slider_time.Value<app.currentStart || app.slider_time.Value>app.currentStart+app.settings.lengthTC
                %set start point to fall into the middle of the time course
                %window (only if no conflict with borders)
                app.currentStart=app.slider_time.Value-app.settings.lengthTC/2;
                app.currentStart=max(0,app.currentStart);
                maxTime=floor((app.nTotalSamples-app.sampleFreq)/app.sampleFreq);
                app.currentStart=min(maxTime-app.settings.lengthTC,app.currentStart);
                
                app.getNewTimeCourse();
                app.refreshTimeCourse();
            end
            app.refreshTimeMarker();
        end

        % Value changed function: text_gotoTime
        function text_gotoTimeValueChanged(app, event)
            app.slider_time.Value = app.text_gotoTime.Value;
            
            app.refreshSnippet(true);
            
            %figure out whether we need to renew the large time course
            if app.slider_time.Value<app.currentStart || app.slider_time.Value>app.currentStart+app.settings.lengthTC
                %set start point to fall into the middle of the time course
                %window (only if no conflict with borders)
                app.currentStart=app.slider_time.Value-app.settings.lengthTC/2;
                app.currentStart=max(0,app.currentStart);
                maxTime=floor((app.nTotalSamples-app.sampleFreq)/app.sampleFreq);
                app.currentStart=min(maxTime-app.settings.lengthTC,app.currentStart);
                
                app.getNewTimeCourse();
                app.refreshTimeCourse();
            end
            app.refreshTimeMarker();
        end

        % Button pushed function: button_nextChan
        function button_nextChanButtonPushed(app, event)
            if app.currentChannel < app.nChannels
                app.currentChannel = app.currentChannel + 1;
                
                if app.currentChannel > app.probe(1).nChannels && app.currentProbeNum == 1
                    app.currentProbeNum = 2;
                    app.refreshProbe();
                end
                
                app.refreshChannel();
                app.refreshSnippet(false);
                app.getNewTimeCourse();
                app.refreshTimeCourse();
            end
        end

        % Button pushed function: button_prevChan
        function button_prevChanButtonPushed(app, event)
            if app.currentChannel > 1
                app.currentChannel = app.currentChannel - 1;
                
                if app.currentChannel <= app.probe(1).nChannels && app.currentProbeNum == 2
                    app.currentProbeNum = 1;
                    app.refreshProbe();
                end
                
                app.refreshChannel();
                app.refreshSnippet(false);
                app.getNewTimeCourse();
                app.refreshTimeCourse();
            end
        end

        % Value changed function: text_channelNum
        function text_channelNumValueChanged(app, event)
            app.currentChannel = app.text_channelNum.Value;
            %limit control
            if app.currentChannel<1
                app.currentChannel=1;
                app.text_channelNum.Value=1;
            end
            if app.currentChannel>app.nChannels
                app.currentChannel=app.nChannels;
                app.text_channelNum.Value=app.currentChannel;
            end
            
            %handle probe changes
            if length(app.probe)>1
                if app.currentChannel>app.probe(1).nChannels
                    oldProbe=app.currentProbeNum;
                    app.currentProbeNum = 2;
                else
                    oldProbe=app.currentProbeNum;
                    app.currentProbeNum = 1;
                end
                if oldProbe~=app.currentProbeNum
                    app.refreshProbe();
                end
            end
            
            app.refreshChannel();
            app.refreshSnippet(false);
            app.getNewTimeCourse();
            app.refreshTimeCourse();
        end

        % Value changed function: check_channel
        function check_channelValueChanged(app, event)
            app.chkChannels(app.currentChannel)=app.check_channel.Value;
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
            app.table_fileInfo.Data{6,2}{1} = num2str(sum(app.badChannels));
            app.refreshSnippet(false)
        end

        % Button pushed function: button_applyRMS
        function button_applyRMSButtonPushed(app, event)
            %compute RMS for each channel, using nr snippets specified in
            %settings
            
            %set lamp to clear computation is being performed
            app.lamp_applyRMS.Color=[1 0 0];
            
            %first part only necessary if RMS has not yet been computed
            if isempty(app.snippetRMS)
                RMSsnippets=zeros(app.nChannels,app.sampleFreq,app.settings.nrSnippetsRMS);
                
                dataFileName = [app.settings.expFolder '/' app.exptName '/' app.exptName '_amplifier.dat'];
                DataFile = fopen(dataFileName,'r');
                
                startSample=linspace(0,app.nTotalSamples-app.sampleFreq,app.settings.nrSnippetsRMS);
                
                for s=1:app.settings.nrSnippetsRMS
                    frewind(DataFile);
                    fseek(DataFile,2*startSample(s)*app.nChannels,'bof');
                    RMSsnippets(:,:,s) = fread(DataFile, [app.nChannels app.sampleFreq], 'int16');
                end
                
                if app.settings.useGPU==1      %use GPU array to speed things up
                    %filter snippets
                    filterSnippets = filter(app.butterParams.b1, app.butterParams.a1, gpuArray(RMSsnippets),[],2);
                    %compute rms
                    app.snippetRMS=gather(rms(filterSnippets,[2 3]));
                else
                    %filter snippets
                    filterSnippets = filter(app.butterParams.b1, app.butterParams.a1, RMSsnippets,[],2);
                    %compute rms
                    app.snippetRMS=rms(filterSnippets,[2 3]);
                end
            end
            
            %now set threshold according to RMS
            app.thresholds=app.spinner_RMS.Value.*app.snippetRMS;
            if app.thresholds(app.currentChannel)<app.slider_threshold.Limits(1)
                app.slider_threshold.Limits(1)=floor(app.thresholds(app.currentChannel)/100)*100;
            end
            if app.thresholds(app.currentChannel)>app.slider_threshold.Limits(2)
                app.slider_threshold.Limits(2)=ceil(app.thresholds(app.currentChannel)/100)*100;
            end
            app.slider_threshold.Value=app.thresholds(app.currentChannel);
            app.plotThreshold();
            
            %update lamp to indicate computation done
            app.lamp_applyRMS.Color=[0 1 0];
            
        end

        % Button pushed function: button_save
        function button_saveButtonPushed(app, event)
            %check that probe is fully set before saving
            if strcmp(app.dropdown_probe1.Value,'select')
                uialert(app.ThresholdGUIUIFigure,'Set probe before saving','Error');
                return;
            end
            
            if strcmp(app.dropdown_area1.Value,'select')
                uialert(app.ThresholdGUIUIFigure,'Set area before saving','Error');
                return;
            end
            
            
            %this needs to be here (even if doesn't have anything to do
            %with spike sorting) because the next step may run multiple instances in parallel,
            %so cannot generate a directory
            if ~exist([app.settings.expFolder '/' app.exptName '/SpikeFiles'],'dir')
                mkdir([app.settings.expFolder '/' app.exptName '/SpikeFiles'])
            end
            
            % for backward compatibility, for now
            %experiment = app.exptName;
            %Th = app.thresholds(~app.badChannels);
            %BadCh = app.badChannels;
            %b1 = app.butterParams.b1;
            %a1 = app.butterParams.a1;
            %isBR = 0;
            %probes = {app.probe.type};
            %if length(app.probe) == 1
            %    CHs = app.probe(1).config;
            %else
            %    CHs = [app.probe(1).config app.probe(2).config];
            %end
            %xPosition=app.probe(1).x;
            %yPosition=app.probe(1).z;
            %save([app.settings.expFolder '/' app.exptName '/experiment.mat'],'experiment','Th','CHs','BadCh','b1','a1','isBR','probes',...
            %    'xPosition','yPosition');
            
            % the new way
            thresholding.thresholds = app.thresholds;
            thresholding.badChannels = app.badChannels;
            thresholding.butter = app.butterParams;
            thresholding.chkChannels = app.chkChannels;
            
            id.exptId = app.exptName;
            id.probes = app.probe;
            id.name = app.settings.name;
            %id.email = app.settings.email;
            id.isBR = false;
            id.sampleFreq = app.sampleFreq;
            
            save([app.settings.expFolder '/' app.exptName '/' app.exptName '_threshold.mat'],'thresholding');
            save([app.settings.expFolder '/' app.exptName '/' app.exptName '_id.mat'],'id');
            
            msgbox('Thresholds and id saved.','Saved','non-modal')
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ThresholdGUIUIFigure and hide until all components are created
            app.ThresholdGUIUIFigure = uifigure('Visible', 'off');
            app.ThresholdGUIUIFigure.AutoResizeChildren = 'off';
            app.ThresholdGUIUIFigure.Position = [100 100 1383 578];
            app.ThresholdGUIUIFigure.Name = 'Threshold GUI';
            app.ThresholdGUIUIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.ThresholdGUIUIFigure);
            app.GridLayout.ColumnWidth = {287, '1x'};
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
            app.panel_location.Position = [7 410 269 151];

            % Create button_getExpt
            app.button_getExpt = uibutton(app.panel_location, 'push');
            app.button_getExpt.ButtonPushedFcn = createCallbackFcn(app, @button_getExptPushed, true);
            app.button_getExpt.Position = [10 98 100 22];
            app.button_getExpt.Text = 'Get Experiment';

            % Create button_save
            app.button_save = uibutton(app.panel_location, 'push');
            app.button_save.ButtonPushedFcn = createCallbackFcn(app, @button_saveButtonPushed, true);
            app.button_save.Enable = 'off';
            app.button_save.Position = [144 98 116 22];
            app.button_save.Text = 'Save Thresholds';

            % Create ProbeAreaLabel_2
            app.ProbeAreaLabel_2 = uilabel(app.panel_location);
            app.ProbeAreaLabel_2.HorizontalAlignment = 'right';
            app.ProbeAreaLabel_2.Position = [8 39 66 22];
            app.ProbeAreaLabel_2.Text = 'Probe/Area';

            % Create dropdown_probe2
            app.dropdown_probe2 = uidropdown(app.panel_location);
            app.dropdown_probe2.Items = {'select', '64D', '64E', '64F', '64G', '64H', '64M', '128P', 'single', 'tetrode', ''};
            app.dropdown_probe2.ValueChangedFcn = createCallbackFcn(app, @dropdown_probe2ValueChanged, true);
            app.dropdown_probe2.Enable = 'off';
            app.dropdown_probe2.Position = [89 39 100 22];
            app.dropdown_probe2.Value = 'select';

            % Create ProbeAreaLabel
            app.ProbeAreaLabel = uilabel(app.panel_location);
            app.ProbeAreaLabel.HorizontalAlignment = 'right';
            app.ProbeAreaLabel.Position = [8 69 66 22];
            app.ProbeAreaLabel.Text = 'Probe/Area';

            % Create dropdown_probe1
            app.dropdown_probe1 = uidropdown(app.panel_location);
            app.dropdown_probe1.Items = {'select', '64D', '64E', '64F', '64G', '64H', '64M', '128P', 'single', 'tetrode'};
            app.dropdown_probe1.ValueChangedFcn = createCallbackFcn(app, @dropdown_probe1ValueChanged, true);
            app.dropdown_probe1.Enable = 'off';
            app.dropdown_probe1.Position = [89 69 100 22];
            app.dropdown_probe1.Value = 'select';

            % Create dropdown_area1
            app.dropdown_area1 = uidropdown(app.panel_location);
            app.dropdown_area1.Items = {'select', 'V1', 'PSS', 'LGN', 'S1', 'other'};
            app.dropdown_area1.ValueChangedFcn = createCallbackFcn(app, @dropdown_area1ValueChanged, true);
            app.dropdown_area1.Enable = 'off';
            app.dropdown_area1.Position = [197 69 63 22];
            app.dropdown_area1.Value = 'select';

            % Create dropdown_area2
            app.dropdown_area2 = uidropdown(app.panel_location);
            app.dropdown_area2.Items = {'select', 'V1', 'PSS', 'LGN', 'S1', 'other'};
            app.dropdown_area2.ValueChangedFcn = createCallbackFcn(app, @dropdown_area2ValueChanged, true);
            app.dropdown_area2.Enable = 'off';
            app.dropdown_area2.Position = [197 39 63 22];
            app.dropdown_area2.Value = 'select';

            % Create SettingsButton
            app.SettingsButton = uibutton(app.panel_location, 'push');
            app.SettingsButton.ButtonPushedFcn = createCallbackFcn(app, @button_settingsButtonPushed, true);
            app.SettingsButton.Position = [10 8 100 22];
            app.SettingsButton.Text = 'Settings';

            % Create panel_fileInfo
            app.panel_fileInfo = uipanel(app.LeftPanel);
            app.panel_fileInfo.AutoResizeChildren = 'off';
            app.panel_fileInfo.Title = 'File Details';
            app.panel_fileInfo.Position = [7 13 269 388];

            % Create table_fileInfo
            app.table_fileInfo = uitable(app.panel_fileInfo);
            app.table_fileInfo.ColumnName = {'Name'; 'Prop'};
            app.table_fileInfo.RowName = {''; ''};
            app.table_fileInfo.Position = [10 15 250 348];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Scrollable = 'on';

            % Create panel_signal
            app.panel_signal = uipanel(app.RightPanel);
            app.panel_signal.Position = [296 13 788 548];

            % Create axis_snippet
            app.axis_snippet = uiaxes(app.panel_signal);
            xlabel(app.axis_snippet, 'Time (s)')
            ylabel(app.axis_snippet, 'Signal (mV)')
            app.axis_snippet.PlotBoxAspectRatio = [1.31331168831169 1 1];
            app.axis_snippet.GridAlpha = 0.15;
            app.axis_snippet.MinorGridAlpha = 0.25;
            app.axis_snippet.Box = 'on';
            app.axis_snippet.Position = [12 155 485 357];

            % Create ThresholdSliderLabel
            app.ThresholdSliderLabel = uilabel(app.panel_signal);
            app.ThresholdSliderLabel.HorizontalAlignment = 'center';
            app.ThresholdSliderLabel.VerticalAlignment = 'top';
            app.ThresholdSliderLabel.Position = [496 514 76 22];
            app.ThresholdSliderLabel.Text = 'Threshold';

            % Create slider_threshold
            app.slider_threshold = uislider(app.panel_signal);
            app.slider_threshold.Limits = [-500 500];
            app.slider_threshold.Orientation = 'vertical';
            app.slider_threshold.ValueChangedFcn = createCallbackFcn(app, @slider_thresholdChanged, true);
            app.slider_threshold.Enable = 'off';
            app.slider_threshold.Position = [501 198 3 282];
            app.slider_threshold.Value = -200;

            % Create button_nextChan
            app.button_nextChan = uibutton(app.panel_signal, 'push');
            app.button_nextChan.ButtonPushedFcn = createCallbackFcn(app, @button_nextChanButtonPushed, true);
            app.button_nextChan.Enable = 'off';
            app.button_nextChan.Position = [187 516 23 23];
            app.button_nextChan.Text = '>';

            % Create button_prevChan
            app.button_prevChan = uibutton(app.panel_signal, 'push');
            app.button_prevChan.ButtonPushedFcn = createCallbackFcn(app, @button_prevChanButtonPushed, true);
            app.button_prevChan.Enable = 'off';
            app.button_prevChan.Position = [119 516 25 23];
            app.button_prevChan.Text = '<';

            % Create button_rejectChan
            app.button_rejectChan = uibutton(app.panel_signal, 'push');
            app.button_rejectChan.ButtonPushedFcn = createCallbackFcn(app, @button_rejectChanPushed, true);
            app.button_rejectChan.Enable = 'off';
            app.button_rejectChan.Position = [220 516 117 22];
            app.button_rejectChan.Text = 'Reject Channel';

            % Create axis_spikePreview
            app.axis_spikePreview = uiaxes(app.panel_signal);
            xlabel(app.axis_spikePreview, 'Time (ms)')
            ylabel(app.axis_spikePreview, 'Signal (mV)')
            app.axis_spikePreview.PlotBoxAspectRatio = [1.31331168831169 1 1];
            app.axis_spikePreview.GridAlpha = 0.15;
            app.axis_spikePreview.MinorGridAlpha = 0.25;
            app.axis_spikePreview.Box = 'on';
            app.axis_spikePreview.Position = [559 231 219 188];

            % Create text_channelNum
            app.text_channelNum = uieditfield(app.panel_signal, 'numeric');
            app.text_channelNum.ValueChangedFcn = createCallbackFcn(app, @text_channelNumValueChanged, true);
            app.text_channelNum.HorizontalAlignment = 'center';
            app.text_channelNum.Enable = 'off';
            app.text_channelNum.Position = [147 516 37 22];
            app.text_channelNum.Value = 1;

            % Create ChannelLabel
            app.ChannelLabel = uilabel(app.panel_signal);
            app.ChannelLabel.HorizontalAlignment = 'center';
            app.ChannelLabel.Position = [52 516 54 22];
            app.ChannelLabel.Text = 'Channel:';

            % Create SegmentstartsLabel
            app.SegmentstartsLabel = uilabel(app.panel_signal);
            app.SegmentstartsLabel.HorizontalAlignment = 'center';
            app.SegmentstartsLabel.Position = [7 25 57 28];
            app.SegmentstartsLabel.Text = {'Segment '; 'start (s)'};

            % Create slider_time
            app.slider_time = uislider(app.panel_signal);
            app.slider_time.ValueChangedFcn = createCallbackFcn(app, @slider_timeValueChanged, true);
            app.slider_time.Enable = 'off';
            app.slider_time.Position = [77 45 457 3];

            % Create text_minThreshold
            app.text_minThreshold = uieditfield(app.panel_signal, 'numeric');
            app.text_minThreshold.ValueChangedFcn = createCallbackFcn(app, @text_minThresholdValueChanged, true);
            app.text_minThreshold.HorizontalAlignment = 'center';
            app.text_minThreshold.Enable = 'off';
            app.text_minThreshold.Position = [504 159 41 22];
            app.text_minThreshold.Value = -500;

            % Create text_maxThreshold
            app.text_maxThreshold = uieditfield(app.panel_signal, 'numeric');
            app.text_maxThreshold.ValueChangedFcn = createCallbackFcn(app, @text_maxThresholdValueChanged, true);
            app.text_maxThreshold.HorizontalAlignment = 'center';
            app.text_maxThreshold.Enable = 'off';
            app.text_maxThreshold.Position = [504 495 41 22];
            app.text_maxThreshold.Value = 500;

            % Create PresetthresholdPanel
            app.PresetthresholdPanel = uipanel(app.panel_signal);
            app.PresetthresholdPanel.Title = 'Preset threshold';
            app.PresetthresholdPanel.Position = [609 444 144 94];

            % Create RMSSpinnerLabel
            app.RMSSpinnerLabel = uilabel(app.PresetthresholdPanel);
            app.RMSSpinnerLabel.HorizontalAlignment = 'right';
            app.RMSSpinnerLabel.Position = [12 43 32 22];
            app.RMSSpinnerLabel.Text = 'RMS';

            % Create spinner_RMS
            app.spinner_RMS = uispinner(app.PresetthresholdPanel);
            app.spinner_RMS.Step = 0.5;
            app.spinner_RMS.Limits = [-6 6];
            app.spinner_RMS.Enable = 'off';
            app.spinner_RMS.Position = [59 43 69 22];
            app.spinner_RMS.Value = -3;

            % Create button_applyRMS
            app.button_applyRMS = uibutton(app.PresetthresholdPanel, 'push');
            app.button_applyRMS.ButtonPushedFcn = createCallbackFcn(app, @button_applyRMSButtonPushed, true);
            app.button_applyRMS.Enable = 'off';
            app.button_applyRMS.Position = [20 11 80 22];
            app.button_applyRMS.Text = 'Apply';

            % Create lamp_applyRMS
            app.lamp_applyRMS = uilamp(app.PresetthresholdPanel);
            app.lamp_applyRMS.Position = [108 12 20 20];

            % Create axis_timecourse
            app.axis_timecourse = uiaxes(app.panel_signal);
            title(app.axis_timecourse, '')
            xlabel(app.axis_timecourse, 'Time (s)')
            ylabel(app.axis_timecourse, 'Signal')
            app.axis_timecourse.PlotBoxAspectRatio = [9.33802816901408 1 1];
            app.axis_timecourse.Position = [12 62 490 90];

            % Create GotosEditFieldLabel
            app.GotosEditFieldLabel = uilabel(app.panel_signal);
            app.GotosEditFieldLabel.HorizontalAlignment = 'center';
            app.GotosEditFieldLabel.Position = [520 83 52 22];
            app.GotosEditFieldLabel.Text = 'Go to (s)';

            % Create text_gotoTime
            app.text_gotoTime = uieditfield(app.panel_signal, 'numeric');
            app.text_gotoTime.ValueChangedFcn = createCallbackFcn(app, @text_gotoTimeValueChanged, true);
            app.text_gotoTime.HorizontalAlignment = 'center';
            app.text_gotoTime.Position = [521 62 51 22];

            % Create check_channel
            app.check_channel = uicheckbox(app.panel_signal);
            app.check_channel.ValueChangedFcn = createCallbackFcn(app, @check_channelValueChanged, true);
            app.check_channel.Text = 'Checked channel';
            app.check_channel.Position = [362 516 115 22];

            % Create panel_probe
            app.panel_probe = uipanel(app.RightPanel);
            app.panel_probe.AutoResizeChildren = 'off';
            app.panel_probe.Position = [10 13 269 548];

            % Create axis_probe
            app.axis_probe = uiaxes(app.panel_probe);
            title(app.axis_probe, '')
            xlabel(app.axis_probe, 'X')
            ylabel(app.axis_probe, 'Y')
            app.axis_probe.PlotBoxAspectRatio = [1 1.53 1];
            app.axis_probe.Position = [10 15 249 493];

            % Create label_selectedProbe
            app.label_selectedProbe = uilabel(app.panel_probe);
            app.label_selectedProbe.HorizontalAlignment = 'center';
            app.label_selectedProbe.Position = [148 513 27 24];
            app.label_selectedProbe.Text = '1';

            % Create button_nextProbe
            app.button_nextProbe = uibutton(app.panel_probe, 'push');
            app.button_nextProbe.ButtonPushedFcn = createCallbackFcn(app, @button_nextProbeButtonPushed, true);
            app.button_nextProbe.Enable = 'off';
            app.button_nextProbe.Position = [182 514 23 23];
            app.button_nextProbe.Text = '>';

            % Create button_prevProbe
            app.button_prevProbe = uibutton(app.panel_probe, 'push');
            app.button_prevProbe.ButtonPushedFcn = createCallbackFcn(app, @button_prevProbeButtonPushed, true);
            app.button_prevProbe.Enable = 'off';
            app.button_prevProbe.Position = [115 514 25 23];
            app.button_prevProbe.Text = '<';

            % Create ProbeLabel
            app.ProbeLabel = uilabel(app.panel_probe);
            app.ProbeLabel.HorizontalAlignment = 'center';
            app.ProbeLabel.Position = [55 514 41 22];
            app.ProbeLabel.Text = 'Probe:';

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