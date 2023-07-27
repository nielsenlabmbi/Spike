# Spike

## Spike processing pipeline:

1) thresholdGui (GUI) and thresholdGuiSettings (GUI):
   - visually set threshold for each channel, mark bad channels
   - uses median approach to calculate a threshold guess for every channel
   - needs read_Intan_Header
   - note on probe configurations: code needs the probeConfig folder. that folder should contain a file for every probe configuration called probeConfig_type, where type is the name used for the probe in the id file (such as 64D). probeConfig files are functions that generate a matrix with 5 columns (column 1: channel number, column 2: x position, column 3: y, column 4: z, column 5: shank number). The threshold gui will automatically list all probes for which probeConfig files exist.
   - note: threshold data contains data for 1 probe only; id file is generated with probe information for all probes present
   *Output: filename_pX_threshold.mat and filename_id.mat (locally and on Z if selected); generates SpikeFiles folder*

2) extractSpikes:
   - run extractSpikes to extract waveforms for each channel for 1 probe (use parfor for speed)
   - implements temporal (and optionally spatial) constraints
   - note: run with parfor loop to speed up processing time
   - note: first job should be 0, not 1
   - note: if you are regenerating files for a previously sorted file, check whether you need to set the legacyFlag to 1 (this determines the amount of overlap between jobs). You can tell whether this needs to be the case by loading a single job file and looking at the settings structure. If it does not contain settings.offsetSamples, or if settings.legacyFlag=1, then you need to set legacyFlag to 1.
   *Output: SpikeFiles\filename_pX_jID_spike.mat; updates to _id.mat file (local and on Z if selected)*

3) extractSpikeProps:
   - extract set of properties for each waveform on one probe
   - note: use parfor to speed up processing time
   *Output: SpikeFiles\filename_jID_pX_spkinfo.mat; updates to _id.mat file (local and on Z if selected)*

4) sortGUI (GUI):
   - sort data
   - load data file after specifying which jobs to load for which probe, as well as the matching ID file
   - Views:
     - 'probe mode' plots 2 waveform parameters against each other for all (or a subset of) channels; the number of channels displayed can be limited using the button (maintains probe display format, but only shows events detected at the selected channels)
     - 'tetrode mode' limits the plot to waveforms detected on one channel, and plots parameters computed at the time of those events based on the simultaneous recordings from different nearby channels (whether or not they crossed the threshold on those channels)
     - 'channel mode' limits data to all events detected on a single channel
   - artefact rejection allows detection of events that occur on many channels simultaneously (specify the number of channels and an optional threshold to set what is considered as an artefact) 
   - note: you can pan, zoom in and out of the main plot using the buttons on the top right
   - when adding units or ROIs, double-click into each roi after drawing to end the drawing
   - info in unit table: unit: unit number and color; ch: main channel at which unit is detected; vis: unit visible or invisible in plot; cat: su, mu or noise; ISIv: percentage of events with ISIs below 1.2ms 
   - waveform plot: choose up to 2 units to display; can add additional waveforms for comparison; nr spikes controls how many spikes are shown per unit
   - waveform reader: click to show the reader (circle in main plot); double click to activate the display of a waveform close to that spot; reader can be dragged around the display
   - unit footprint: distribution of detection channels for a selected unit
   - output is spkSort structure, which in addition to general info contains *unitid*: vector with unit assignment for each time stamp (-1 artefact, no distinction between SU, MU and noise), *spktimes*: vector with all spike time stamps, *unitinfo*: cell array with type assignemt (SU, MU, noise, none) for each unitid
   *Output: filename_pX_spkSort.mat (local and on Z); updates to _id.mat file (local and on Z) and database*



## Further processing for sorted spike data:

1) extractTrials (in spkFunctions)
   - needs analyzer git repository in path (needs AnalyzerUtils)
   - extracts trial info from Analyzer: stim condition in each trial
   - extracts events and their timing from digital file
   - output is trialInfo structure
   - note that eventId=0 corresponds to end of trials only
   *Output: filename_trialInfo.mat*

2) SUTrialData (in spkFunctions)
   - only uses data for sorted single units
   - note: units marked 'none' or 'noise' are excluded
   - extracts spike times in a window around a selected event
   - time of spikes is reported relative to event in ms 
   - also computes mean number of spikes before and after event
   - adds trial info from Analyzer (copy from trialInfo structure) 
   *Output: filename_pX_SUTrial.mat* 

3) computePSTH (in spkFunctions)
   - compute PSTH for a unit
   - uses the data computed from SUTrialData 
   - note: the same funciton can be used for MUThresh data as well by supplying the data for a channel rather than a single unit (as provided by MUTrialData)


## Processing of multi-unit data
There are 2 kinds of multi-unit data that can be analyzed:
- MUThresh: Thresholded but not spike sorted data
- MUEnv: 'Nikos' style analysis using the downsampled, rectified continuous signal
Lastly, MURaw will refer to the underlying original multi-unit signal (the bandpass filtered signal from every channel that is used to compute the other 2 MU signals).

### MUThresh analysis pipeline

1) computeMUThreshold (in MU analysis)
- computes automated threshold for every channel
 *Output: filename_MUthreshold.mat and updates to _id.mat file*

2) extractSpikes, extractSpikeProps
- run usual pipeline of extractSpikes and extractSpikeProps, setting the MUflag to 1 so that it uses the MUthresholding file and generates the right output
 *Output: SpikeFiles\filename_pX_jID_MUspike.mat, SpikeFiles\filename_pX_jID_MUspkinfo.mat, updates to _id.mat file*

3) mergeMUspkInfo (in MU analysis)
- merges the separate MUspkinfo jobs, extracting spike times and detection channels for further analysis
*Output: filename_MUspkMerge.mat*

4) extractTrials (in spkFunctions)
- see above for SU analysis
- note: this file only needs to be generated once and then works for SU and MU data

5) MUThreshTrialData (in MU analysis)
- takes the trialInfo and MUspkMerge data files
- reorganizes into a structure MUThresh, containing information relative to a specific event for every channel and every trial
- computes number of spikes and firing rates in baseline, stimulus period and entire trial
- adds trial info from Analyzer
*Output: filename_pX_MUThreshTrial.mat* 


### MUEnv analysis pipeline

1) MUEnvTrialData (in MU analysis)
   - processes data for all channels, in a time window around a selected event
   - output is the rectified and downsampled signal for every channel for each event
   - also computes z transform of the signal
*Output: filename_pX_MUEnvTrial.mat* 

### MURaw analysis pipeline
1) MURawTrialData (in MU analysis)
   - processes data for single channel only
   - extracts the filtered data for that channel in a time window around a selected event
   - note: channel number provided as input refers to the original recording channel
*Output: filename__cCh_pX_MURawTrial.mat* 


## Using the new spike sorting pipeline with data sorted using the old pipeline:
The new spike sorting GUI can be used with data processed with the old spike sorting process by following these steps (all of these files can be found in util):
1) Make sure a spkSort file has been generated for the old data
2) Compute the raw (Intan) detection channel:
   - function: detChOldFormat 
   - In the old system, channels were sorted according to their position on the probe before extracting spikes. The tranformation from recording channel to detection channel was hard coded, and differs slightly from our current method. Depending on the time at which the data was sorted, the spike/Spike file may contain the recording channel in addition to the sorting channel, but the detCh entry always refers to the sorted channels. In addition, it also removes channels that were marked as bad, which may be recorded in the experiment.mat file (which doesn't exist for all files, and was overwritten in the case of sorting files with multiple probes).
   - To compute the original detection channel therefore requires a few assumptions, and differs based on what is available. It is important to check that this step is done correctly by plotting the data.
   - The function adds a removes 'detChSort' from spkSort to avoid conflicts with the detChSort assigned in spike Sorting.
 3) Extract spike waveforms:
   - function: extractSpikesOldFormat
   - This function uses the detCh info added to spkSort in step 2 to extract waveforms from the amplifier file at the time stamps saved in spkSort (which are the ones contained in the original spike or Spike file). It uses the same parameters as extractSpikes (radius and spikeSamples), and organizes channels the same way as used in the new pipeline.
   - Before running the code, make sure the SpikeFiles folder exists (usually generated during thresholding)
 4) Extract spike properties:
   - function: extractSpikeProps
   - same function used normally
 5) Reorganize spkSort order
   - function: resortSpkSort
   - this step is crucial to avoid having units matched to the wrong timestamps; it reorganizes the entries in spkSort according to job, channel and timestamp rather than timestamp
   - also adds a few other fields that are needed by the sortGui
 6) Sort:
   - use same sortGui as always
   - Note: While we are keeping the same time stamps and therefore extracting the same basic waveforms as in the original implementation, there are some differences in how things are computed, most notably anything that has to do with channel position (such as the center of mass for the channel)
   - Also note that this obviously does not change how the channels were assigned during initial sorting (with a minimum applied across channels instead of keeping channels indepedent as in the new version)
 7) Divide into trials:
   - function: extractTrialsBR
   - This function addresses all of the issues with the digital pulses (hopefully), and saves the trial events in the same way as for the intan data
    

## Dealing with blackrock files
To sort files acquired with Blackrock, they need to be converted into files that follow the same format as the intan files; at that point, the same code as always can be applied. 
1) convertBrIntan
   - provide name of blackrock file
   - generates a fake header file, and an amplifier file
*Output: filename_info.rhd, filename_amplifier.dat*


