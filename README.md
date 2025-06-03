# Spike

## Spike processing pipeline:

1) probe configuration file
   - need a configuration file for every probe, stored in the probeConfig folder
   - file are called probeConfig_type, where type is the name of the probe
   - probeConfig files are functions that generate a matrix with 5 columns (column 1: channel number, column 2: x position, column 3: y, column 4: z, column 5: shank number)
   - convention for Z: 0 = bottom; positive numbers = channels higher on probe
   - can visualize any configuration file by calling probeViewer(type)
   _Output: probeConfig_XXX.mat_


2) thresholdGui (GUI) and thresholdGuiSettings (GUI):
   - visually set threshold for each channel, mark bad channels
   - uses median approach to calculate a threshold guess for every channel
   - needs read_Intan_Header
   - needs configuration file (the threshold gui will automatically list all probes for which probeConfig files exist)
   - note: threshold data contains data for 1 probe only; id file is generated with probe information for all probes present
   - multiple versions are possible and will be marked with a suffix for the threshold file; each will generate its own SpikeFiles folder
   _Output: filename_pX_threshold*.mat and filename_id.mat (locally and on Z if selected); generates SpikeFiles* folder_

3) extractSpikes:
   - run extractSpikes to extract waveforms for each channel for 1 probe (use parfor for speed)
   - implements temporal (and optionally spatial) constraints
   - takes 2 optional inputs:
     - id: for batch processing, the id file can be read outside the parfor loop and then provided as an input argument to the function (otherwise processing may stop because multiple processes attempt to read the id file)
     - tSuffix: suffix for threshold and SpikeFiles folder 
   - note: run with parfor loop to speed up processing time
   - note: first job should be 0, not 1
   - note: if you are regenerating files for a previously sorted file, check whether you need to set the legacyFlag to 1 (this determines the amount of overlap between jobs). You can tell whether this needs to be the case by loading a single job file and looking at the settings structure. If it does not contain settings.offsetSamples, or if settings.legacyFlag=1, then you need to set legacyFlag to 1.
   - note: different versions of spike files are only indicated by the folder they are in, not by their filename
   _Output: SpikeFiles*\filename_pX_jID_spike.mat; filename_pX_extractSpk*.mat file (local and on Z if selected)_

4) extractSpikeProps:
   - extract set of properties for each waveform on one probe
   - takes 2 optional inputs:
     - id: for batch processing, the id file can be read outside the parfor loop and then provided as an input argument to the function (otherwise processing may stop because multiple processes attempt to read the id file)
     - tSuffix: suffix for threshold and SpikeFiles folder 
   - note: use parfor to speed up processing time
   - note: different versions of spike files are only indicated by the folder they are in, not by their filename
   _Output: SpikeFiles\filename_jID_pX_spkinfo.mat; filename_pX_extractSpkProp*.mat file (local and on Z if selected)_

5) sortGUI (GUI):
   - sort data
   - load data file after specifying an id file or a spkSort file
   - for the id file, can limit the number of jobs and the range of jobs
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
   - refreshjobs reloads a random sample of jobs if less than 100% of jobs are loaded 
   - undo can undo the last action applied
   - apply history applies all of the steps applied to a previously sorted file to the current ones
   - output is spkSort structure, which in addition to general info contains *unitid*: vector with unit assignment for each time stamp (-1 artefact, no distinction between SU, MU and noise), *spktimes*: vector with all spike time stamps, *unitinfo*: cell array with type assignemt (SU, MU, noise, none) for each unitid
   - note: different versions are possible by loading different spkinfo files (determined when loading from id file); different spkSort files can also be saved by changing the suffix for the spkSort file (does not have to match the suffix for the spkinfo files)
   - note: if there is an error about a missing field in spkSort, it's due to using the old id file format. in this case, first run splitIdFile to convert the old id file format to the new one (adds entries to spkSort, e.g.)
   _Output for 100% jobs: filename_pX_spkSort*.mat (local and on Z); filename_pX_sortHist*.mat (local and on Z); updates to database_
   _Output for <100% jobs: filename_pX_partSpkSort*.mat (local and on Z); filename_pX_partSortHist*.mat (local and on Z); updates to database_



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
   - optional input argument allows different versions of input/output 
   _Output: filename_pX_SUTrial*.mat_

3) computePSTH (in spkFunctions)
   - compute PSTH for a unit
   - uses the data computed from SUTrialData 
   - note: the same funciton can be used for MUThresh data as well by supplying the data for a channel rather than a single unit (as provided by MUTrialData)


## Partial spike sorting pipeline (for long files)

1) sortGUI 
   - First sort using less than 100% of jobs, making sure by clicking refreshjobs that the sorting is representative. 
   - Also assign categories to all units. 
   _Output: filename_pX_partSpkSort*.mat, filename_pX_partSortHist*.mat_

2) applySortFast
   - apply partSortHist to each spkinfo job file
   - optional argument: file suffix
   - can be run in a parfor loop
   - creates spkSortP structure with fields unitid, spktimes, detCh, detChSort and info
   _Output: SpikeFiles*\filename_jID_pX_spkSort.mat_

3) mergeJobSort
   - merge the individual job spkSort files
   - creates spkSort structure containing the fields unitid, spktimes, detCh, detChSort, unitinfo, spkProps, info
   - note that spkProps (waveform info) is copied from the partSpkSort file and therefore only reflects properties of that random sample of spikes
   - unit categories are copied form the partSpkSort file as well (as set in the sortGui); SU that now have ISI violations are downgraded to MU in this process
   *Output: SpikeFiles\filename__pX_spkSort.mat; updates to _id.mat file (local and on Z if selected)*


## Spike sorting across separate files (all functions in util)

1) mergeIntan
   - merges intan amplifier files
   - merging occurs independently of probes
   - no id file needed for merging (will be created in next step)
   - also creates header for the merged file
   - collects all relevant merging information in structure mergeInfo
   - by default, the unit of the output file is set to uMMM to indicate a merged file; the experiment ID is specified by the user
    *Output: filename_uMMM_expId_amplifier.dat, filename_uMMM_expId_info.rhd, filename_uMMM_expId_mergeInfo.mat*

2) regular processing
   - proceed through the regular pipeline, using the threshold GUI to generate an ID file for the merged amplifier data
   - should end with spkSort file for merged data

3) splitIntan
   - uses the info saved in mergeInfo to split the sorted data into separate sort files
   - spkSort files will be generated in the correct folders for the original (pre-merged) data
   - it will overwrite existing spkSort files for the same probe only
   - also adds information to the id file (or generates id file) for those data files
   - unit assignment will be copied from the merged file
   *Output: for every file that was merged, filenames_pX_spkSort.mat and filenames_id.mat*


## Processing of multi-unit data

There are 2 kinds of multi-unit data that can be analyzed:
- MUThresh: Thresholded but not spike sorted data
- MUEnv: 'Nikos' style analysis using the downsampled, rectified continuous signal
Lastly, MURaw will refer to the underlying original multi-unit signal (the bandpass filtered signal from every channel that is used to compute the other 2 MU signals).

### MUThresh analysis pipeline

1) computeMUThreshold (in MU analysis)
- computes automated threshold for every channel
- optional argument: suffix (to allow different versions)
 _Output: filename_MUthreshold*.mat; generates SpikeFiles* folder_

2) extractSpikes, extractSpikeProps
- run usual pipeline of extractSpikes and extractSpikeProps, setting the MUflag to 1 so that it uses the MUthresholding file and generates the right output
- same optional arguments as for the SU pipeline to allow preloading of ID files and multiple versions of output files
 _Output: SpikeFiles*\filename_pX_jID_MUspike.mat, SpikeFiles*\filename_pX_jID_MUspkinfo.mat, filename_pX_extractSpk*.mat file, filename_pX_extractSpkProp*.mat file_

3) mergeMUspkInfo (in MU analysis)
- merges the separate MUspkinfo jobs, extracting spike times and detection channels for further analysis
- optional argument: suffix to indicate version
_Output: filename_MUspkMerge*.mat_

4) extractTrials (in spkFunctions)
- see above for SU analysis
- note: this file only needs to be generated once and then works for SU and MU data

5) MUThreshTrialData (in MU analysis)
- takes the trialInfo and MUspkMerge data files
- reorganizes into a structure MUThresh, containing information relative to a specific event for every channel and every trial
- computes number of spikes and firing rates in baseline, stimulus period and entire trial
- adds trial info from Analyzer
- optional input argument allows different versions of input/output 
_Output: filename_pX_MUThreshTrial*.mat_


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



