# Spike

## Spike processing pipeline:

1) thresholdGui (GUI) and thresholdGuiSettings (GUI):
   - visually set threshold for each channel, mark bad channels
   - uses median approach to calculate a threshold guess for every channel
   - needs microprobe_wiring directory and read_Intan_Header
   - note: threshold data contains data for 1 probe only; id file is generated with probe information for all probes present
   *Output: filename_pX_threshold.mat and filename_id.mat (locally and on Z if selected); generates SpikeFiles folder*

2) extractSpikes:
   - run extractSpikes to extract waveforms for each channel for 1 probe (use parfor for speed)
   - implements temporal (and optionally spatial) constraints
   - note: run with parfor loop to speed up processing time
   - note: first job should be 0, not 1
   *Output: SpikeFiles\filename_pX_jID_spike.mat; updates to _id.mat file (local and on Z if selected)*

3) extractSpikeProps:
   - extract set of properties for each waveform on one probe
   - note: use parfor to speed up processing time
   *Output: SpikeFiles\filename_jID_pX_spkinfo.mat; updates to _id.mat file (local and on Z if selected)*

4) sortGUI (GUI):
   - sort data
   - load data file after specifying which jobs to load for which probe, as well as the matching ID file
   - 'probe mode' plots 2 waveform parameters against each other for all (or a subset of) channels; 'tetrode mode' limits the plot to waveforms detected on one channel, and plots parameters computed for those waveforms based on the recordings from different nearby channels 
   - note: you can limit the number of channels for the 'probe mode' display (maintains probe display format, but only shows events detected at the selected channels)
   - artefact rejection allows detection of events that occur on many channels simultaneously (specify the number of channels and an optional threshold to set what is considered as an artefact) 
   - note: you can pan, zoom in and out of the main plot using the buttons on the top right
   - when adding units or ROIs, double-click into each roi after drawing to end the drawing
   - info in unit table: unit: unit number and color; ch: main channel at which unit is detected; vis: unit visible or invisible in plot; cat: su, mu or noise; ISIv: percentage of events with ISIs below 1.2ms 
   - waveform plot: choose up to 2 units to display; can add additional waveforms for comparison; nr spikes controls how many spikes are shown per unit
   - unit footprint: distribution of detection channels for a selected unit
   - output is spkSort structure, which in addition to general info contains *unitid*: vector with unit assignment for each time stamp (-1 artefact, no distinction between SU, MU and noise), *spktimes*: vector with all spike time stamps, *unitinfo*: cell array with type assignemt (SU, MU, noise, none) for each unitid
   *Output: filename_pX_spkSort.mat (local and on Z); updates to _id.mat file (local and on Z) and database*



## Further processing:
1) for unsorted (MUA) analysis - mergeSpkInfo (in spkFunctions):
   - merges the individual spkInfo files into one file for MUA analysis
   - keeps spktimes, detCh, detChSort in structure spkMerge 
   *Output: filename_pX_spkMerge.mat*

2) extractTrials (in spkFunctions)
   - needs analyzer git repository in path (needs AnalyzerUtils)
   - extracts trial info from Analyzer: stim condition in each trial
   - extracts events and their timing from digital file
   - output is trialInfo structure
   - note that eventId=0 corresponds to end of trials only
   *Output: filename_trialInfo.mat*

3) SUTrialData (in spkFunctions)
   - only uses data for sorted single units
   - note: units marked 'none' or 'noise' are excluded
   - extracts spike times in a window around a selected event
   - time of spikes is reported relative to event in ms 
   - also computes mean number of spikes before and after event
   - adds trial info from Analyzer (copy from trialInfo structure) 
   *Output: filename_pX_SUTrial.mat* 


## Using the new spike sorting pipeline with old data:
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


