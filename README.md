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



