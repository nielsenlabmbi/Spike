# Spike

## Spike processing pipeline:

1) thresholdGui (GUI) and thresholdGuiSettings (GUI):
   - generate settings file (saved locally) if necessary 
   - visually set threshold for each channel, mark bad channels
   - needs microprobe_wiring directory and read_Intan_Header
   - note: runs on all channels (across probes)
   *Output: filename_threshold.mat and filename_id.mat (locally and on Z if selected); generates SpikeFiles folder*

2) extractSpikes and extractSpikesSettings (GUI):
   - first generate settings file (saved locally) with general settings using GUI
   - then run extractSpikes to extract waveforms for each channel (use parfor for speed)
   - implements temporal (and optionally spatial) constraints
   - to run, use *extractSpikes(animalId,unitId,expId,settings,parts,jobId)*. Specify *unitId* without the 'u'. *settings* is the global settings structure generated with the GUI (make sure to load this into the workspace first), *parts* the number of segments into which to divice the data file, and *jobID* the current job
   - note: run with parfor loop to speed up processing time
   - note: first job should be 0, not 1
   - note: this runs on all channels, independent of probe
   *Output: SpikeFiles\filename_jID_spike.mat; updates to _id.mat file (local and on Z if selected)*

3) extractSpikeProps:
   - extract set of properties for each waveform
   - to run, use *extractSpikeProps(expFolder,animalId,unitId,ExpId,jobId,name)*, where *expFolder* is the base folder for physiology files, and *name* is the initials of the person executing the command
   - note: use parfor to speed up processing time
   - note: splits data into separate files for each probe
   *Output: SpikeFiles\filename_jID_pX_spkinfo.mat; updates to _id.mat file (local and on Z if selected)*

4) sortGUI (GUI):
   - sort data
   - load data file after specifying which jobs to load for which probe, as well as the matching ID file
   - 'probe mode' plots 2 waveform parameters against each other for all (or a subset of) channels; 'tetrode mode' limits the plot to waveforms detected on one channel, and plots parameters computed for those waveforms based on the recordings from different nearby channels 
   - note: you can limit the number of channels for the 'probe mode' display (maintains probe display format, but only shows events detected at the selected channels)
   - artefact rejection allows detection of events that occur on many channels simultaneously (specify the number of channels and an optional threshold to set what is considered as an artefact) 
   - note: you can pan, zoom in and out of the main plot using the buttons on the top right, but make sure to de-select them all before starting to add a unit (they interfere with the roi drawing)
   - when adding units or ROIs, double-click into each roi after drawing to end the drawing
   - info in unit table: unit: unit number and color; ch: main channel at which unit is detected; vis: unit visible or invisible in plot; cat: su, mu or noise; ISIv: percentage of events with ISIs below 1.2ms 
   - waveform plot: choose up to 2 units to display; can add additional waveforms for comparison; nr spikes controls how many spikes are shown per unit
   - unit footprint: distribution of detection channels for a selected unit
   - output is spkSort structure, which in addition to general info contains *unitid*: vector with unit assignment for each time stamp (-1 artefact, no distinction between SU, MU and noise), *spktimes*: vector with all spike time stamps, *unitinfo*: cell array with type assignemt (SU, MU, noise, none) for each unitid
   *Output: filename_pX_spkSort.mat (local and on Z); updates to _id.mat file (local and on Z) and database*

5) for unsorted (MUA) analysis - mergeSpkInfo:
   - merges the individual spkInfo files into one file for MUA analysis
   - keeps spktimes, detCh, detChSort in structure spkMerge 
   *Output: filename_pX_spkMerge.mat*

## Further processing:
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



