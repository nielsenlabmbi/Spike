# Spike

## Spike processing pipeline:

1) thresholdGui (GUI) and thresholdGuiSettings (GUI):
   - generate settings file (saved locally) if necessary 
   - visually set threshold for each channel, mark bad channels
   - needs microprobe_wiring directory and read_Intan_Header
   *Output: filename_threshold.mat and filename_id.mat; generates SpikeFiles folder*

2) extractSpikes and extractSpikesSettings (GUI):
   - first generate settings file (saved locally) with general settings using GUI
   - then run extractSpikes to extract waveforms for each channel (use parfor for speed)
   - implements temporal (and optionally spatial) constraints
   - note: first job should be 0, not 1
   *Output: SpikeFiles\filename_jID_spike.mat; updates to _id.mat file*

3) extractSpikeProps:
   - extract set of properties for each waveform
   - runs on individual job files (use parfor for speed)
   - splits data according to probes
   *Output: SpikeFiles\filename_jID_pX_spkinfo.mat; updates to _id.mat file*

4) mergeSpikeProps:
   - merges individual property files
   - works on data from 1 probe
   *Output: filename_spk_pX.mat; updates to _id.mat file*

5) sortGUI



