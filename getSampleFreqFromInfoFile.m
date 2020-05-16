function sample_rate = getSampleFreqFromInfoFile(filename)
    % read_Intan_RHD2000_file
    %
    % Version 2.01, 11 October 2017
    %
    % Reads Intan Technologies RHD2000 data file generated by evaluation board
    % GUI or Intan Recording Controller.  Data are parsed and placed into
    % variables that appear in the base MATLAB workspace.  Therefore, it is
    % recommended to execute a 'clear' command before running this program to
    % clear all other variables from the base workspace.
    %
    % Example:
    % >> clear
    % >> read_Intan_RHD2000_file
    % >> whos
    % >> amplifier_channels(1)
    % >> plot(t_amplifier, amplifier_data(1,:))

    % [file, path, filterindex] = ...
    %     uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'off');
    % 
    % if (file == 0)
    %     return;
    % end

    % Read most recent file automatically.
    % path = 'C:\Users\Reid\Documents\RHD2132\testing\';
    % d = dir([path '*.rhd']);
    % file = d(end).name;

    tic;
    fid2 = fopen(filename, 'r');

    % Check 'magic number' at beginning of file to make sure this is an Intan
    % Technologies RHD2000 data file.
    magic_number = fread(fid2, 1, 'uint32');
    if magic_number ~= hex2dec('c6912702')
        error('Unrecognized file type.');
    end

    % Read version number.
    data_file_main_version_number = fread(fid2, 1, 'int16');
    data_file_secondary_version_number = fread(fid2, 1, 'int16');

    % Read information of sampling rate and amplifier frequency settings.
    sample_rate = fread(fid2, 1, 'single');
end