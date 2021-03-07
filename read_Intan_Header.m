function header=read_Intan_Header(filename)

% based on read_Intan_RHD2000_file
%
% Version 3.0, 8 February 2021
%
% reads only the header from an _info.rhd file
% currently does not return all information but only most relevant


fid = fopen(filename, 'r');


% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
header.magic_number = fread(fid, 1, 'uint32');
if header.magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
header.data_file_main_version_number = fread(fid, 1, 'int16');
header.data_file_secondary_version_number = fread(fid, 1, 'int16');

% Read information of sampling rate and amplifier frequency settings.
header.sample_rate = fread(fid, 1, 'single');
header.dsp_enabled = fread(fid, 1, 'int16');
header.actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
header.actual_lower_bandwidth = fread(fid, 1, 'single');
header.actual_upper_bandwidth = fread(fid, 1, 'single');
header.desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
header.desired_lower_bandwidth = fread(fid, 1, 'single');
header.desired_upper_bandwidth = fread(fid, 1, 'single');


% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
header.notch_filter_mode = fread(fid, 1, 'int16');

%impedance test info
desired_impedance_test_frequency = fread(fid, 1, 'single');
actual_impedance_test_frequency = fread(fid, 1, 'single');

%Read notes (if present)
note1=fread_QString(fid);
note2=fread_QString(fid);
note3=fread_QString(fid);


% If data file is from GUI v1.1 or later, see if temperature sensor data
% was saved.
num_temp_sensor_channels = 0;
if ((header.data_file_main_version_number == 1 && header.data_file_secondary_version_number >= 1) ...
    || (header.data_file_main_version_number > 1))
    num_temp_sensor_channels = fread(fid, 1, 'int16');
end

% If data file is from GUI v1.3 or later, load board mode.
board_mode = 0;
if ((header.data_file_main_version_number == 1 && header.data_file_secondary_version_number >= 3) ...
    || (header.data_file_main_version_number > 1))
    board_mode = fread(fid, 1, 'int16');
end

% If data file is from v2.0 or later (Intan Recording Controller),
% load name of digital reference channel.
if (header.data_file_main_version_number > 1)
    reference_channel = fread_QString(fid);
end


% Read signal summary from data file header.

header.number_of_signal_groups = fread(fid, 1, 'int16');

for g = 1:header.number_of_signal_groups
    header.signal_group_name{g} = fread_QString(fid);
    header.signal_group_prefix{g} = fread_QString(fid);
    header.signal_group_enabled(g) = fread(fid, 1, 'int16');
    header.signal_group_num_channels(g) = fread(fid, 1, 'int16');
    header.signal_group_num_amp_channels(g) = fread(fid, 1, 'int16');

    header.channel_signal_type{g}=[];
    header.channel_enabled{g}=[];
    
    if (header.signal_group_num_channels(g) > 0 && header.signal_group_enabled(g) > 0)
        
        %we're only saving certain info for now
        for s = 1:header.signal_group_num_channels(g)
            native_channel_name = fread_QString(fid);
            custom_channel_name = fread_QString(fid);
            native_order = fread(fid, 1, 'int16');
            custom_order = fread(fid, 1, 'int16');
            
            header.channel_signal_type{g}(s) = fread(fid, 1, 'int16'); %signifies signal type
            header.channel_enabled{g}(s) = fread(fid, 1, 'int16');
            
            chip_channel = fread(fid, 1, 'int16');
            board_stream = fread(fid, 1, 'int16');
            voltage_trigger_mode = fread(fid, 1, 'int16');
            voltage_threshold = fread(fid, 1, 'int16');
            digital_trigger_channel = fread(fid, 1, 'int16');
            digital_edge_polarity = fread(fid, 1, 'int16');
            electrode_impedance_magnitude = fread(fid, 1, 'single');
            electrode_impedance_phase = fread(fid, 1, 'single');           
        end
    end
    
    %summary data for channel
    N_amp_enabled=sum(header.channel_signal_type{g}==0 & header.channel_enabled{g}==1);
    header.signal_group_num_amp_enabled(g)=N_amp_enabled;
end

end


function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return
end
