%generate wiring information for tetrode
%standard format for columns in probewiring:
%column 1: channel number (starts with 0)
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_tetrode

probewiring=[0 -10 0 -10 1
    1 -10 0 10 1
    2 10 0 -10 1
    3 10 0 10 1];