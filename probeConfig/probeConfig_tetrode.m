%generate wiring information for tetrode
%standard format for columns in probewiring:
%column 1: channel number
%column 2: x
%column 3: y
%column 4: z
%column 5: shaft
function probewiring=probeConfig_tetrode

probewiring=[1 -20 0 -20 1
    2 -20 0 20 1
    3 20 0 -20 1
    4 20 0 20 1];