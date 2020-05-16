A.  	Instructions for Initial Setup.

For new users:
You must select the correct head stage developer in the file get_headstage_source.m located in the sub-folder microprobe wiring/headstage and pcb wiring.

Selection only needs to take place once per computer; no further changes are necessary after that.


************************************************************
Instructions:

1.	Open get_headstage_source.m in the Matlab editor.

2a.	If you are using an Intan Technologies headstage, type headstage_source='Intan';
2b.	If you are using a UCLA head stage, type headstage_source='UCLA';

3.  	Save the file.  No further changes are necessary.

*************************************************************


B.	Contents of Microprobe Wiring Folder

The microprobe wiring folder contains Matlab executable scripts that allow users to determine the relative Cartesian coordinates of every recording site.

Running a "probe_XYZ" script (e.g., probe_64D, probe_256A) will generate an array called s, containing the following information:

s.channels	lists the amplifier channel numbers. 
s.x		lists the relative x coordinate of the recording site in units of microns.
s.y		lists the relative y coordinate of the recording site in units of microns.  The default value is zero unless the probe is 3D.
s.z		lists the relative z coordinate of the recording site in units of microns.  The default is for the deepest electrode to have z=0.
s.shaft		lists which silicon prong the recording site is on.  For single-prong probes this value is 1, for multi-prong probesthis value is 1, 2, 3...
s.tipelectrode	lists the distance between the tip of the probe and the edge of the nearest electrode, in units of microns.


Additionally, running a script will display the amplifier channels numbers at their respective location on the microprobe array.


C. 	If your Molex Slimstack connector is attached on the bottom of the printed circuit board (usually only if you needed to attach an optical fiber), then
select the corresponding script ending in "_bottom".  The default is for the connector to be attached on the top.


D.	Questions and bugs: please contact Sotiris Masmanidis (smasmanidis at ucla dot edu).