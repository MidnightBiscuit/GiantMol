# GiantMol
Repository for GiantMol data analysis programs. Created and maintained by AP since 2020 02 27.

## Programs

[TOC]

### FortranSourceCode

The Fortran Source Codes used for simulations in different contexts. There is code for simulating the molecular detection of a GiantMolecule passing through a cooled ion cloud. There is also two codes for simulating the behaviour of the ion cloud when submitted only to RF fields (no laser) : one cuts the laser just after laser cooling, one cuts the laser after the injection phase.

### SimAnalysisPointbyPoint

This program treats data from interaction simulations of GMol with Ca+. The simulation produces one directory for each condition (i.e point with unique paramters). This structured is used since 2020 01 31.


### PlotOscilloLabview

This program allows for plotting data taken from oscilloscopes or labview. I set parameters for different oscillo brands and versions (Siglent 1102CML+, Lecroy Waveace 2022, Lecroy Wavesurfer 3034). I talk about GiantMol Labview. This is useful for any oscillo-based purpose, in particular I need this for MCP reading.



### H2MAnalysis

This program treats data from H2M mass spectrometry. I can test molecules on their setup and then process the MCP measurment with this program.



### MCPAnalysis

This program treats data from MCP. Often recorded with oscilloscopes.



### ImageAnalysis

This program treats images from camera. Can also be useful for any other image analysis, such as beam profile.



