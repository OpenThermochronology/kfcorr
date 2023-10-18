---
### ABOUT kfcorr

Calculation of correlation coefficient between age and log1o(r/ro) spectra.

Version of Oscar Lovera's corrfft code. Changes made by Peter Zeitler to permit simpler input and output, and to plot results if the user has the gmt package installed.

---
### AUTHOR

Peter Zeitler, Lehigh University, Bethlehem, PA USA

---
### COMPILATION

Compile `kfcorr` like this (you MUST do it this way, using the `fno-automatic` and the `fallow-argument-mismatch` flags):

`gfortran kfcorr210.f90 -o kfcorrm2 -fno-automatic -O2  -fallow-argument-mismatch -w`

(substitute your local source file name for `kfcorr210.f90`, and your preferred executable name for `kfcorrm2`).

Using O2 optimization seems to work reliably. The two '-f...' flags are required to make this complex legacy code compile - newer compilers are more strict about F77
code that used to generate warnings. These now throw errors that block compilation, so the two flags are required to work around that.

For MacOS, a good source for gcc and gfortran installer packages can be found at:

[hpc.sourceforge.net](https://hpc.sourceforge.net)

---
### REQUIREMENTS FOR PLOTTING

Although kfcorr can be used without the plotting option, it's useful to a have a visual repord.

For the plotting option to work, you need to have an installation of either gmt 5 or gmt 6. Installation packages and instructions can be found at:

[www.generic-mapping-tools.org/download/](https://www.generic-mapping-tools.org/download/)

---
### INPUT and OUTPUT

#### Input

Two files are required, both with UNIX-style line breaks and both having the same length. The files should have as many entries as there are heating steps in the age spectrum. For both files, the first column is cumulative fractional loss, and the second column is either age or log10(r/r<sub>o</sub>).

#### Output

Text output is minimal - it's just the cross-correlation coefficient, reported to the console. If you requested the plot option, you should also get a gmt-generated pdf plot having the name `samplename-correlation.pdf`. On this plot, the age spectrum is red and the log10(r/r<sub>o</sub> spectrum is blue. Regions that are not included in the correlation are grayed out. The correlation over the specified region is also printed on the plot.

---
### USAGE

`./kfcorrm2 spectrumfilename  RRofilename  plotflag  samplename`

- `spectrumfilename` refers to a file of format fractional loss and age 

- `RRofilename` refers to a file of format fractional loss and log10(r/r<sub>o</sub>)

- `plotflag` is an integer: 1 results in use of gmt plotting commands to display output

- `samplename` is used to tag the output filename. This can be omitted for a generic result

*NOTE: You can place place a copy of your kfcorr executable in any directory in your PATH, and then do work with this code in any other directory containing data (i.e., you won't need to drag around copies of the executable). In this case, in any directory containing your datafile, USAGE becomes simply:*

`kfcorrm2 spectrumfilename  RRofilename  plotflag  samplename`


---
### HELP and EXAMPLES

There is no help. You are alone in a vast and uncaring universe (but you do have the Open Thermochronology GitHub site as a friend).

There is an EXAMPLES directory including two pairs of input files and examples of the output plots for them.

---