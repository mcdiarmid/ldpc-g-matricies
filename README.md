# IRIG 106 LDPC Generator Matricies

This repo exists solely for anyone interested in creating 
and using the LDPC generator matricies outlined in the 
[IRIG 106 Telemetry Standards](http://www.irig106.org/docs/106-20/106-20_Telemetry_Standards.pdf)
pdf file.  LDPC matricies are commonly stored in the [alist file format](http://www.inference.org.uk/mackay/codes/alist.html).
For example, [GNURadio](https://github.com/gnuradio/gnuradio) LDPC encoders and decoders utilize this file format.
Therefore, part of the script is a function for writing these generator matricies to an alist file.
For convenience sake, the six alist files for the combinations of 
k=1024,4096 and r=1/2,3/4,4/5 exist in the [alist folder](./alist)
for anyone who just wants to download these and plug them into an LDPC encoder or decoder.
