libs = -L ~/Software/cfitsio/ \
-L/Users/gullo/Software/heasoft-6.28/x86_64-apple-darwin18.7.0/lib \
-lcfitsio -lXSFunctions -lXSModel -lXSUtil -lXS #-lCCfits_2.5

flags = -Wall

line_limit = -ffree-line-length-none

checks= -pedantic -Wall

# main = check_table.f90
main = check_table.f90
main_single_spectrum = plot_single_table_spectrum.f90
routines = read_fits.f90 compton_hump_estimate.f90 donthcomp.f xsdskb.f 

compile:
	gfortran $(libs) $(line_limit) -c $(routines) $(main)

comp_single:
	gfortran $(libs) $(line_limit) -c $(routines) $(main_single_spectrum)

run:  compile
	rm -rf read
	gfortran  $(flags) $(libs) -o read *.o
	rm -rf *.o *~

single: comp_single
	rm -rf single_spectrum
	gfortran $(flags) $(libs) -o single_spectrum *.o 
	rm -rf *.o *~


clean:
	rm -rf *.o *~ read single_spectrum fort.* 

