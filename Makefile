libs = -L ~/Software/cfitsio/ \
-L/Users/gullo/Software/heasoft-6.28/x86_64-apple-darwin18.7.0/lib \
-lcfitsio -lXSFunctions -lXSModel -lXSUtil -lXS -lCCfits_2.5

flags = -Wall

line_limit = -ffree-line-length-none

checks= -pedantic -Wall

main = check_table.f90
routines = read_fits.f90

compile:
	gfortran $(libs) $(line_limit) -c $(routines) $(main)

run:  compile 
	rm -rf read
	gfortran  $(flags) $(libs) -o read *.o

# go:
# 	gfortran $(main) $(flags) $(libs) -o read 


clean:
	rm -rf *.o *~ read
