FITSDIR=/home/18ts91r02/mylib

visfitsgrid: fitsprog.c visfitsgrid.c read_fits_func.c 
	gcc -g -o visfitsgrid fitsprog.c visfitsgrid.c read_fits_func.c -I$(FITSDIR)/include -L$(FITSDIR)/lib -lcfitsio -lfftw3 -lm
	rm -rf *~

clean:
	rm -rf *~
	rm -rf *.o

