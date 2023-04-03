FITSDIR=/home/asif/mylib
read_fits_funcmfs.o:	read_fits_funcmfs.c
	gcc -c -I$(FITSDIR)/include -L$(FITSDIR)/lib read_fits_funcmfs.c -o read_fits_funcmfs.o

tgefuncs.o:	tgefuncs.c
	gcc -c -I$(FITSDIR)/include -L$(FITSDIR)/lib tgefuncs.c -o tgefuncs.o
		

main.o:	main.c
	gcc -c -I$(FITSDIR)/include -L$(FITSDIR)/lib main.c -o main.o

tge:		read_fits_funcmfs.o tgefuncs.o main.o
	gcc read_fits_funcmfs.o tgefuncs.o main.o -I$(FITSDIR)/include -L$(FITSDIR)/lib -lm -lcfitsio -lfftw3 -lgsl -lgslcblas -o tge

clean:
	rm -rf *.o
	rm -rf *~

