FITSDIR=/home/18ts91r02/mylib
target: ngrf_2d clean

ngrf_2d:  ngrf_2d.c fitsprog_2d.c ngrf_2d_funcs.c beam.c
	gcc -g -o ngrf_2d   ngrf_2d.c ngrf_2d_funcs.c fitsprog_2d.c beam.c -I$(FITSDIR)/include -L$(FITSDIR)/lib -lfftw3 -lcfitsio -lm -lgsl -lgslcblas
	rm -rf *~ *.o
	
read_ngrf_3d:  read_ngrf_3d.c fitsprog_2d.c beam.c
	gcc -g -o read_ngrf_3d   read_ngrf_3d.c  fitsprog_2d.c beam.c -I$(FITSDIR)/include -L$(FITSDIR)/lib -lfftw3 -lcfitsio -lm -lgsl -lgslcblas
	rm -rf *~ *.o
clean:
	rm -rf ngrf_2d
