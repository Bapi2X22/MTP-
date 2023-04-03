# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include "tgefuncs.h"

char INFITS[128], input[128], GV[128];


int main(int argc, char *argv[])
{

  if(argc!=4) 
    {
      printf("Usage: %s <input FITS file> <input parameters file> <output GV file>\n", argv[0]);
      return 1;
    }
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",input); 
  sscanf(argv[3],"%s",GV); 
  
  read_inputs(input); 
  initialize(); 

  group_loop();
  WriteGVnAC();
  printf("\nEnd.\n");
  return 0;   
} 
