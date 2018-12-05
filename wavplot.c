#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){
  FILE* fpDAT;
  int nskip;
  int i;
  int framelen;


  if( ( fpDAT = fopen( argv[1], "rb" ) ) == NULL )  exit( 1 );
  if( ( nskip    = atoi( argv[2] ) ) < 0 )  exit( 1 );
  if( ( framelen = atoi( argv[3] ) ) < 0 )  exit( 1 );

  fprintf( stderr, "# DATfile = %s\n", argv[1] );
  fprintf( stderr, "# %d samples are skipped.\n", nskip );
  fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );

  short *sdata = (short*)calloc( framelen, sizeof(short) );

  if ( sdata == NULL) exit (1);

  fseek( fpDAT, nskip*sizeof(short), SEEK_SET );
  fread( sdata, sizeof(short), framelen, fpDAT);
  fclose(fpDAT);

  // for ( i = 0 ; i < framelen; i++)
  //   {
  //     sdata[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i];
  //   }

  for(i = 0; i < framelen; i++){
    fprintf(stdout,"%d\n",sdata[i]);
  }
}
