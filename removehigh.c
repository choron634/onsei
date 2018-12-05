#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
  FILE* fpDAT;
  int threshold;
  int i;
  int framelen;


  if( ( fpDAT = fopen( argv[1], "rb" ) ) == NULL )  exit( 1 );
  if( ( threshold    = atoi( argv[2] ) ) < 0 )  exit( 1 );
  if( ( framelen = atoi( argv[3] ) ) < 0 )  exit( 1 );

  fprintf( stderr, "# DATfile = %s\n", argv[1] );
  fprintf( stderr, "# threshold is %d.\n", threshold );
  fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );

  short *sdata = (short*)calloc( framelen, sizeof(short) );

  if ( sdata == NULL) exit (1);

  fread( sdata, sizeof(short), framelen, fpDAT);
  fclose(fpDAT);

  for ( i = threshold ; i < framelen; i++)
    {
      sdata[i] = 0;
    }

  for(i = 0; i < framelen; i++){
    printf( "%lf\n", sdata[i] );
  }
}
