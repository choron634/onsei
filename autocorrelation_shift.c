#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char **argv)
{
  /* check the format of input */
  if (argc != 5)
  {
    fprintf(stderr, "Usage: %s DATfile skip[sample] frame_length[sample] outputfile\n", argv[0]);
    exit(1);
  }
  FILE* fpDAT;
  int nskip;
  int framelen;
  int count = 0;
  while(1){
    /* check the validity of input */
    if( ( fpDAT = fopen( argv[1], "rb" ) ) == NULL )  exit( 1 );
    if( ( nskip    = atoi( argv[2] ) ) < 0 )  exit( 1 );
    if( ( framelen = atoi( argv[3] ) ) < 0 )  exit( 1 );

    // fprintf( stderr, "# DATfile = %s\n", argv[1] );
    // fprintf( stderr, "# %d samples are skipped.\n", nskip + count*framelen/2);
    // fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );


    short *sdata = (short*)calloc( framelen, sizeof(short) );
    double *correlation = (double*)calloc(framelen,sizeof(double));
    double *xr = (double*)calloc(framelen,sizeof(double));
    if ( sdata == NULL) exit (1);

    fseek( fpDAT, (nskip + count * framelen/2)*sizeof(short), SEEK_SET );
    int size = fread( sdata, sizeof(short), framelen, fpDAT);
    if(size == 0){
      break;
    }
    fclose( fpDAT );

    /* windowing */
    for (int i = 0 ; i < framelen; i++)
    {
      xr[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i];
    }


    /* criculate correlation */
    for (int r = 0; r < framelen; r++){
      for (int t = 0; t < framelen-1-r ;t++){
        //correlation[r] += sdata[t]*sdata[t+r];
        correlation[r] += xr[t]*xr[t+r];
      }
    }

    /* normalize */
    double r0=correlation[0];
    double max = 0;
    int max_index = 0;
    for (int i = 10; i < framelen; i++){
      correlation[i]=correlation[i]/r0;
      if(max < correlation[i] && i != 0)
      {
        max = correlation[i];
        max_index = i;
      }
    }

    int ff = max_index;

    FILE *fp;
    if( ( fp = fopen( argv[4], "a" ) ) == NULL )  exit( 1 );

    fprintf(fp, "%lf\n", max);

    fclose(fp);

    count++;
  }
  return 0;
}
