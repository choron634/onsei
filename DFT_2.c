#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void DFT( double *xr, double *xi, double *Xr, double *Xi, int N )
{
  int k, n;
  double wr, wi;

  for( k = 0 ; k < N ; k++ ) {
    Xr[k] = Xi[k] = 0.0;
    for( n = 0 ; n < N ; n++ ) {
      wr = cos( 2.0*M_PI*n*k/N );  wi = sin( 2.0*M_PI*n*k/N );
      Xr[k] += xr[n]*wr + xi[n]*wi; /* should be modified. (hint using wr & wi) */
      Xi[k] += -xr[n]*wi + xi[n]*wr; /* should be modified. (hint using wr & wi) */
    }
  }
  return;
}

void IDFT( double *Xr, double *Xi, double *xr, double *xi, int N )
{
  int k, n;
  double wr, wi;

  for( n = 0 ; n < N ; n++ ) {
    xr[n] = xi[n] = 0.0;
    for( k = 0 ; k < N ; k++ ) {
      wr = cos( 2.0*M_PI*n*k/N );  wi = sin( 2.0*M_PI*n*k/N );
      xr[n] += Xr[n]*wr - Xi[n]*wi; /* should be modified. (hint using wr & wi) */
  	  xi[n] += Xr[n]*wi + Xi[n]*wr; /* should be modified. (hint using wr & wi) */
    }
    xr[n] /= N;  xi[n] /= N;
  }
  return;
}

int main( int argc, char *argv[] )
{
  int  i, nskip, framelen;
  double spec;
  short  *sdata;
  double *xr, *xi, *Xr, *Xi;
  FILE *fpDAT;

  if( 4 != argc ) {
    fprintf( stderr, "Usage : %s DATfile skip[sample] frame_length[sample]\n", argv[0] );
    exit( 1 );
  }

  if( NULL == ( fpDAT = fopen( argv[1], "r" ) ) )  exit( 1 );
  if( 0 > ( nskip    = atoi( argv[2] ) ) )  exit( 1 );
  if( 0 > ( framelen = atoi( argv[3] ) ) )  exit( 1 );

  fprintf( stderr, "# DATfile = %s\n", argv[1] );
  fprintf( stderr, "# %d samples are skipped.\n", nskip );
  fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );

  /* memory allocation */
  sdata = (short*)calloc( sizeof(short), framelen );
  xr = (double*)calloc( sizeof(double), framelen );
  xi = (double*)calloc( sizeof(double), framelen );
  Xr = (double*)calloc( sizeof(double), framelen );
  Xi = (double*)calloc( sizeof(double), framelen );
  if( NULL == sdata || NULL == xi || NULL == xi || NULL == Xr || NULL == Xi )  exit( 1 );

  fseek( fpDAT, nskip*sizeof(short), SEEK_SET );
  fread( sdata, sizeof(short), framelen, fpDAT );
  fclose( fpDAT );

  for( i = 0 ; i < framelen ; i++ ) {
    xr[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1)))*sdata[i];
    xi[i] = 0.0;
  }

  /* for ( i = 0 ; i < 1000 ; i++ ) { */
  DFT( xr, xi, Xr, Xi, framelen );
  /* } */

  /*
  IDFT( Xr, Xi, xr, xi, framelen );
  */

  for( i = 0 ; i < framelen ; i++ ) {
    spec = log10( (1.0/framelen)*( Xr[i]*Xr[i]+Xi[i]*Xi[i] ) );
    printf( "%d %lf\n", i, spec );
  }

  /*
  for( i = 0 ; i < framelen ; i++ ) {
    printf( "%d %lf\n", i, xr[i] );
  }
  */

  return( 0 );
}
