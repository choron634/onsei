#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FFT (double *xr, double *xi, double *Xr, double *Xi, int N);

void FFT (double *xr, double *xi, double *Xr, double *Xi, int N)
{
  int i,j,k,n,n2;
  double theta, wr, wi;

  static double *rbuf, *ibuf;
  static int bufsize = 0;

  /* memory allocation for buffers */
  if ( bufsize != N )
    {
      bufsize = N;
      rbuf = (double*)calloc(  bufsize, sizeof(double) );
      ibuf = (double*)calloc(  bufsize, sizeof(double)  );
    }

  /* bit reverse of xr[] & xi[] --> store to rbuf[] and ibuf[] */
  i = j = 0 ;
  rbuf[j] = xr[j];  ibuf[j] = xi[j];
  for( j = 1 ; j < N-1 ; j++ )
    {
      for( k = N/2 ; k <= i ; k /= 2 )  i -= k;
      i += k;
      rbuf[j] = xr[i];  ibuf[j] = xi[i];
    }
  rbuf[j] = xr[j];  ibuf[j] = xi[j];

   /* butterfly calculation */
  theta = -2.0*M_PI;
  for( n = 1 ; ( n2 = n*2 ) <= N ; n = n2 )
    {
      theta *= 0.5;
      for ( i = 0 ; i < n ; i++ )
	{
	  wr = cos(theta*i);  wi = sin(theta*i);
	  for ( j = i ; j < N ; j += n2 )
	    {
	      k = j + n;
	      /*
		X[j] = 1*buf[j] + W*buf[k];
		X[k] = 1*buf[j] - W*buf[k];
		Note : X[], buf[], and W are complex numbers.
		Re{ X[n] } = Xr[n], Im{ X[n] } = Xi[n];
		Re{ buf[n] } = rbuf[n], Im{ buf[n] } = ibuf[n];
		Re{ W } = wr, Im{ W } = wi;
	      */
	      Xr[j] = rbuf[j] + wr*rbuf[k] - wi*ibuf[k];  /* ??????????, using wr, wi, rbuf, and ibuf */
	      Xi[j] = ibuf[j] + wr*ibuf[k] + wi*rbuf[k];  /* ??????????, using wr, wi, rbuf, and ibuf */
	      Xr[k] = rbuf[j] - wr*rbuf[k] + wi*ibuf[k];  /* ??????????, using wr, wi, rbuf, and ibuf */
	      Xi[k] = ibuf[j] - wr*ibuf[k] - wi*rbuf[k];  /* ??????????, using wr, wi, rbuf, and ibuf */
	    }
	}
      for( i = 0 ; i < N ; i++ )
	{
	  rbuf[i] = Xr[i];
	  ibuf[i] = Xi[i];
	}
    }
  return;
}

int main (int argc, char **argv)
{
  /* check the format of input */
  if (argc != 4)
    {
      fprintf(stderr, "Usage: %s DATfile skip[sample] frame_length[sample]\n", argv[0]);
      exit(1);
    }
  FILE* fpDAT;
  int nskip;
  int framelen;

  /* check the validity of input */
  if( ( fpDAT = fopen( argv[1], "rb" ) ) == NULL )  exit( 1 );
  if( ( nskip    = atoi( argv[2] ) ) < 0 )  exit( 1 );
  if( ( framelen = atoi( argv[3] ) ) < 0 )  exit( 1 );

  fprintf( stderr, "# DATfile = %s\n", argv[1] );
  fprintf( stderr, "# %d samples are skipped.\n", nskip );
  fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );

  short *sdata = (short*)calloc( framelen, sizeof(short) );
  double *correlation = (double*)calloc(framelen,sizeof(double));
  double *xr = (double*)calloc(framelen,sizeof(double));
  if ( sdata == NULL) exit (1);

  fseek( fpDAT, nskip*sizeof(short), SEEK_SET );
  fread( sdata, sizeof(short), framelen, fpDAT);
  fclose( fpDAT );

  /* windowing */
  for (int i = 0 ; i < framelen; i++)
    {
      xr[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i];
    }


  /* cauculate correlation */
  for (int r = 0; r < framelen; r++){
    for (int t = 0; t < framelen-r ;t++){
      //correlation[r] += sdata[t]*sdata[t+r];
      correlation[r] += xr[t]*xr[t+r];
    }
    for(int t = framelen-r; t < framelen; t++)
    {
      correlation[r] += xr[t]*xr[t-framelen+r];
    }
  }

  /* normalize */
  double r0=correlation[0];
  for (int i = 0; i < framelen; i++){
    correlation[i]=correlation[i]/r0;
  }

  double *rcorrelation = (double*)calloc( framelen, sizeof(double) );
  double *icorrelation = (double*)calloc( framelen, sizeof(double) );
  double *Sr = (double*)calloc( framelen, sizeof(double) );
  double *Si = (double*)calloc( framelen, sizeof(double) );

  for(int i = 0; i < framelen; i++)
  {
    rcorrelation[i] = correlation[i];
  }

  FFT(rcorrelation, icorrelation, Sr, Si, framelen);

  for (int i = 0 ; i < framelen ; i++)
    {
      double spec = log10(Sr[i]);
      printf("%lf\n",spec);
    }

  return 0;
}
