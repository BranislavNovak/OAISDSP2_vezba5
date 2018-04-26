#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "obrada.h"
#include "mdct.h"

#define B 12	// 10, 8, 6

double time_buffer[MDCT_SIZE];
double mdct_buffer[MDCT_SIZE/2];
double in_delay[MDCT_SIZE/2];
double out_delay[MDCT_SIZE/2];
int histogram[1<<B];
double p[1 << B];
double b[1 << B];

extern double window[MDCT_SIZE];

void obrada(double *in, double *out, int N)
{
  int i;

  for (i=0; i<N; i++)
    out[i] = in[i];

  // zad1
  for (int i = 0; i < N; i++)
  {
	  time_buffer[i] = in_delay[i];
  }

  for (int i = 0; i < N; i++)
  {
	  time_buffer[N + i] = in[i];
  }

  mdct(time_buffer, mdct_buffer);

  // zad2		 koder	\/
  for (i = 0; i < N; i++)
  {
	  if (mdct_buffer[i] < 0) {
		  mdct_buffer[i] = -1 * sqrt((-1)*mdct_buffer[i]);
	  }
	  else {
		  mdct_buffer[i] = sqrt(mdct_buffer[i]);
	  }
  }

  // zad3		 kvantizacija	\/
  for (i = 0; i < N; i++)
  {
	  //mdct_buffer[i] = round((mdct_buffer[i] / (1 << 12))*(1 << (B - 1)));	// kvantizacija iz [-2^12,2^12) u opseg [-2^B-1, 2^B-1)
	  mdct_buffer[i] = round((mdct_buffer[i] / (1 << (13 - B))));			// dva deljenja stavljena u jedno
  }


  // zad4		entropijsko kodovanje	   \/
  for (i = 0; i < N; i++) 
  {
	  histogram[(int)mdct_buffer[i] + (1 << (B - 1))]++;				// posto bi niz isao u negativno, nemamo negativne indekse, pa onda pomeramo u desno niz za 2^b-1
  }
  
  // zad4		entropijsko kodovanje	   /\

  for (i = 0; i < N; i++) 
  {
	  //mdct_buffer[i] = round((mdct_buffer[i] * (1 << 12))*(1 << (B - 1)));	// inverzna kvantizacija iz  [-2^B-1, 2^B-1) u opseg [-2^12,2^12)
	  mdct_buffer[i] = round(mdct_buffer[i] * (1 << (13 - B)));				// dva mnozenja stavljeno u jedno
  }

  // zad3		 inverzna kvantizacija		/\

  for (i = 0; i < N; i++)
  {
	  if (mdct_buffer[i] < 0) {
		  mdct_buffer[i] = (-1)*mdct_buffer[i] * mdct_buffer[i];
	  }
	  else {
		  mdct_buffer[i] = mdct_buffer[i] * mdct_buffer[i];
	  }
  }
  // zad2		  dekoder  /\
  
  imdct(mdct_buffer, time_buffer);

  for (int i = 0; i < N; i++)
  {
	  out[i] = (time_buffer[i] + out_delay[i]) / 2;
  }

  for (int i = 0; i < N; i++) 
  {
	  in_delay[i] = in[i];
  }

  for (int i = 0; i < N; i++)
  {
	  out_delay[i] = time_buffer[N + i];
  }
  // zad1
}

void statistika()
{
	int suma_histograma = 0;
	int i;
	int broj_bita = 0;

	for (i = 0; i < (1 << B); i++) {
		suma_histograma += histogram[i];
	}

	for (i = 0; i < (1 << B); i++) {
		p[i] = (double)histogram[i] / suma_histograma;
	}

	for (i = 0; i < (1 << B); i++) {
		if (p[i] == 0) {
			b[i] = 0;
		}
		else {
			b[i] = (double)log2(1 / p[i]);
		}
	}

	for (i = 0; i < (1 << B); i++) {
		broj_bita += histogram[i] * b[i];
	}

	printf("Broj bita: %d\n", broj_bita);
	printf("Broj mega bita: %f\n", broj_bita*0.000000125);
	printf("Stepen kompresije: %f\n", ((16.0 * suma_histograma)/broj_bita));
}

