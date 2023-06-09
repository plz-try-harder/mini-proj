// implementation from STM DT0089 Design Tip “The Goertzel algorithm to
// compute individual terms of the discrete Fourier transform (DFT)”

#include <dtmf_detection.h>

#include <math.h>
#include <stdint.h>

#define N	 200	// points for Goertzel
#define Fs	8000	// audio sampling freq.

#define b	12		// scaling bits
#define th	3000000	// detection threshold - nneds tuning

static const float S = (1<<b);	// scaling factor

// DTMF frequencies
static const int frow[4] = { 697, 770, 852, 941 }; // 1st tone
static const int fcol[4] = { 1209, 1336, 1477, 1633 }; // 2nd tone

// DTMF symbols
static const char sym[16] = {
		'1', '4', '7', '*',
		'2', '5', '8', '0',
		'3', '6', '9', '#',
		'A', 'B', 'C', 'D'
};

// DTMF decoding matrix
static const int symmtx[4][4] = {
		{ 0, 4, 8, 12 },
		{ 1, 5, 9, 13 },
		{ 2, 6, 10, 14 },
		{ 3, 7, 11, 15 }
};

static int win[N]; // Window

// Goertzel

// we detect 8 frequencies: 4 for column + 4 for row

static int c[8], cw[8], sw[8]; 	// Goertzel constants
static int z1[8], z2[8]; 		// Goertzel status registers
static int I[8], Q[8], M2[8]; 	// Goertzel output: real, imag, squared magnitude


void dtmf_detection_init()
{
	for(int i=0;i<N;i++) { // init window (Hamming)
		win[i] = (int)round(S*(0.54 -0.46*cosf(2.0*M_PI*(float)i/(float)(N-1))));
	}

	for(int i=0;i<4;i++) { // init Goertzel constants

		float u = 2.0*M_PI*round((float)N*(float)frow[i]/(float)Fs)/(float)N;
		cw[i] = (int)round(S*cosf(u));
		 c[i] = cw[i]<<1;
		sw[i] = (int)round(S*sinf(u));

		float v = 2.0*M_PI*round((float)N*(float)fcol[i]/(float)Fs)/(float)N;
		cw[i+4] = (int)round(S*cosf(v));
		 c[i+4] = cw[i+4]<<1;
		sw[i+4] = (int)round(S*sinf(v));
	}
}

// len must be at least N. Only the first N samples are processed
char dtmf_detection_detect(uint16_t buffer[], int len)
{
	if(len<N) {
		return -1;
	}

	// Goertzel reset
	for(int i=0;i<8;i++) {
		z1[i]=0;
		z2[i]=0;
	}

	for(int n=0; n<N; ++n ) {

		int x = ((buffer[n]*win[n])>>b); // windowing

		for(int i=0; i<8; i++) {
			float z0 = x + ((c[i]*z1[i])>>b) - z2[i]; // Goertzel iteration
			z2[i] = z1[i]; 	// Goertzel status update
			z1[i] = z0;
		}
	}

	// finalize and decode
	int i1 = -1;
	int i2 = -1;
	for(int i=0;i<8;i++) {
		I[i] = ((cw[i]*z1[i])>>b) - z2[i]; 	// Goertzel final I
		Q[i] = ((sw[i]*z1[i])>>b);			// Goertzel final Q
		M2[i] = I[i]*I[i] + Q[i]*Q[i]; 		// magnitude squared

		// DTMF decoding
		if(M2[i]>th) {
			if(i<4) {
				// find 1st tone, one peak allowed
				if(i1==-1)
					i1=i;
				else
					i1=4;
			} else {
				// find 2nd tone, one peak allowed
				if(i2==-1)
					i2=i-4;
				else
					i2=4;
			}
		}
	}
	if((i1>-1)&&(i1<4)&&(i2>-1)&&(i2<4))
		return sym[symmtx[i1][i2]];
	else
		return ' '; // nothing detected
}
