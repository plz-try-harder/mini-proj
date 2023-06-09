#ifndef DTMF_DETECTION_H
#define DTMF_DETECTION_H

// implementation from STM DT0089 Design Tip “The Goertzel algorithm to
// compute individual terms of the discrete Fourier transform (DFT)”

#include <stdint.h>

// must be called first to init internal data structures
void dtmf_detection_init();

// input is buffer of unsigned 16-bit PCM values
// buffer length must be large enough (200) to run the algorithm
// return value is detected symbol, or ' ' (space for none) or <0 on error
char dtmf_detection_detect(uint16_t buffer[], int len);


#endif // #ifdef DTMF_DETECTION_H
