#ifndef mcufft_h
#define mcufft_h

#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif


/* Custom constants */
#define FFT_FORWARD 0x01
#define FFT_REVERSE 0x00

/* Windowing type */
#define FFT_WIN_TYP_RECTANGLE 			0x00  /* rectangle (Box car) */
#define FFT_WIN_TYP_HAMMING 			0x01  /* hamming */
#define FFT_WIN_TYP_HANN 				0x02  /* hann */
#define FFT_WIN_TYP_TRIANGLE 			0x03  /* triangle (Bartlett) */
#define FFT_WIN_TYP_NUTTALL 			0x04  /* nuttall */
#define FFT_WIN_TYP_BLACKMAN 			0x05  /* blackman */
#define FFT_WIN_TYP_BLACKMAN_NUTTALL 	0x06  /* blackman nuttall */
#define FFT_WIN_TYP_BLACKMAN_HARRIS 	0x07  /* blackman harris*/
#define FFT_WIN_TYP_FLT_TOP 			0x08  /* flat top */
#define FFT_WIN_TYP_WELCH 				0x09  /* welch */


void mcufft_init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency);
void mcufft_windowing(uint8_t windowType, uint8_t dir);
void mcufft_compute(uint8_t dir);
void mcufft_complex_to_magnitude(void);
void mcufft_major_peak(float* mag_out, float* freq_out, float magFact);


uint8_t	mcufft_lib_revision(void);


#ifdef __cplusplus
}
#endif


#endif  // mcufft_h
