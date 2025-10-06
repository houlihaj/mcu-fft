#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "mcufft.h"


/* Mathematial constants */
#define twoPi 	6.28318531
#define fourPi 	12.56637061
#define sixPi 	18.84955593

#define FFT_LIB_REV 0x10

/* Private Variables */
static uint16_t _samples;
static float 	_sampling_frequency;
static float*	_vReal;
static float*	_vImag;
static uint8_t 	_power;

/* Private Functions and Macros */
static void 	Swap(float* x, float* y);
static uint8_t 	Exponent(uint16_t value);
#define sq(x) 	((x)*(x))


/**
 * @brief Initialize InputSignal structure passed by pointer.
 *
 * @param self
 * @param arr_size
 * @return
 */
uint8_t InputSignal_create(InputSignal* self, const uint8_t arr_size)
{
    self->arr_size = arr_size;

    self->real_arr = malloc(arr_size * sizeof(*self->real_arr));
    self->imag_arr = malloc(arr_size * sizeof(*self->imag_arr));

    if (!self->real_arr || !self->imag_arr) {
        fprintf(stderr, "Error allocating memory for at least one array\n");
        free(self->real_arr);
        free(self->imag_arr);
        return 1;  // return error code
    }

    return 0;
}


/**
 * @brief Free memory allocated for InputSignal and its members.
 *
 * @param self
 * @return
 */
uint8_t InputSignal_destroy(InputSignal* self)
{
    free(self->real_arr);
    free(self->imag_arr);
    return 0;
}


/**
 * @brief Initialize InputSignal structure passed by pointer.
 *
 * @param self
 * @param fft_size
 * @return
 */
uint8_t FTSignal_create(InputSignal* self, const uint8_t fft_size)
{
    self->fft_size = fft_size;

    self->real_arr = malloc(arr_size * sizeof(*self->real_arr));
    self->imag_arr = malloc(arr_size * sizeof(*self->imag_arr));

    if (!self->real_arr || !self->imag_arr) {
        fprintf(stderr, "Error allocating memory for at least one array\n");
        free(self->real_arr);
        free(self->imag_arr);
        return 1;  // return error code
    }

    return 0;
}


/**
 * @brief Free memory allocated for InputSignal and its members.
 *
 * @param self
 * @return
 */
uint8_t FTSignal_destroy(InputSignal* self)
{
    free(self->real_arr);
    free(self->imag_arr);
    return 0;
}


/**
 * @brief MCU init
 *
 * @param real_arr
 * @param imag_arr
 * @param num_samples
 * @param sampling_frequency
 * @return
 */
uint8_t mcufft_init(
	float* real_arr, float* imag_arr, uint16_t num_samples, float sampling_frequency
)
{
	_vReal = vReal;
	_vImag = vImag;
	_samples = samples;
	_sampling_frequency = sampling_frequency;
	_power = Exponent(samples);
	return 0;
}


/**
 * @brief FFT compute
 *
 * @param dir
 * @return
 */
void mcufft_compute(uint8_t dir)
{// Computes in-place complex-to-complex FFT /
	// Reverse bits /
	uint16_t j = 0;
	for (uint16_t i = 0; i < (_samples - 1); i++) {
		if (i < j) {
			Swap(&_vReal[i], &_vReal[j]);
			if (dir==FFT_REVERSE)
				Swap(&_vImag[i], &_vImag[j]);
		}
		uint16_t k = (_samples >> 1);
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	// Compute the FFT  /
	float c1 = -1.0;
	float c2 = 0.0;
	uint16_t l2 = 1;
	for (uint8_t l = 0; (l < _power); l++) {
		uint16_t l1 = l2;
		l2 <<= 1;
		float u1 = 1.0;
		float u2 = 0.0;
		for (j = 0; j < l1; j++) {
			 for (uint16_t i = j; i < _samples; i += l2) {
					uint16_t i1 = i + l1;
					float t1 = u1 * _vReal[i1] - u2 * _vImag[i1];
					float t2 = u1 * _vImag[i1] + u2 * _vReal[i1];
					_vReal[i1] = _vReal[i] - t1;
					_vImag[i1] = _vImag[i] - t2;
					_vReal[i] += t1;
					_vImag[i] += t2;
			 }
			 float z = ((u1 * c1) - (u2 * c2));
			 u2 = ((u1 * c2) + (u2 * c1));
			 u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		c1 = sqrt((1.0 + c1) / 2.0);
		if (dir == FFT_FORWARD) {
			c2 = -c2;
		}
	}
	// Scaling for reverse transform /
	if (dir != FFT_FORWARD) {
		for (uint16_t i = 0; i < _samples; i++) {
			 _vReal[i] /= _samples;
			 _vImag[i] /= _samples;
		}
	}
}


/**
 * @brief Remove the DC term from the input signal
 *
 * @param signal
 * @param signal_size
 * @return
 */
uint8_t mcufft_dc_removal(float* signal, uint16_t signal_size)
{
	// Compute the mean
	float signal_sum = 0;
	for (uint16_t i = 0; i < signal_size; i++) {
		signal_sum += signal[i];
	}
	float signal_mean = signal_sum / signal_size;

	// Subtract the mean
	for (uint16_t i = 0; i < signal_size; i++) {
		signal[i] -= signal_mean;
	}

	return 0;
}


/**
 * @brief Compute the magnitude of the complex number
 *
 * @return
 */
uint8_t mcufft_complex_to_magnitude(void)
{ // vM is half the size of vReal and vImag
	for (uint16_t i = 0; i < _samples; i++) {
		_vReal[i] = sqrt(sq(_vReal[i]) + sq(_vImag[i]));
	}
	return 0;
}


/**
 * @brief Compute the magnitude of the complex number
 * using the alpha max plus beta min algorithm
 *
 * Pythagorean addition approximation using the
 * alpha max plus beta min algorithm
 *
 * Approximate the magnitude of a complex number z = a + b * i
 * given the real (a) and imaginary (b) parts.
 *
 * https://github.com/andrestumpf/coregis_public/blob/46479157cc8ab44299b1f61a31e7a575f76909a1/mpic/src/func.c#L1424
 *
 * @param real_arr
 * @param image_arr
 * @param arr_size
 * @return
 */
uint8_t mcufft_complex_to_magnitude_approximation(
	float* real_arr, float* image_arr, uint16_t arr_size
)
{
	// implementation of alpha max plus beta min algorithm
	// largest error: 3.96%, mean error: 2.41% 
	float min;
  float max;
  const float alpha = 0.96043387F;
  const float beta = 0.39782473F;

  for (uint16_t i = 0; i < arr_size; i++) {
    float re = real_arr[i];
    float im = imag_arr[i];

    if (re < 0) re = -re;
    if (im < 0) im = -im;
    if (re > im) {
			max = re;
			min = im;
    } else {
			min = re;
			max = im;
    }

    real_arr[i] = alpha * max + beta * min;
  }

  return 0;
}


/**
 * @brief Apply the window function
 *
 * @param windowType
 * @param dir
 * @return
 */
uint8_t mcufft_windowing(uint8_t windowType, uint8_t dir)
{// Weighing factors are computed once before multiple use of FFT
// The weighing function is symetric; half the weighs are recorded
	float samplesMinusOne = ((float)(_samples) - 1.0);
	for (uint16_t i = 0; i < (_samples >> 1); i++) {
		float indexMinusOne = (float)(i);
		float ratio = (indexMinusOne / samplesMinusOne);
		float weighingFactor = 1.0;
		// Compute and record weighting factor
		switch (windowType) {
		case FFT_WIN_TYP_RECTANGLE: // rectangle (box car)
			weighingFactor = 1.0;
			break;
		case FFT_WIN_TYP_HAMMING: // hamming
			weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_HANN: // hann
			weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_TRIANGLE: // triangle (Bartlett)
			weighingFactor = 1.0 - ((2.0 * abs(indexMinusOne - (samplesMinusOne / 2.0))) / samplesMinusOne);
			break;
		case FFT_WIN_TYP_NUTTALL: // nuttall
			weighingFactor = 0.355768 - (0.487396 * (cos(twoPi * ratio))) + (0.144232 * (cos(fourPi * ratio))) - (0.012604 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN: // blackman
			weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) + (0.07922 * (cos(fourPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_NUTTALL: // blackman nuttall
			weighingFactor = 0.3635819 - (0.4891775 * (cos(twoPi * ratio))) + (0.1365995 * (cos(fourPi * ratio))) - (0.0106411 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_HARRIS: // blackman harris
			weighingFactor = 0.35875 - (0.48829 * (cos(twoPi * ratio))) + (0.14128 * (cos(fourPi * ratio))) - (0.01168 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_FLT_TOP: // flat top
			weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) + (0.1980399 * cos(fourPi * ratio));
			break;
		case FFT_WIN_TYP_WELCH: // welch
			weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) / (samplesMinusOne / 2.0));
			break;
		}
		if (dir == FFT_FORWARD) {
			_vReal[i] *= weighingFactor;
			_vReal[_samples - (i + 1)] *= weighingFactor;
		}
		else {
			_vReal[i] /= weighingFactor;
			_vReal[_samples - (i + 1)] /= weighingFactor;
		}
	}
	return 0;
}


/**
 * @brief MCU Library Revision
 *
 * @return
 */
uint8_t mcufft_lib_revision(void)
{
	return(FFT_LIB_REV);
}


/**
 * @brief Exponent
 *
 * @param value
 * @return
 */
uint8_t Exponent(uint16_t value)
{
	uint8_t result = 0;
	while (((value >> result) & 1) != 1) result++;
	return(result);
}


/**
 * @brief Swap
 *
 * @param x
 * @param y
 * @return
 */
void Swap(float* x, float* y)
{
	float temp = *x;
	*x = *y;
	*y = temp;
}

