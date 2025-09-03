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
 * @brief MCU init
 *
 * @param vReal
 * @param vImag
 * @param samples
 * @param sampling_frequency
 * @return
 */
void mcufft_init(float* vReal, float* vImag, uint16_t samples, float sampling_frequency)
{
	_vReal = vReal;
	_vImag = vImag;
	_samples = samples;
	_sampling_frequency = sampling_frequency;
	_power = Exponent(samples);
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
 * @brief Compute the magnitude of the complex number
 *
 * @return
 */
void mcufft_complex_to_magnitude(void)
{ // vM is half the size of vReal and vImag
	for (uint16_t i = 0; i < _samples; i++) {
		_vReal[i] = sqrt(sq(_vReal[i]) + sq(_vImag[i]));
	}
}


/**
 * @brief Apply the window function
 *
 * @param windowType
 * @param dir
 * @return
 */
void mcufft_windowing(uint8_t windowType, uint8_t dir)
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
}


/**
 * @brief Major Peak
 *
 * @param mag_out
 * @param freq_out
 * @param mag_fact
 * @return
 */
void mcufft_major_peak(float* mag_out, float* freq_out, float mag_fact)
{
	float maxY = 0;
	uint16_t IndexOfMaxY = 0;

	for (uint16_t i = 1; i < ((_samples >> 1) + 1); i++) {
		if ((_vReal[i - 1] < _vReal[i]) && (_vReal[i] > _vReal[i + 1])) {
			if (_vReal[i] > maxY) {
				maxY = _vReal[i];
				IndexOfMaxY = i;
			}
		}
	}
	float delta = 0.5 * ((_vReal[IndexOfMaxY - 1] - _vReal[IndexOfMaxY + 1]) / (_vReal[IndexOfMaxY - 1] - (2.0 * _vReal[IndexOfMaxY]) + _vReal[IndexOfMaxY + 1]));
	float interpolatedX = ((IndexOfMaxY + delta)  * _sampling_frequency) / (_samples - 1);
	if(IndexOfMaxY == (_samples >> 1))
		interpolatedX = ((IndexOfMaxY + delta)  * _sampling_frequency) / (_samples);

	*mag_out = _vReal[IndexOfMaxY] / mag_fact;
	*freq_out = interpolatedX;
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
