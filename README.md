# mcu-fft
Fast Fourier Transform implementation for Arduino and MCUs

## Features

- Remotely invoke functions on device.
- Ability to process function parameters.
- Statically allocated memory.
- Backspace to remove unintentional keypresses.

## Introduction

This package contains files to implement a simple command-line interface.
The package includes mcufft.h, and mcufft.c.

## Integration details

- Integrate mcufft.h and mcufft.c files into your project.
- Include the mcufft.h header file in your code like below.

```c
#include "mcufft.h"
```

## File information

- mcufft.h : This header file contains the definitions of the cli user API.
- mcufft.c : This source file contains the implementation of the CLI.

## References
- [arduinoFFT](https://github.com/kosme/arduinoFFT)
- [microFFT](https://github.com/marlonValerio/microFFT)
