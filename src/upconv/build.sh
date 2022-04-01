#!/bin/sh

gcc -std=c99 -static -Wl,--stack,10485760 -Os -msse2 -fopenmp -ftree-vectorize -o upconv.exe upconv.c fft_filter.c fileio.c mconv.c wav2raw.c start_exec.c ./../PLG_AUDIO_IO/PLG_AUDIO_IO.c -lfftw3-3 -liconv
