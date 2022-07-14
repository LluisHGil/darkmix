# Basic makefile for project 'darkmix'
#
# Author: Pablo Arnalte-Mur (OAUV)
#

CC = g++
CCC = $(CC) -c

LIBS =
LIBS_FFTW = -lfftw3 -lm

CFLAGS = -Wall -O3

.SUFFIXES: .o .c .h

OBJ =

FILES = kernel_lambda kernel_absolute_residuals kernel_fourier_residuals

.PHONY: all clean cleanall

all: $(FILES)

.c.o:

	$(CCC)  -o $*.o $<


kernel_lambda: kernel_lambda.cpp $(OBJ)
	@echo "Compiling kernel_lambda"
	$(CC) $(CFLAGS) kernel_lambda.cpp -o kernel_lambda $(OBJ) $(LIBS)
	@echo "Done"

kernel_absolute_residuals: kernel_absolute_residuals.cpp $(OBJ)
	@echo "Compiling kernel_absolute_residuals"
	$(CC) $(CFLAGS) kernel_absolute_residuals.cpp -o kernel_absolute_residuals $(OBJ) $(LIBS)
	@echo "Done"

kernel_fourier_residuals: kernel_fourier_residuals.cpp $(OBJ)
	@echo "Compiling kernel_fourier_residuals"
	$(CC) $(CFLAGS) kernel_fourier_residuals.cpp -o kernel_fourier_residuals $(OBJ) $(LIBS) $(LIBS_FFTW)
	@echo "Done"

clean:
	@rm -f *~ "#*" *.o

cleanall:
	@rm -f *~ "#*" *.o kernel_lambda kernel_absolute_residuals kernel_fourier_residuals

	
