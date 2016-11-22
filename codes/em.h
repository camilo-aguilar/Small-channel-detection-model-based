#include <math.h>
#include <stdlib.h>
#include "tiff.h"
#include "allocate.h"
#include "random.h"

#pragma warning(disable : 4996)

//#define DEF_COVARIANCE
#undef DEF_COVARIANCE
#define GAUSSIAN_BLUR
//#undef GAUSSIAN_BLUR
//#define EXPONETIAL_BLUR
#undef EXPONETIAL_BLUR

#define MAX_CLASSES 15
#define CLASSES 2
#define MAX_LEVEL 20
#define MAX_MPMITER 100
#define MAX_EMITER 100
//#define EM_IMG_OUTPUT
#undef EM_IMG_OUTPUT
#define REPET 1	// # of repetition

void emmpm(unsigned char **y, unsigned char **xt, double *beta, double *gamma, int emiter, int mpmiter,
		int rows, int cols, int classes, double **blur, int blur_size, int enable_blur);
void emmpm_betaimg(unsigned char **y, unsigned char **xt, double **beta[], double *gamma, int emiter, int mpmiter,
		int rows, int cols, int classes, double **blur, int blur_size, int enable_blur);
