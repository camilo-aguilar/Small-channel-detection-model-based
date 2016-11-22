#include <math.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <stdlib.h>
#include "allocate.h"
#include "tiff.h"
#include "random.h"
//#define PI  3.14159265358979323846
#define MAX_CLASSES 15
#undef MPM_IMG_OUTPUT
#define DEF_GAUSSIAN 		0
//#define DIFF_SUM

#ifdef DEF_COVARIANCE
void mpm(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double *beta, double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double delta, double **diff, double *s_p, double ****kmat,
		double *vmat);
void mpm_betaimg(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double **beta[], double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double delta, double **diff, double *s_p, double ****kmat,
		double *vmat);
#else
void mpm(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double *beta, double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double **diff, double *s_p);
void mpm_betaimg(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double **beta[], double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double **diff, double *s_p);
#endif
double gaussianRandom(double average, double stdev);
