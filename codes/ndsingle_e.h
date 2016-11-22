#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "allocate.h"
#include "random.h"
#include "tiff.h"

int Single_Object(unsigned char **y, Curve *omega, double GD_step, double delta_t, 
				  double lambda_g, double lambda_G, double epsilon,
				  double mu0, double mu1, double vari0, double vari1,
				  int cols, int rows, IplImage *image, const char* win_name);
