#ifndef _UTIL_
#define _UTIL_

#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
//#include "cv.h"
//#include "cxcore.h"
//#include "highgui.h"
#include "allocate.h"
#include "random.h"
#include "tiff.h"

#define INF			1000000000
#define SCALE		3

typedef struct point
{
	int		x;
	int		y;
} Point;

typedef struct dpoint
{
	double		x;
	double		y;
} DPoint;

//#define delta(a)	(double)((a)-(int)(a))
//#define real_coord(f,a,b)	(double)(f[(int)(a)+1][(int)(b)+1])*(delta(a))*(delta(b)) -(double)(f[(int)(a)][(int)(b)+1])*(delta(a)-1)*(delta(b)) -(double)(f[(int)(a)+1][(int)(b)])*(delta(b)-1)*(delta(a)) +(double)(f[(int)(a)][(int)(b)])*(delta(b)-1)*(delta(a)-1)
#define delta(a)	((a)-floor(a))
#define real_coord(f,a,b)	((double)(f[(int)(a)+1][(int)(b)+1])*(delta(a))*(delta(b)) -(double)(f[(int)(a)][(int)(b)+1])*(delta(a)-1)*(delta(b)) -(double)(f[(int)(a)+1][(int)(b)])*(delta(b)-1)*(delta(a)) +(double)(f[(int)(a)][(int)(b)])*(delta(b)-1)*(delta(a)-1))
#define rotatex(x,y,t) (double)(x)*cos(t) - (double)(y)*sin(t)
#define rotatey(x,y,t) (double)(x)*sin(t) + (double)(y)*cos(t)
#define distance2(x1,y1,x2,y2) ((double)(x1)-(double)(x2))*((double)(x1)-(double)(x2))+((double)(y1)-(double)(y2))*((double)(y1)-(double)(y2))
#define distance1(x1,y1,x2,y2) sqrt(distance2(x1,y1,x2,y2))

// generate (min,max]
#define ubnd_random(min,max) (random2()*((double)(max)-(double)(min))+(double)(min))
// generate [min,max)
#define lbnd_random(min,max) ((1.-random2())*((double)(max)-(double)(min))+(double)(min))
// generate [min,max]
#define bbnd_rand(min,max) ((double)(min) + ((double)rand()*((double)(max)-(double)(min))/RAND_MAX))
// generate integer [min,max]
#define int_random(min,max)  (int)((1.-random2())*((double)((max)+1.)-(double)(min))+(double)(min))

#ifndef max
#define max(a,b)	(((a) > (b))? (a) : (b))
#endif

#ifndef min
#define min(a,b)	(((a) < (b))? (a) : (b))
#endif

//void LoadImageFromMemory(IplImage *image, unsigned char **y);
//void LoadImageFromMemory(IplImage *image, int **y);
//void LoadImageFromMemory(IplImage *image, double **y);
//void clearImage(IplImage *image);
double dist(DPoint a, DPoint b);
void double_to_uchar(double **lm, unsigned char **uclm, int cols, int rows);


#endif
