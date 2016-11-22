#ifndef _SYNTHETIC_
#define _SYNTHETIC_

#include "ndmulti_e.h"
#include "tiff.h"
#include "util.h"

#define SYN_CH_NUM	200

typedef struct ellipse_param
{
	DPoint			center;
	double			a;
	double			b;
	double			n;	   // super ellipse curvature
	int				color; // 1 = white, 0 = black
} Ellipse_Param;

typedef struct neckdent_syn
{
	int				num;
	ND_Type			type;		// denting area, necking area
	DPoint			center;
	double			syn_width;	// synthetic width
	double			syn_length;	// synthetic length
	double			width;		// channel width
	double			length;		// channel length
	double			theta;		// direction
	double			curvature;	// super ellipse curvature
	double			amp;
	double			offset;
	DPoint			r[4];
	double			radius;
	DPoint			lines[4];	// x = m, y = b : y = mx+b
	Ellipse_Param	ellipse[4];
} NeckDent_Syn;

void gen_syn_image(unsigned char **img, double *mean, double stdev, NeckDent_Syn *syn_ch, int *syn_num,
				   MPP_Parameters mpp, int cols, int rows, double **blur, int blur_size, char *infileName, char *gtfileName);
void gen_syn_image2(unsigned char **img, double *mean, double stdev, NeckDent_Syn *syn_ch, int *syn_num,
				   MPP_Parameters mpp, int cols, int rows, double **blur, int blur_size, char *infileName, char *gtfileName);
#endif
