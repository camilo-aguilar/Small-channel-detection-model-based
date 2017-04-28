#ifndef _NDMULTI_
#define _NDMULTI_

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
#include "util.h"
#include "em.h"

#include "Parameters.h"

#define USE_BIRTHMAP
#define MAX_MKPNT_NUM		1000//220		// Max marked point number
//#define SIMPLE_CH
#ifdef SIMPLE_CH
	#define RJMCMC
#endif
#define RJMCMC
#define RJMCMC_ITER_DIV		1000
//#define TEST_SINGLE_E
#define	IMG_SYNTHETIC	0// synthetic
#define	IMG_BSE			1// BSEImage
#define IMG_SLICE001	2// slice_crop017
#define IMG_SLICE017	3// slice_crop017
#define IMG_SLICE170	4// slice_crop170
#define IMG_IMG0		5// img0_crop
#define IMG_IMG1		6// img0_crop
#define IMG_IMG1_T		7// img0_crop
#define INPUT_IMAGE		IMG_IMG1_T
#define DEBUG_SYN_CH_EN
#define FIXED_ANGLE	// fixed angle

#define HARDCORE_REPULSION	4 // 2
#define MAX_MPP_ITER_NUM	100000
#define OBJECT_IN_RATIO		0.9

#define STEP_DX					0.5
#define STEP_DY					0.5

#define JMP_BIRTH			0
#define JMP_DEATH			1
#define JMP_TRANSLATION		2
#define JMP_DILATION		3
#define JMP_ROTATION		4
#define JMP_SWITCHING		5
#define NUM_JUMP_TYPE		6

#define DELTA_TRANSLATION	1.5//2
#define DELTA_DILATION		1.5//2
#define DELTA_ROTATION		0.02//0.5

//#define DENT_NECK_S_E_F		0.8

#define TEXT_NONE			-1
#define TEXT_SINGLE_E		0
#define TEXT_MULTIPLE_E		1
#define TEXT_TOTAL_E		2
#define TEXT_BOTH_E			3
#define TEXT_E0_E1			4
#define TEXT_SINGLE_E_DIST	5
#define TEXT_SINGLE_TEST	6
#define TEXT_NUM_SINGLE_E	7
#define TEXT_NUMBER			8
#define RUNNGIN_TEXT		TEXT_SINGLE_E_DIST//TEXT_NONE
#define SAVE_TEXT			TEXT_NONE

// for avg_Gaussian_neck and avg_Gaussian_dent
#define WEIGHT_AMP	4
#define WEIGHT_OFFSET 1
#define NECK_DISCON_TH	-0.08
#define DENT_DISCON_TH	2 //-0.07

#define f_exp(a, b, y2) (a)*(1-exp(y2))+(b)

typedef enum
{
	STATE_NON_EXIST,
	STATE_NEW_BORN,
	STATE_EXIST,
} MP_State;

typedef enum
{
	NDTYPE_BOTH,
	NDTYPE_NECKING,
	NDTYPE_DENTING,
	NDTYPE_CANDY,
} ND_Type;

#define ERR_NUM 30
typedef struct neckdent
{
	int				num;
	MP_State		state;
	ND_Type			type;		// denting area, necking area
	DPoint			center;
	Point			enda;
	Point			endb;
	double			width;
	double			length;
	double			theta;		// direction
//	double			length; 
//	double			width;
	DPoint			r[12];
	double			single_E;	// single energy
	double			multiple_E;	// multiple energy
	int				sort_idx;
	double			e[ERR_NUM];
} NeckDent;



int nd_mpp_rjmcmc(unsigned char **yimg, double **lm, double *mean, double *vari, 
		   double variance, NeckDent *mp, 
			MPP_Parameters mpp, double *total_e, int *mp_num, 
		   int cols, int rows);

/*
int nd_mpp_multiple_birth_n_death(unsigned char **yimg, double **lm, double *mean,
			double *vari, double variance, NeckDent *mp,
			MPP_Parameters mpp, double *total_e, int *mp_num, 
			int cols, int rows, IplImage *image, const char* win_name);



void test_single_object(unsigned char **yimg, double *mean, double *vari, double variance,
						MPP_Parameters mpp, 
						int cols, int rows, IplImage *image, const char* win_name);


void draw_all_nds(NeckDent *mp, int np_num, IplImage *image, int energy_type, double alpha,
				  double lambda_int, int line_thickness);

*/
void calculate_betaimg(double **beta_img[], double *beta, NeckDent *mp, 
					   int np_num, MPP_Parameters mpp, int cols, int rows);


void draw_all_nds_2D(NeckDent *mp, MPP_Parameters mpp, int np_num, unsigned char **image,
					 unsigned char **image2, int cols, int rows);
void nd_GenCornerPoint(NeckDent *mp);


int nd_mpp_multiple_birth_n_death(unsigned char **yimg, double **lm, double *mean,
			double *vari, double variance, NeckDent *mp,
			MPP_Parameters mpp, double *total_e, int *mp_num, 
			int cols, int rows);
			
#endif
