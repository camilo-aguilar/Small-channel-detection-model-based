#ifndef _GUI_FUNCTIONS_
#define _GUI_FUNCTIONS_

#include <opencv2/opencv.hpp>
#include "util.h"
#include "QualityCandy.h"


void display_image_double(double **y,int rows, int cols, Candy *C);
void display_only_one_double(double **y, int rows, int cols, lineObj *mp, int color);


void save_neck_dent_char(unsigned char **yimg, NeckDent *mp, MPP_Parameters mpp, char *outfilePrefix, int rows, int cols, int np_num);

void draw_all_nds(NeckDent *mp, int np_num, IplImage *image, int energy_type, double alpha,
				  double lambda_int, int line_thickness);
void draw_both(NeckDent *mp, IplImage *image, int text_type, double alpha, double lambda_int, int line_thickness);

void draw_neck_dent(unsigned char** yimg, NeckDent *mp, MPP_Parameters mpp, int rows, int cols, int np_num);

#endif