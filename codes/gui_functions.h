#ifndef _GUI_FUNCTIONS_
#define _GUI_FUNCTIONS_

#include <opencv2/opencv.hpp>
#include "util.h"
#include "qualityCandy.h"


void display_image_double(double **y,int rows, int cols, Candy *C);
void display_only_one_double(double **y, int rows, int cols, lineObj *mp);


#endif