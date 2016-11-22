#include "util.h"
#define CAMILO_OPENCV 0

#if CAMILO_OPENCV
// Load Image to the Image memory
void LoadImageFromMemory(IplImage *image, unsigned char **y)
{
	int i,j,k;
	int step, channel, height, width;
	uchar *imgdata;
	double di,dj;

	step = image->widthStep;
	imgdata = (uchar *)image->imageData;
	channel = image->nChannels;
	height = image->height;
	width = image->width;
	for(i = 0; i < height-SCALE; i++){
		for(j = 0; j < width-SCALE; j++){
			di = (double)i/SCALE;
			dj = (double)j/SCALE;
			for(k = 0; k < channel; k++){
				imgdata[i*step + j*channel + k] = (uchar)(real_coord(y,di,dj));
			}
		}
	}
	for(i = height-SCALE; i < height; i++){
		for(j = 0; j < width-SCALE; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[(height-SCALE-1)*step + j*channel + k];
		}
	}
	for(i = 0; i < height-SCALE; i++){
		for(j = width-SCALE; j < width; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[i*step + (width-SCALE-1)*channel + k];
		}
	}
	for(i = height-SCALE; i < height; i++){
		for(j = width-SCALE; j < width; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[(height-SCALE-1)*step + (width-SCALE-1)*channel + k];
		}
	}
}


void LoadImageFromMemory(IplImage *image, int **y)
{
	int i,j,k;
	int step, channel, height, width;
	uchar *imgdata;
	double di,dj;

	step = image->widthStep;
	imgdata = (uchar *)image->imageData;
	channel = image->nChannels;
	height = image->height;
	width = image->width;
	for(i = 0; i < height-SCALE; i++){
		for(j = 0; j < width-SCALE; j++){
			di = (double)i/SCALE;
			dj = (double)j/SCALE;
			for(k = 0; k < channel; k++){
				imgdata[i*step + j*channel + k] = (uchar)(real_coord(y,di,dj));
			}
		}
	}
	for(i = height-SCALE; i < height; i++){
		for(j = 0; j < width-SCALE; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[(height-SCALE-1)*step + j*channel + k];
		}
	}
	for(i = 0; i < height-SCALE; i++){
		for(j = width-SCALE; j < width; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[i*step + (width-SCALE-1)*channel + k];
		}
	}
	for(i = height-SCALE; i < height; i++){
		for(j = width-SCALE; j < width; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[(height-SCALE-1)*step + (width-SCALE-1)*channel + k];
		}
	}
}

void LoadImageFromMemory(IplImage *image, double **y)
{
	int i,j,k;
	int step, channel, height, width;
	uchar *imgdata;
	double di,dj;

	step = image->widthStep;
	imgdata = (uchar *)image->imageData;
	channel = image->nChannels;
	height = image->height;
	width = image->width;
	for(i = 0; i < height-SCALE; i++){
		for(j = 0; j < width-SCALE; j++){
			di = (double)i/SCALE;
			dj = (double)j/SCALE;
			for(k = 0; k < channel; k++){
				imgdata[i*step + j*channel + k] = (uchar)(real_coord(y,di,dj));
			}
		}
	}
	for(i = height-SCALE; i < height; i++){
		for(j = 0; j < width-SCALE; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[(height-SCALE-1)*step + j*channel + k];
		}
	}
	for(i = 0; i < height-SCALE; i++){
		for(j = width-SCALE; j < width; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[i*step + (width-SCALE-1)*channel + k];
		}
	}
	for(i = height-SCALE; i < height; i++){
		for(j = width-SCALE; j < width; j++){
			for(k = 0; k < channel; k++)
				imgdata[i*step + j*channel + k] = imgdata[(height-SCALE-1)*step + (width-SCALE-1)*channel + k];
		}
	}
}

void clearImage(IplImage *image)
{
	int i,j,k;
	int step, channel, height, width;
	uchar *imgdata;
	double di,dj;

	step = image->widthStep;
	imgdata = (uchar *)image->imageData;
	channel = image->nChannels;
	height = image->height;
	width = image->width;
	for(i = 0; i < height-SCALE; i++){
		for(j = 0; j < width-SCALE; j++){
			di = (double)i/SCALE;
			dj = (double)j/SCALE;
			for(k = 0; k < channel; k++){
				imgdata[i*step + j*channel + k] = 0;
			}
		}
	}
}
#endif

double dist(DPoint a, DPoint b)
{
	return sqrt((a.x - b.x)*(a.x - b.x)+(a.y - b.y)*(a.y - b.y));
}


void double_to_uchar(double **lm, unsigned char **uclm, int cols, int rows)
{
	double min, max;
	int i,j;

	min = 10000;
	max = 0;
	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			if(lm[i][j]>max) max = lm[i][j];
			if(lm[i][j]<min) min = lm[i][j];
		}
	}
	if(min==max){
		min = 0;
	}

	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			uclm[i][j] = (unsigned char)((lm[i][j] - min)/(max-min)*255.);
		}
	}
}
