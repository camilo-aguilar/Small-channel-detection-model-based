
#include "gui_functions.h"
#include "qualityCandy.h"

/* In case OPENCV is installed. Heps for Debugging*/
	// Load Image to the Image memory
void LoadImageFromMemoryChar(IplImage *image, unsigned char **y)
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
			for(k = 0; k < channel; k++){
				di = (double)i/SCALE;
				dj = (double)j/SCALE;
				imgdata[i*step + j*channel + k] = (uchar)(real_coord(y,di,dj));
			}
		}
	}
}



void LoadImageFromMemoryDouble(IplImage *image, double **y)
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
	for(i = 0; i < height-1; i++){
		for(j = 0; j < width-1; j++){
			for(k = 0; k < channel; k++){
				di = (double)i;
				dj = (double)j;
				imgdata[i*step + j*channel + k] = (uchar)(real_coord(y,di,dj));
			}
		}
	}
}

//KDW WriteText
void DrawCandyLine(IplImage *image, lineObj *mp, int i, CvScalar color, int text_type, int line_thickness)
{
	CvPoint enda, endb, pt1;
	char text[200];
	double hscale = 0.4;
	double vscale = 0.4;
	double shear = 0.2;
	int line_type = 8;
	CvFont font1;
//	CvScalar red = CV_RGB(255,0,0);
	CvScalar yellow = CV_RGB(255,255,0);
	CvScalar white = CV_RGB(255,255,255);
	CvScalar green = CV_RGB(0,255,0);
	CvScalar yellowgreen = CV_RGB(125,255,0);
	CvScalar orange = CV_RGB(255,160,0);
	CvScalar cyan = CV_RGB(0,255,255);
	int font_thickness = 0;
	int connectivity = 8;

	cvInitFont(&font1,CV_FONT_HERSHEY_DUPLEX,hscale,vscale,shear,font_thickness,line_type);

	enda.x = mp->enda.x*SCALE;
	enda.y = mp->enda.y*SCALE;
	endb.x = mp->endb.x*SCALE;
	endb.y = mp->endb.y*SCALE;
	cvLine( image, enda, endb, color,line_thickness,8,0);

	/*
	enda.x = (mp->enda.x - mp->width/2*sin(mp->theta))*SCALE;
	enda.y = (mp->enda.y - mp->width/2*cos(mp->theta))*SCALE;
	endb.x = (mp->endb.x - mp->width/2*sin(mp->theta))*SCALE;
	endb.y = (mp->endb.y - mp->width/2*cos(mp->theta))*SCALE;
//	cvLine( image, enda, endb, orange,line_thickness,8,0);
	enda.x = (mp->enda.x + mp->width/2*sin(mp->theta))*SCALE;
	enda.y = (mp->enda.y + mp->width/2*cos(mp->theta))*SCALE;
	endb.x = (mp->endb.x + mp->width/2*sin(mp->theta))*SCALE;
	endb.y = (mp->endb.y + mp->width/2*cos(mp->theta))*SCALE;
//	cvLine( image, enda, endb, orange,line_thickness,8,0);
	enda.x = mp->enda.x*SCALE;
	enda.y = mp->enda.y*SCALE;
	endb.x = mp->endb.x*SCALE;
	endb.y = mp->endb.y*SCALE;
//	cvCircle(image, enda, NEIGHBOORHOOD*SCALE, cyan, line_thickness,8,0);
//	cvCircle(image, endb, NEIGHBOORHOOD*SCALE, cyan, line_thickness,8,0);
	*/
	pt1.x = (int)((float)enda.x + (float)endb.x)/2.0;
	pt1.y = (int)((float)enda.y + (float)endb.y)/2.0;

	double dataterm = mp->dataterm;
	
	if(text_type != TEXT_NONE){
		switch(text_type)
		{
			case TEXT_MULTIPLE_E:
				sprintf(text, "%1.2f",mp->engergy_for_transition);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_SINGLE_TEST:
				sprintf(text, "%1.2f",dataterm);
				cvPutText(image,text,pt1,&font1,cyan);

				sprintf(text, "(%d, %d)",(int)mp->x, mp->y);
				pt1.x -= 75; pt1.y += 10;
				cvPutText(image,text,pt1,&font1,orange);

				sprintf(text, "(%d, %d,%1.2f)",(int)mp->len, mp->width, mp->theta);
				pt1.x -= 0; pt1.y += 12;
				cvPutText(image,text,pt1,&font1,orange);
				break;
			case TEXT_SINGLE_E:
				sprintf(text, "%1.2f",dataterm);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_NUMBER:
				sprintf(text, "%d",i);
				cvPutText(image,text,pt1,&font1,white);
				break;
			case TEXT_BOTH_E:
			default:
				sprintf(text, "%.2f, %.2f",dataterm, mp->engergy_for_transition);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
		}
	}
}

//KDW DrawAllLines
void DrawAllLines(Candy *C, IplImage *image, int text_type, int line_thickness)
{
	int n_f = C->n_f;
	int n_s = C->n_s;
	int n_d = C->n_d;
	CvScalar red = CV_RGB(255,0,0);
	CvScalar yellow = CV_RGB(255,255,0);
//	CvScalar white = CV_RGB(255,255,255);
	CvScalar green = CV_RGB(0,255,0);
	CvScalar yellowgreen = CV_RGB(125,255,0);
	CvScalar orange = CV_RGB(255,160,0);
	CvScalar cyan = CV_RGB(0,255,255);
	CvPoint pt1;

	Node *p = C->link_f;

	for (int i = 0;i<n_f;i++)
	{
		p= p->next;
		DrawCandyLine(image, p->index, i, red, text_type, line_thickness);
	}

	p = C->link_s;
	for (int i = 0;i<n_s;i++)
	{
		p= p->next;
		DrawCandyLine(image, p->index, i, yellow, text_type, line_thickness);

	}

	p = C->link_d;
	for (int i = 0;i<n_d;i++)
	{
		p= p->next;
		DrawCandyLine(image, p->index, i, green, text_type, line_thickness);

	}


}



void display_image_char(unsigned char **y, int rows, int cols)
{
	IplImage *img = NULL;
	
	LoadImageFromMemoryChar(img, y);
	cvNamedWindow("ChannelMpp");
	cvShowImage("ChannelMpp", img);
	cvWaitKey(0);
	cvReleaseImage(&img);

}

void display_image_double(double **y, int rows, int cols, Candy *C)
{
	IplImage *img = NULL;
	int channel = 3;

	img = cvCreateImage(cvSize(cols,rows),IPL_DEPTH_8U, channel);

	LoadImageFromMemoryDouble(img, y);

	DrawAllLines(C, img, TEXT_NONE, 1);


	cvNamedWindow("ChannelMpp");
	cvShowImage("ChannelMpp", img);
	cvWaitKey(0);
	cvReleaseImage(&img);

}

void display_only_one_double(double **y, int rows, int cols, lineObj *mp)
{
	IplImage *img = NULL;
	int channel = 3;
	CvScalar green = CV_RGB(0,0,255);

	img = cvCreateImage(cvSize(cols,rows),IPL_DEPTH_8U, channel);

	LoadImageFromMemoryDouble(img, y);

	DrawCandyLine(img, mp, 0, green,  TEXT_SINGLE_E, 1);


	cvNamedWindow("FreeSeg");
	cvShowImage("FreeSeg", img);
	cvWaitKey(0);
	cvReleaseImage(&img);

}
