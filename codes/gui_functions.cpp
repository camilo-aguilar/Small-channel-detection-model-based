
#include "gui_functions.h"
#include "QualityCandy.h"

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

	
	CvScalar white = CV_RGB(255,255,255);
	
	
	CvScalar orange = CV_RGB(255,160,0);
	CvScalar cyan = CV_RGB(0,255,255);
	int font_thickness = 0;
	

	cvInitFont(&font1,CV_FONT_HERSHEY_DUPLEX,hscale,vscale,shear,font_thickness,line_type);

	enda.x = mp->enda.x*SCALE;
	enda.y = mp->enda.y*SCALE;
	endb.x = mp->endb.x*SCALE;
	endb.y = mp->endb.y*SCALE;
	cvLine( image, enda, endb, color,line_thickness,8,0);


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
	cvWaitKey(DELAY_FOR_DISPLAY);
	cvReleaseImage(&img);
}

void display_image_double(double **y, int rows, int cols, Candy *C)
{
	IplImage *img = NULL;
	int channel = 3;

	img = cvCreateImage(cvSize(cols,rows),IPL_DEPTH_8U, channel);

	LoadImageFromMemoryDouble(img, y);

	DrawAllLines(C, img, DISPAY_TEXT, 1);


	cvNamedWindow("ChannelMpp");
	cvShowImage("ChannelMpp", img);
	cvWaitKey(DELAY_FOR_DISPLAY);
	cvReleaseImage(&img);
}

void display_only_one_double(double **y, int rows, int cols, lineObj *mp, int color)
{
	IplImage *img = NULL;
	int channel = 3;
	CvScalar green = CV_RGB(0,0,255);
	CvScalar red = CV_RGB(255,0,0);

	img = cvCreateImage(cvSize(cols,rows),IPL_DEPTH_8U, channel);

	LoadImageFromMemoryDouble(img, y);

	if(color == 1)
		DrawCandyLine(img, mp, 0, green,  TEXT_SINGLE_E, 1);
	else
		DrawCandyLine(img, mp, 0, red,  TEXT_SINGLE_E, 1);

	cvNamedWindow("FreeSeg");
	cvShowImage("FreeSeg", img);
	cvWaitKey(DELAY_FOR_DISPLAY);
	cvReleaseImage(&img);
}

void draw_all_nds(NeckDent *mp, int np_num, IplImage *image, int energy_type, double alpha, double lambda_int, int line_thickness)
{
	int k;
	for(k = 0; k<np_num; k++)
		draw_both(&(mp[k]), image, energy_type, alpha, lambda_int, line_thickness);
}


void draw_both(NeckDent *mp, IplImage *image, int text_type, double alpha, double lambda_int, int line_thickness)
{
	CvPoint data_pnts[3][4];
    CvPoint* ppt[3] = { data_pnts[0],data_pnts[1],data_pnts[2] };
	int npts[3];
	double x1, x2, y1, y2, x3, y3, x4, y4;
	// for text
	char text[200];
	double hscale = 0.4;
	double vscale = 0.4;
	double shear = 0.2;
	int line_type = 8;
	CvFont font1;
	CvPoint pt1;

	CvScalar yellow = CV_RGB(255,255,0);
	CvScalar green = CV_RGB(0,255,0);
	CvScalar yellowgreen = CV_RGB(125,255,0);
	CvScalar orange = CV_RGB(255,160,0);
	CvScalar cyan = CV_RGB(0,255,255);
	int font_thickness = 0;
	int connectivity = 8;

	cvInitFont(&font1,CV_FONT_HERSHEY_DUPLEX,hscale,vscale,shear,font_thickness,line_type);

	if(mp->type == NDTYPE_CANDY){
		data_pnts[0][0].x = mp->enda.x*SCALE;
		data_pnts[0][0].y = mp->enda.y*SCALE;
		data_pnts[0][1].x = mp->endb.x*SCALE;
		data_pnts[0][1].y = mp->endb.y*SCALE;
		cvLine( image, data_pnts[0][0], data_pnts[0][1], yellow,line_thickness,8,0);
	}
	else{
		data_pnts[0][0].x = (int)(mp->r[ 0].x*SCALE);	data_pnts[0][0].y = (int)(mp->r[ 0].y*SCALE);
		data_pnts[0][1].x = (int)(mp->r[ 1].x*SCALE);	data_pnts[0][1].y = (int)(mp->r[ 1].y*SCALE);
		data_pnts[0][2].x = (int)(mp->r[ 3].x*SCALE);	data_pnts[0][2].y = (int)(mp->r[ 3].y*SCALE);
		data_pnts[0][3].x = (int)(mp->r[ 2].x*SCALE);	data_pnts[0][3].y = (int)(mp->r[ 2].y*SCALE);
		npts[0] = 4;
		if(mp->type == NDTYPE_DENTING){
			cvPolyLine( image, ppt, npts, 1, 1, yellow,line_thickness,8,0);
			x3 = mp->r[ 0].x + (mp->r[ 1].x - mp->r[ 0].x)/4.;
			y3 = mp->r[ 0].y + (mp->r[ 1].y - mp->r[ 0].y)/4.;
			x4 = mp->r[ 2].x + (mp->r[ 3].x - mp->r[ 2].x)/4.;
			y4 = mp->r[ 2].y + (mp->r[ 3].y - mp->r[ 2].y)/4.;
			x1 = x3 + (x4 - x3)/3.;
			x2 = x3 + (x4 - x3)*2./3.;
			y1 = y3 + (y4 - y3)/3.;
			y2 = y3 + (y4 - y3)*2./3.;
			data_pnts[0][0].x = (int)(x1*SCALE);	data_pnts[0][0].y = (int)(y1*SCALE);
			data_pnts[0][3].x = (int)(x2*SCALE);	data_pnts[0][3].y = (int)(y2*SCALE);
			x1 = mp->r[ 1].x + (mp->r[ 3].x - mp->r[ 1].x)/3.;
			x2 = mp->r[ 1].x + (mp->r[ 3].x - mp->r[ 1].x)*2./3.;
			y1 = mp->r[ 1].y + (mp->r[ 3].y - mp->r[ 1].y)/3.;
			y2 = mp->r[ 1].y + (mp->r[ 3].y - mp->r[ 1].y)*2./3.;
			data_pnts[0][1].x = (int)(x1*SCALE);	data_pnts[0][1].y = (int)(y1*SCALE);
			data_pnts[0][2].x = (int)(x2*SCALE);	data_pnts[0][2].y = (int)(y2*SCALE);
			npts[0] = 4;																	 	
			cvPolyLine( image, ppt, npts, 1, 1, yellow,line_thickness,8,0);
		}
		else if(1)// if(mp->type == NDTYPE_NECKING)
		{
			if(mp->type == NDTYPE_NECKING)
				cvPolyLine( image, ppt, npts, 1, 1, green,line_thickness,8,0);
			else if(mp->type == NDTYPE_BOTH)
				cvPolyLine( image, ppt, npts, 1, 1, yellowgreen,line_thickness,8,0);
			else // NDTYPE_CANDY
				cvPolyLine( image, ppt, npts, 1, 1, cyan,line_thickness,8,0);
			x1 = mp->r[ 0].x + (mp->r[ 2].x - mp->r[ 0].x)/3.;
			x2 = mp->r[ 0].x + (mp->r[ 2].x - mp->r[ 0].x)*2./3.;
			y1 = mp->r[ 0].y + (mp->r[ 2].y - mp->r[ 0].y)/3.;
			y2 = mp->r[ 0].y + (mp->r[ 2].y - mp->r[ 0].y)*2./3.;
			data_pnts[0][0].x = (int)(x1*SCALE);	data_pnts[0][0].y = (int)(y1*SCALE);
			data_pnts[0][3].x = (int)(x2*SCALE);	data_pnts[0][3].y = (int)(y2*SCALE);
			x1 = mp->r[ 1].x + (mp->r[ 3].x - mp->r[ 1].x)/3.;
			x2 = mp->r[ 1].x + (mp->r[ 3].x - mp->r[ 1].x)*2./3.;
			y1 = mp->r[ 1].y + (mp->r[ 3].y - mp->r[ 1].y)/3.;
			y2 = mp->r[ 1].y + (mp->r[ 3].y - mp->r[ 1].y)*2./3.;
			data_pnts[0][1].x = (int)(x1*SCALE);	data_pnts[0][1].y = (int)(y1*SCALE);
			data_pnts[0][2].x = (int)(x2*SCALE);	data_pnts[0][2].y = (int)(y2*SCALE);
			npts[0] = 4;																	 	
			if(mp->type == NDTYPE_NECKING)
				cvPolyLine( image, ppt, npts, 1, 1, green,line_thickness,8,0);
			else if(mp->type == NDTYPE_NECKING)
				cvPolyLine( image, ppt, npts, 1, 1, yellowgreen,line_thickness,8,0);
			else // NDTYPE_CANDY
				cvPolyLine( image, ppt, npts, 1, 1, cyan,line_thickness,8,0);
		}
		else{
			x1 = mp->r[ 0].x + (mp->r[ 2].x - mp->r[ 0].x)/2.;
			y1 = mp->r[ 0].y + (mp->r[ 2].y - mp->r[ 0].y)/2.;
			data_pnts[0][0].x = (int)(x1*SCALE);	data_pnts[0][0].y = (int)(y1*SCALE);
			x1 = mp->r[ 1].x + (mp->r[ 3].x - mp->r[ 1].x)/2.;
			y1 = mp->r[ 1].y + (mp->r[ 3].y - mp->r[ 1].y)/2.;
			data_pnts[0][1].x = (int)(x1*SCALE);	data_pnts[0][1].y = (int)(y1*SCALE);
			cvLine( image, data_pnts[0][0], data_pnts[0][1], orange,line_thickness,8,0);
		}
	}
	pt1 = cvPoint((int)((mp->center.x-7)*SCALE), (int)(mp->center.y*SCALE));//-14
	if(text_type != TEXT_NONE){
		switch(text_type)
		{
			case TEXT_SINGLE_E:
				sprintf(text, "%1.2f",mp->single_E);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_MULTIPLE_E:
				sprintf(text, "%1.2f",mp->multiple_E);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_TOTAL_E:
				sprintf(text, "%1.2f",alpha*mp->single_E + lambda_int*mp->multiple_E);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_BOTH_E:
				sprintf(text, "(%1.2f, %1.2f)",mp->single_E, mp->multiple_E);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_E0_E1:
				if(mp->type == NDTYPE_BOTH)
					sprintf(text, "%1.2f,%1.2f",mp->single_E, mp->e[3]);
				else
					sprintf(text, "%1.2f,%1.2f",mp->single_E, mp->e[4]);
				cvPutText(image,text,pt1,&font1,cyan);
				sprintf(text, "(%1.1f,%1.1f)",mp->e[0], mp->e[1]);
				pt1.x -= 4; pt1.y += 10;
				cvPutText(image,text,pt1,&font1,orange);
				break;
			case TEXT_SINGLE_TEST:
				sprintf(text, "%1.2f",mp->single_E);
				cvPutText(image,text,pt1,&font1,cyan);

				sprintf(text, "(%d,%1.2f,%1.1f)",(int)mp->e[0], mp->e[1], mp->e[2]);
				pt1.x -= 75; pt1.y += 10;
				cvPutText(image,text,pt1,&font1,orange);

				sprintf(text, "(%d,%1.2f,%1.1f)",(int)mp->e[3], mp->e[4], mp->e[5]);
				pt1.x -= 0; pt1.y += 12;
				cvPutText(image,text,pt1,&font1,orange);

				sprintf(text, "(%1.2f,%1.1f)", mp->e[24], mp->e[25]);
				pt1.x -= 0; pt1.y += 12;
				cvPutText(image,text,pt1,&font1,orange);
				break;
			case TEXT_SINGLE_E_DIST:
				sprintf(text, "(%1.2f, %1.2f)",mp->single_E, mp->e[25]); // single energy and error
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_NUM_SINGLE_E:
				sprintf(text, "(%d, %1.2f)",mp->num, mp->single_E);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
			case TEXT_NUMBER:
			default:
				sprintf(text, "%d",mp->num);
				cvPutText(image,text,pt1,&font1,cyan);
				break;
		}
	}
}

void save_neck_dent_char(unsigned char **yimg, NeckDent *mp, MPP_Parameters mpp, char *outfilePrefix, int rows, int cols, int np_num)
{
	char filename[100];
	sprintf(filename, "%s_w_neck_dent_channels.png",outfilePrefix);
	
	IplImage *image = 0;
	int channel = 3;

	image = cvCreateImage(cvSize(cols, rows), IPL_DEPTH_8U, channel);

	LoadImageFromMemoryChar(image, yimg);
	draw_all_nds(mp, np_num, image, SAVE_TEXT, mpp.alpha, mpp.lambda_int, 2);//TEXT_E0_E1, TEXT_NONE

	if(!cvSaveImage(filename, image, 0)){
		printf("Could not save: %s\n",filename);
	}
	cvReleaseImage(&image);
}


void draw_neck_dent(unsigned char** yimg, NeckDent *mp, MPP_Parameters mpp, int rows, int cols, int np_num)
{
	IplImage *image = 0;
	int channel = 3;
	image = cvCreateImage(cvSize(cols, rows), IPL_DEPTH_8U, channel);

	LoadImageFromMemoryChar(image, yimg);
	draw_all_nds(mp, np_num, image, 0, mpp.alpha, mpp.lambda_int, 1);
	cvShowImage("NECK_DENT", image);
	cvWaitKey(DELAY_FOR_DISPLAY);
	cvReleaseImage(&image);
}	