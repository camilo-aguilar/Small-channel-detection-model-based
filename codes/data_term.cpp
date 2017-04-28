#include "QualityCandy.h"
//KDW #include "graph.h"
#include "tiff.h"
#include "allocate.h"
#include "util.h" //KDW


#if INCLUDE_OPENCV 
	#include "gui_functions.h"
#endif


#define GAUSS_FILTER_LENGTH 5



/*
	Filters Image at (i,j) with Rotated Combination of Filter_x, Filter_y as:
	[cos(theta) sin(theta); -sin(theta) cos(theta) ][Filter_x Filter_y]^T
*/
double im_filter(int **img_seg, int rows, int cols, double theta, int i_coord, int j_coord)
{
#if 0
	static double test_Gaus_x[GAUSS_FILTER_LENGTH][GAUSS_FILTER_LENGTH] = 
{ {0.0058 ,   0.0261,   0.0431,   0.0261,   0.0058},
  {0.0131,   0.0585,   0.0965,   0.0585,   0.0131},
  {  0  ,       0 ,       0 ,        0,        0},
  {-0.0131,   -0.0585,   -0.0965,   -0.0585,   -0.0131},
   {-0.0058,   -0.0261,   -0.0431,   -0.0261,   -0.0058}
};

static double test_Gaus_y[GAUSS_FILTER_LENGTH][GAUSS_FILTER_LENGTH] = 
{   {0.0058  ,  0.0131    ,     0  , -0.0131 ,  -0.0058},
    {0.0261  ,  0.0585    ,     0  , -0.0585 ,  -0.0261},
    {0.0431  ,  0.0965    ,     0  , -0.0965 ,  -0.0431},
    {0.0261  ,  0.0585    ,     0  , -0.0585 ,  -0.0261},
    {0.0058  ,  0.0131    ,     0  , -0.0131 ,  -0.0058}
};
#endif
static double test_Gaus_xx[GAUSS_FILTER_LENGTH][GAUSS_FILTER_LENGTH] = 
{
	{0.0087 ,   0.0392   , 0.0646 ,   0.0392  ,  0.0087},
    {     0   ,      0    ,     0   ,      0    ,     0},
    {-0.0215 ,  -0.0965  , -0.1592  , -0.0965  , -0.0215},
    {     0  ,       0  ,       0   ,      0    ,     0},
    {0.0087  ,  0.0392 ,   0.0646   , 0.0392   , 0.0087}
};

static double test_Gaus_yy[GAUSS_FILTER_LENGTH][GAUSS_FILTER_LENGTH] = 
{
{    0.0087   ,      0 ,   -0.0215   ,      0 ,   0.0087},
{    0.0392   ,      0 ,  -0.0965    ,     0  ,  0.0392},
{    0.0646   ,      0 ,  -0.1592    ,     0  ,  0.0646},
{    0.0392   ,      0 ,  -0.0965    ,     0  ,  0.0392},
{    0.0087   ,      0 ,  -0.0215    ,     0  ,  0.0087}
};

static double test_Gays_xy[GAUSS_FILTER_LENGTH][GAUSS_FILTER_LENGTH] =
{
{ 0.0117    ,0.0261       ,  0 ,  -0.0261  , -0.0117},
{    0.0261  ,  0.0585     ,    0 ,  -0.0585 ,  -0.0261},
{         0  ,       0     ,    0  ,       0 ,        0},
{   -0.0261  , -0.0585     ,    0  ,  0.0585 ,   0.0261},
{   -0.0117  , -0.0261    ,     0  ,  0.0261 ,   0.0117}
};

	int i=0,j=0;
	int x=0,y=0;
	int filter_l = (GAUSS_FILTER_LENGTH-1.0)/2.0;
	double sum = 0;
	double result = 0;

	double cos_t = cos(theta);
	double sin_t = sin(theta);

	for (i=-filter_l; i <= filter_l; i++)
	{
		for(j=-filter_l; j<= filter_l; j++)
		{
			x = i_coord-i;
			if(x < 0)
				x = rows -1 - i;
			else if(x > rows-1)
				x = -i -1;
			
			y =j_coord-j;
			if(y < 0)
				y = cols -1 -j;
			else if(y > cols-1)
				y = -j -1;

			//sum += cos_t * test_Gaus_x[i+2][j+2]*(double)img_seg[x][y] + sin_t * test_Gaus_y[i+2][j+2]*(double)img_seg[x][y];		
			sum += (255.00 - (double)img_seg[x][y]) * (cos_t*cos_t*test_Gaus_xx[i+2][j+2] - 4*sin_t*cos_t*test_Gays_xy[i+2][j+2] + sin_t*sin_t*test_Gaus_yy[i+2][j+2]);
		}
	}
	result = sum;
	return result;

}


site drawRect(double **patch, int patch_len, int y, int x, int L, double theta, int upper_w, int lower_w)
{
	site offset;

	int y2 = (int)floor(y-(double)L*sin(theta)+0.5);
	int x2 = (int)floor(x+(double)L*cos(theta)+0.5);

	int x_upper = (int)floor(x-(double)upper_w/2.0*sin(theta)+0.5);
	int y_upper = (int)floor(y-(double)upper_w/2.0*cos(theta)+0.5);

	int x_upper2 = (int)floor(x2-(double)upper_w/2.0*sin(theta)+0.5);
	int y_upper2 = (int)floor(y2-(double)upper_w/2.0*cos(theta)+0.5);

	int x_upper_down = (int)floor(x+(double)upper_w/2.0*sin(theta)+0.5);
	int y_upper_down = (int)floor(y+(double)upper_w/2.0*cos(theta)+0.5);

	int x_upper_down2 = (int)floor(x2+(double)upper_w/2.0*sin(theta)+0.5);
	int y_upper_down2 = (int)floor(y2+(double)upper_w/2.0*cos(theta)+0.5);

	int x_lower = (int)floor(x_upper_down+(double)lower_w*sin(theta)+0.5);
	int y_lower = (int)floor(y_upper_down+(double)lower_w*cos(theta)+0.5);

	int x_lower2 = (int)floor(x_upper_down2+(double)lower_w*sin(theta)+0.5);
	int y_lower2 = (int)floor(y_upper_down2+(double)lower_w*cos(theta)+0.5);

	int x_left = (int)floor(x_upper-(double)lower_w*sin(theta)+0.5);
	int y_left = (int)floor(y_upper-(double)lower_w*cos(theta)+0.5);

	int x_left2 = (int)floor(x_upper2-(double)lower_w*sin(theta)+0.5);
	int y_left2 = (int)floor(y_upper2-(double)lower_w*cos(theta)+0.5); 

	int x_mid = int((x_left + x_lower2)/2.0);
	int y_mid = int((y_left + y_lower2)/2.0);

	offset.x = (int)(x_mid - patch_len/2.0-1);
	offset.y = (int)(y_mid - patch_len/2.0-1);
	int offsetX = offset.x;
	int offsetY = offset.y;

	line(x_upper_down-offsetX,y_upper_down-offsetY,x_upper-offsetX,y_upper-offsetY,128,patch);
	line(x_upper_down-offsetX,y_upper_down-offsetY,x_upper_down2-offsetX,y_upper_down2-offsetY,128,patch);
	line(x_upper-offsetX,y_upper-offsetY,x_upper2-offsetX,y_upper2-offsetY,128,patch);
	line(x_upper_down2-offsetX,y_upper_down2-offsetY,x_upper2-offsetX,y_upper2-offsetY,128,patch);


	subfill(patch,int((x_upper_down2+x_upper)/2.0)-offsetX,int((y_upper_down2+y_upper)/2.0)-offsetY,128);

	line(x_upper_down-offsetX,y_upper_down-offsetY,x_lower-offsetX,y_lower-offsetY,255,patch);
	line(x_upper_down-offsetX,y_upper_down-offsetY,x_upper_down2-offsetX,y_upper_down2-offsetY,255,patch);
	line(x_lower-offsetX,y_lower-offsetY,x_lower2-offsetX,y_lower2-offsetY,255,patch);
	line(x_upper_down2-offsetX,y_upper_down2-offsetY,x_lower2-offsetX,y_lower2-offsetY,255,patch);

	subfill(patch,int((x_upper_down+x_lower2)/2.0)-offsetX,int((y_upper_down+y_lower2)/2.0)-offsetY,255);

	line(x_upper-offsetX,y_upper-offsetY,x_left-offsetX,y_left-offsetY,100,patch);
	line(x_upper-offsetX,y_upper-offsetY,x_upper2-offsetX,y_upper2-offsetY,100,patch);
	line(x_left-offsetX,y_left-offsetY,x_left2-offsetX,y_left2-offsetY,100,patch);
	line(x_left2-offsetX,y_left2-offsetY,x_upper2-offsetX,y_upper2-offsetY,100,patch);

	subfill(patch,int((x_upper+x_left2)/2.0)-offsetX,int((y_upper+y_left2)/2.0)-offsetY,100);

	return offset;
}



void subfill (double **tmpmask, int x, int y,int new_c )
{
	if (tmpmask[y][x] != new_c)
	{	tmpmask[y][x] = new_c;
	subfill (tmpmask, x+1, y,new_c );
	subfill (tmpmask, x-1, y,new_c );
	subfill (tmpmask, x, y+1,new_c );
	subfill (tmpmask, x, y-1,new_c );
	}

}
 


#define STEP_DX					0.5
#define STEP_DY					0.5
#define SYM_TH					0.5//0.30//0.09
#define OBJECT_IN_RATIO		0.9
/***********************************************************************/
//
//	n1n1n1n1n1n1n1n1n1n1n1
//						
//	c2c2c2c2c2c2c2c2c2c2c2
//						
//	n2n2n2n2n2n2n2n2n2n2n2
//
/***********************************************************************/
double dataterm_rec(int **img_seg, lineObj *line, double **patch, int patch_len, int line_w,double ***Matrix, double prior_pen1, double prior_pen2, int height,int width)
{
	int **y;
	int i, num_out = 0, num_s = 0, num;
	int num_ch = 0, num_nch = 0;
	int num_1ch = 0, num_2ch = 0, num_3ch = 0;
	int num_n1ch = 0, num_n2ch = 0;
	double di, dj, w, l;
	double w_5, l_5;
	double dtmp, dtmp2;
	double sum_all = 0, sum_all2 = 0, sum_s = 0;
	double sum_ch = 0, sum_ch2 = 0;
	double sum_1ch = 0, sum_1ch2 = 0;
	double sum_2ch = 0, sum_2ch2 = 0;
	double sum_3ch = 0, sum_3ch2 = 0;
	double sum_nch = 0, sum_nch2 = 0;
	double sum_n1ch = 0, sum_n1ch2 = 0;
	double sum_n2ch = 0, sum_n2ch2 = 0;
	double likely=-10000;
	double t_test = 0, t_test1 = 0;//, t_test2; 
	double sin_t, cos_t, di_sin_t, di_cos_t;
	double cx, cy, dj_cos_t[100], dj_sin_t[100];
	DPoint pt;
	int flag;
	int cols = width, rows = height;
	double e[30];
	double t = -line->theta; //KDW Huixi's definition of theta is different from mine !!!
	double error_th = prior_pen1;

	y  = img_seg;
	sin_t = sin(t);
	cos_t = cos(t);
	l = (double)line->len;
	w = (double)line_w;
	w_5 = w/2.;
	l_5 = l/2.;
	cx = (double)line->x;
	cy = (double)line->y;

	double derivative_energy =0;
	double derivative_threshold = DERIVATIVE_THRESHOLD;
	double derivative_likely;

	num = 0;

	/*	Calculate length x and y components at angle theta */
	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++)
	{
		dj_cos_t[i] = dj*cos_t;
		dj_sin_t[i] = dj*sin_t;
	}
	
	/*	Calculate width x and y components at angle theta */
	di = -w_5;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;

	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++)
	{
		flag = 0;
		
		//Find sum of intensity, sum of itensity^2 and N in  non channel 1
		pt.x = cx + dj_cos_t[i] + di_sin_t;
		pt.y = cy + dj_sin_t[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
		{
			dtmp = real_coord(y,pt.y,pt.x);
			sum_n1ch += dtmp;
			sum_n1ch2 += dtmp*dtmp;
			num_n1ch++;
			num++;
			flag = 1;
		}
		else
			num_out++;
		
		//Find intensity and N in  non channel 2
		pt.x = cx + dj_cos_t[i] - di_sin_t;
		pt.y = cy + dj_sin_t[i] + di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
		{
			dtmp2 = real_coord(y,pt.y,pt.x);
			sum_n2ch += dtmp2;
			sum_n2ch2 += dtmp2*dtmp2;
			num_n2ch++;
			num++;
			if(flag)
			{
				sum_s += fabs(dtmp-dtmp2);
				num_s++;
			}
		}
		else
			num_out++;

	}

	//Find Intensyty in Channel 
	di = 0;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++)
	{
		pt.x = cx + dj_cos_t[i] + di_sin_t;
		pt.y = cy + dj_sin_t[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
		{
			derivative_energy += im_filter(y, rows, cols, t, floor(pt.y+0.5), floor(pt.x + 0.5));
			dtmp = real_coord(y,pt.y,pt.x);
			sum_2ch += dtmp;
			sum_2ch2 += dtmp*dtmp;
			num_2ch++;
			num++;
		}
		else
			num_out++;
	}

	// If enough of the object is out of the image, ignore object
	if(((double)num)/(double)(num_out+num)<OBJECT_IN_RATIO)
	{
		e[0] = 0;
		e[1] = 0;
		e[2] = 0;
		e[3] = 0;
		e[4] = 0;
		e[5] = 0;
		return INF;
	}

	num_nch += num_n1ch;
	sum_nch += sum_n1ch;

	sum_nch2 += sum_n1ch2;
	num_nch += num_n2ch;
	sum_nch += sum_n2ch;
	sum_nch2 += sum_n2ch2;

	num_ch += num_2ch;
	sum_ch += sum_2ch;
	sum_ch2 += sum_2ch2;

	sum_all += sum_nch;
	sum_all2 += sum_nch2;
	sum_all += sum_ch;
	sum_all2 += sum_ch2;

	/* Sample Mean and Sample STD */
	/* Non Channel 1 */
	if(num_n1ch){
		sum_n1ch = sum_n1ch/(double)num_n1ch;
		sum_n1ch2 = sum_n1ch2/(double)num_n1ch;
	}
	/* Non Channel 2 */
	if(num_n2ch){
		sum_n2ch = sum_n2ch/(double)num_n2ch;
		sum_n2ch2 = sum_n2ch2/(double)num_n2ch;
	}
	/*Channel 2 */
	if(num_2ch){
		sum_2ch = sum_2ch/(double)num_2ch;
		sum_2ch2 = sum_2ch2/(double)num_2ch;
	}

	/* Symmetry Error*/
	if(num_s){
		sum_s = sum_s/(double)num_s;
	}

	sum_2ch2 = sum_2ch2 - sum_2ch*sum_2ch; // variance of a channel 2

	sum_n1ch2 = sum_n1ch2 - sum_n1ch*sum_n1ch; // variance of the outside of a chanel 1
	sum_n2ch2 = sum_n2ch2 - sum_n2ch*sum_n2ch; // variance of the outside of a chanel 2

	if(num_ch){
		sum_ch = sum_ch/(double)num_ch;
		sum_ch2 = sum_ch2/(double)num_ch;
	}
	if(num_nch){
		sum_nch = sum_nch/(double)num_nch;
		sum_nch2 = sum_nch2/(double)num_nch;
	}
	if(num){
		sum_all = sum_all/(double)(num_ch+num_nch);
		sum_all2 = sum_all2/(double)(num_ch+num_nch);
	}

	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel
	sum_all2 = sum_all2 - sum_all*sum_all; // variance of the outside of a chanel
	if(num !=0)
	{
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		if((sum_ch2==0)&&(sum_nch2==0))
			t_test1 = fabs(sum_ch-sum_nch);
		else
			t_test1 = fabs(sum_ch-sum_nch)/(sqrt(sum_ch2/(double)num_ch+sum_nch2/(double)num_nch));
		//t_test2 = fabs(sum_n1ch-sum_n2ch)/(sqrt(sum_n1ch2/(double)num_n1ch+sum_n2ch2/(double)num_n2ch));
		

		t_test = t_test1/(2*fmax(SYM_TH,sum_s/sum_all));//*sqrt(sum_all2)
		//printf("T_test: %.2f Symmetry Error: %.2f \n", t_test,sum_s/sum_all);



		if(derivative_energy <= derivative_threshold)
			derivative_likely = 1- derivative_energy/derivative_threshold;
		else
			derivative_likely = exp(-(derivative_energy-derivative_threshold)/(derivative_energy))-1;


		if (sum_nch<= sum_ch)
		{
			if(t_test < error_th)
				likely = 1-t_test/error_th;
			else
				likely = exp(-(t_test-error_th)/(t_test))-1;
		}
		else
		{
			likely = 1;
		}

		#if DERIVATIVE_LIKELY
			likely = (0.95*likely + 0.05*derivative_likely)/2.0;
		#endif

		#if 0
			lineObj freeSeg;
			freeSeg.enda.x = cx + dj_cos_t[0] + di_sin_t;
			freeSeg.enda.y = cy + dj_sin_t[0] - di_cos_t; 
			int temp_coord = floor(2.0*l_5/STEP_DX + 0.5);
 			freeSeg.endb.x = cx + dj_cos_t[temp_coord] + di_sin_t;
			freeSeg.endb.y = cy + dj_sin_t[temp_coord] - di_cos_t;

			double **yimg;
			yimg = (double **)get_img(width, height, sizeof(unsigned char));
			for(int iii = 0; iii< height; iii++)
				for(int jjj=0; jjj< width; jjj++)
				{
					yimg[iii][jjj] = y[iii][jjj];
				}		 

			//display_only_one_double(yimg, height, width, &freeSeg);
			printf("std_channel:%.2f std_non_channel:%.2f \n",sum_ch2,sum_nch2);
			free_img((void **)yimg);
		#endif

		if (likely >= 0)
			likely = pow(likely,1.0/2.0);
		else
		{
			likely = -1*pow(-1*likely,1.0/2.0);
		}

	}
	else
	{
		likely = INF;
	}
	
	e[ 0] = (double)num_ch;
	e[ 1] = sum_ch;
	e[ 2] = sum_ch2;
	e[ 3] = (double)num_nch;
	e[ 4] = sum_nch;
	e[ 5] = sum_nch2;
	e[ 6] = (double)num_1ch;
	e[ 7] = sum_1ch;
	e[ 8] = sum_1ch2;
	e[ 9] = (double)num_2ch;
	e[10] = sum_2ch;
	e[11] = sum_2ch2;
	e[12] = (double)num_3ch;
	e[13] = sum_3ch;
	e[14] = sum_3ch2;
	e[15] = (double)num_n1ch;
	e[16] = sum_n1ch;
	e[17] = sum_n1ch2;
	e[18] = (double)num_n2ch;
	e[19] = sum_n2ch;
	e[20] = sum_n2ch2;
	e[21] = sum_s;
	e[22] = 0;
	e[23] = 0;
	e[24] = sum_all2;
	e[25] = t_test;
	e[26] = sum_s;
	e[27] = 0;

	patch[0][0] = derivative_energy;
	return likely;

}



double dataterm_line(int x1, int y1, int x2, int y2, double **img,double T)
{
    int dx, dy, inx, iny, e;
	double N1,N2;

	N1=N2=0;

	dx = x2 - x1;
	dy = y2 - y1;
	inx = dx > 0 ? 1 : -1;
	iny = dy > 0 ? 1 : -1;

	dx = ABS(dx);
	dy = ABS(dy);

	if(dx >= dy) {
	dy <<= 1;
	e = dy - dx;
	dx <<= 1;
	while (x1 != x2) {
		N2++;
		if (img[y1][x1] > T)
			N1++;
	//img[y1][x1] = color;
	if(e >= 0) {
	y1 += iny;
		e-= dx;
	}
	e += dy; x1 += inx;
	}
	} else {
	dx <<= 1;
	e = dx - dy;
	dy <<= 1;
	while (y1 != y2) {
	N2++;
	if (img[y1][x1] >T)
		N1++;
	if(e >= 0) {
	x1 += inx;
	e -= dy;
	}
	e += dx; y1 += iny;
	}
	}
	N2++;
	if (img[y1][x1] > T)
		N1++;

	return (N1/N2);
}


void K_means(double **img1, double *mu,double *cov, int height,int width,int k)
{
  int *hist;
  int *tag;
  int i,j,tmp,num[MAX_K];
  double max_pixel,min,sum;
   
  
  hist = (int *)malloc(256*sizeof(int)); //histogram of image
  tag = (int *)malloc(256*sizeof(int)); //tag of image
  
  
 // k = 18; //mean_k value
  max_pixel = 0;
  min = 255;

  for(i = 0;i<k;i++)
	  num[i] = 0;
  
  for ( i = 0; i < height; i++ )
  for ( j = 0; j < width; j++ ) {
	if (max_pixel < img1[i][j])
		max_pixel = img1[i][j];
	if (min > img1[i][j])
		min = img1[i][j];
  }
    for (i = 0;i<k;i++)
  {
	  mu[i] = min + (i+1)*(max_pixel-min)/(k+1);

  }

  for (i=0;i<256;i++)
	  {
		  hist[i] = 0;
		  tag[i] = 0;
  }

  for ( i = 0; i <height; i++ )
  for ( j = 0; j <width; j++ ) {
	  tmp = (int)img1[i][j];
	  hist[tmp] ++;
  }

  for (i=0;i<256;i++)
  {
	  if(hist[i]!=0)
	  {
		  min = 255;
		for (j = 0;j<k;j++)
		{
			if(ABS(i-mu[j])< min)
			{
				min = ABS(i-mu[j]);
				tmp = j;
			}
		}
		tag[i] = tmp;
	  }
  }

  for(i = 0;i<k;i++)
  {
	  tmp = 0;
	  sum = 0;
	  for(j = 0;j<256;j++)
	  {
		  if (tag[j] == i)
		  {
			  tmp += hist[j] * j; 
			  sum += hist[j];
		  }
	  }
	  mu[i] = (double)tmp/sum;
  }
 
    for (i=0;i<256;i++)
  {
	  if(hist[i]!=0)
	  {
		  min = 255;
		for (j = 0;j<k;j++)
		{
			if(ABS(i-mu[j])< min)
			{
				min = ABS(i-mu[j]);
				tmp = j;
			}
		}
		tag[i] = tmp;
	  }
  }


  for ( i = 0; i < height; i++ )
  for ( j = 0; j < width; j++ ) {
	   tmp = tag[(int)img1[i][j]];
	   num[tmp]++;
	  
  }

 /* for ( i = 0; i < height; i++ )
  for ( j = 0; j < width; j++ ) {
		tmp = tag[img1[i][j]];
	    mu[tmp] += (double)img1[i][j]/(double)num[tmp];
  }*/

  for(tmp = 0;tmp < k; tmp++)
  { double tmp_cov = 0;
	  for (i = 0;i< height;i++)
		 for (j = 0;j< width;j++)
		 {
			 if(tag[(int)img1[i][j]]== tmp)
			 { tmp_cov += (mu[tmp]-img1[i][j])*(mu[tmp]-img1[i][j]);
			 }
		  }

		 cov[tmp] = tmp_cov/num[tmp];


		
  }

  free((void *)hist);
  free((void *)tag);
  

}
