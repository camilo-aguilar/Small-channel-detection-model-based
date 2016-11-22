#include "QualityCandy.h"
//KDW #include "graph.h"
#include "tiff.h"
#include "allocate.h"
#include "util.h" //KDW


void likelyhoodMap (double **input_img, double **out_img_mpp, double **out_img_seg, double *mu, double *cov,double **patch, int L, int upper, int lower, double theta, int height, int width,int patch_len,int flag)
{
	int **indexmask=(int **)get_img(patch_len,patch_len,sizeof(int));
	printf("theta = %f\n",theta);
	for (int i = 0+THICK_W*2;i<height-THICK_W*2;i++)
	for (int j = 0+THICK_W*2; j < width-THICK_W*2;j++)
	{
		out_img_mpp[i][j] = GetLikelyhood(i,j,input_img,L,theta,upper,lower,patch, patch_len, height, width);
	 
#if 1
		if (out_img_mpp[i][j] > 0)
		{
		out_img_seg[i][j] = GetLikelyhood_seg(i,j,input_img,L,theta,upper,lower,patch,patch_len,mu,cov,indexmask,flag);
		if (out_img_seg[i][j]>0)
			out_img_seg[i][j] = 0;
		}
		//printf("d = %f,seg=%f,theta = %f\n",out_img_mpp[i][j],out_img_seg[i][j],theta);
#endif
	}

	free_img( (void**)indexmask );

}

#if 1 //KDW
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

#if 0
	//for test
	 FILE *fp;
  struct TIFF_img input_img;
  get_TIFF ( &input_img,  patch_len, patch_len, 'g' );
  for(int k = 0;k<patch_len;k++)
	  for(int l = 0;l<patch_len;l++)
		  input_img.mono[k][l] = patch[k][l];
   if ( ( fp = fopen ( "patch_test.tiff", "wb" ) ) == NULL ) {
    fprintf ( stderr, "cannot open file MBD_img.tif\n");
    exit ( 1 );
  }

  /* write image */
  if ( write_TIFF ( fp, &input_img ) ) {
    fprintf ( stderr, "error writing TIFF file" );
    exit ( 1 );
  }

  /* close image file */
  fclose ( fp );
	//end test
#endif
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
#else
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

	int x_mid = int((x_left + x_lower2)/2.0);
	int y_mid = int((y_left + y_lower2)/2.0);

	offset.x = (int)(x_mid - patch_len/2.0-1);
	offset.y = (int)(y_mid - patch_len/2.0-1);
	int offsetX = offset.x;
	int offsetY = offset.y;

	line(x-offsetX,y-offsetY,x2-offsetX,y2-offsetY,128,patch);
	line(x_upper-offsetX,y_upper-offsetY,x_upper2-offsetX,y_upper2-offsetY,100,patch);
	line(x_upper_down-offsetX,y_upper_down-offsetY,x_upper_down2-offsetX,y_upper_down2-offsetY,255,patch);

#if 0
	//for test
	 FILE *fp;
  struct TIFF_img input_img;
  get_TIFF ( &input_img,  patch_len, patch_len, 'g' );
  for(int k = 0;k<patch_len;k++)
	  for(int l = 0;l<patch_len;l++)
		  input_img.mono[k][l] = patch[k][l];
   if ( ( fp = fopen ( "patch_test.tiff", "wb" ) ) == NULL ) {
    fprintf ( stderr, "cannot open file MBD_img.tif\n");
    exit ( 1 );
  }

  /* write image */
  if ( write_TIFF ( fp, &input_img ) ) {
    fprintf ( stderr, "error writing TIFF file" );
    exit ( 1 );
  }

  /* close image file */
  fclose ( fp );
	//end test
#endif

	return offset;
}
#endif

double GetLikelyhood(int y, int x, double **input_img, int L, double theta, int upper_w, int lower_w,  double **patch,int patch_len, int height, int width)
{
	site offset;
	double mu1,mu2,mu3;
	double sum1,sum2,sum3,sum4,sum5,sum6;
	double N1,N2,N3;
	double cov1,cov2,cov3;

	sum3=sum4=sum1=sum2=sum5=sum6=N3=N1=N2 = 0;

	for (int i =0;i<patch_len;i++)
		for (int j = 0;j<patch_len;j++)
			patch[i][j] = 0;

	offset = drawRect(patch,patch_len,y,x,L,theta,upper_w,lower_w);

	for (int i = 0;i < patch_len;i++)
		for (int j = 0; j< patch_len;j++)
	{
		if (patch[i][j]== 128)
		{
			//output_img[i + offset.y][j+offset.x] = 128;
			if ((i + offset.y > 0 )&& (i + offset.y < height) && (j+offset.x > 0) &&(j+offset.x< width))
			{
			sum1+=input_img[i + offset.y][j+offset.x];
			sum3+=input_img[i + offset.y][j+offset.x]*input_img[i + offset.y][j+offset.x];
			N1++;
			}
		}
		if (patch[i][j]== 100)
		{
			//output_img[i + offset.y][j+offset.x] = 128;
			if ((i + offset.y > 0 )&& (i + offset.y < height) && (j+offset.x > 0) &&(j+offset.x< width))
			{
			sum5+=input_img[i + offset.y][j+offset.x];
			sum6+=input_img[i + offset.y][j+offset.x]*input_img[i + offset.y][j+offset.x];
			N3++;
			}
		}
		if (patch[i][j]== 255)
		{
			//output_img[i + offset.y][j+offset.x] = 255;
			if ((i + offset.y > 0 )&& (i + offset.y < height) && (j+offset.x > 0) &&(j+offset.x< width))
			{
			sum2+=input_img[i + offset.y][j+offset.x];
			sum4+=input_img[i + offset.y][j+offset.x]*input_img[i + offset.y][j+offset.x];
			N2++;
			}
		}
	}
	mu1 = sum1/N1;
	mu2 = sum2/N2;
	mu3 = sum5/N3;
	cov1 = sum3/N1-mu1*mu1;
	cov2=sum4/N2-mu2*mu2;
	cov3 = sum6/N3-mu3*mu3;

	double student_t1 = 10*ABS(mu1-mu2)/sqrt(cov1/N1+cov2/N2); 
	double student_t2 = 10*ABS(mu1-mu3)/sqrt(cov1/N1+cov3/N3); 
 
	double student_t = (student_t1+student_t2);
	if ( mu1 < mu2 || mu1 < mu3|| ABS(mu1-mu2)<5 || ABS(mu1-mu3)<5)//10)
		student_t = 0;
 
	if (student_t <50)//100)
		student_t = 0;
 
	return student_t;
}

#if 1
double GetLikelyhood_seg(int y, int x, double **input_img, int L, double theta, int upper_w, int lower_w,  double **patch,int patch_len,double *mu, double *cov, int **indexmask,int flag)
{
	site offset;

	for (int i =0;i<patch_len;i++)
		for (int j = 0;j<patch_len;j++)
			patch[i][j] = 0;

	offset = drawRect(patch,patch_len,y,x,L,theta,upper_w,lower_w);
#if 0
	typedef Graph<double,double,double> GraphType;
	GraphType *g = new GraphType(/*estimated # of nodes*/ patch_len*patch_len, /*estimated # of edges*/patch_len*patch_len*2); 
	GraphType *g_basic = new GraphType(/*estimated # of nodes*/ patch_len*patch_len, /*estimated # of edges*/patch_len*patch_len*2); 
	
	int index = 0;

	for (int i = 0;i < patch_len;i++)
		for (int j = 0; j< patch_len;j++)
	{
		if (patch[i][j]== 128)
		{ 
			if ((i + offset.y > 0 )&& (i + offset.y < 728) && (j+offset.x > 0) &&(j+offset.x< width))
			{

				g -> add_node();
				g_basic->add_node();

				int yy =i+offset.y;
				int xx = j+offset.x;
				indexmask[i][j] = index;

				double aa = (input_img[yy][xx]-mu[0])*(input_img[yy][xx]-mu[0])/(2*cov[0])+log(sqrt(2*_PI*cov[0]));
				double bb = (input_img[yy][xx]-mu[1])*(input_img[yy][xx]-mu[1])/(2*cov[1])+log(sqrt(2*_PI*cov[1]));
				// assume mu1,cov1 are for foreground, mu2,cov2 are for background
				g -> add_tweights( index,   /* capacities */  aa-SEG_OFFSET, bb+SEG_OFFSET );
				g_basic -> add_tweights( index,   /* capacities */  aa, bb);
				index++;
			}
		}
		if (patch[i][j]== 255 || patch[i][j] == 100)
		{
			if ((i + offset.y > 0 )&& (i + offset.y < 728) && (j+offset.x > 0) &&(j+offset.x< width))
			{
				g -> add_node();
				g_basic->add_node();

				int yy =i+offset.y;
				int xx = j+offset.x;
				indexmask[i][j] = index;

				double aa = (input_img[yy][xx]-mu[0])*(input_img[yy][xx]-mu[0])/(2*cov[0])+log(sqrt(2*_PI*cov[0]));
				double bb = (input_img[yy][xx]-mu[1])*(input_img[yy][xx]-mu[1])/(2*cov[1])+log(sqrt(2*_PI*cov[1]));
				// assume mu1,cov1 are for foreground, mu2,cov2 are for background
				g -> add_tweights( index,   /* capacities */  aa+SEG_OFFSET, bb-SEG_OFFSET );
				g_basic -> add_tweights( index,   /* capacities */  aa, bb);
				index++;
			}
			}
		}

		// check horizontal links
	for (int i = 0;i < patch_len;i++)
		for (int j = 0; j< patch_len-1;j++)
		{
				int yy =i+offset.y;
				int xx = j+offset.x;
			
				if ((yy > 0 )&& (yy < 728) && (xx > 0) &&(xx< width-1)&&patch[i][j]!=0&&patch[i][j+1]!=0)
				{
					if (patch[i][j]!=patch[i][j+1])
						 g->add_edge(indexmask[i][j],indexmask[i][j+1],SMALL_BETA,SMALL_BETA);
				
					else
						g->add_edge(indexmask[i][j],indexmask[i][j+1],BIG_BETA,BIG_BETA);

					g_basic->add_edge(indexmask[i][j],indexmask[i][j+1],BIG_BETA,BIG_BETA);
					
			}
		}
		//check vertical links
	for (int i = 0;i < patch_len-1;i++)
		for (int j = 0; j< patch_len;j++)
		{
				int yy =i+offset.y;
				int xx = j+offset.x;
			
				if ((yy > 0 )&& (yy < 728-1) && (xx > 0) &&(xx< width)&&patch[i][j]!=0&&patch[i+1][j]!=0)
				{
					if (patch[i][j]!=patch[i+1][j])
						 g->add_edge(indexmask[i][j],indexmask[i+1][j],SMALL_BETA,SMALL_BETA);
				
					else
						g->add_edge(indexmask[i][j],indexmask[i+1][j],BIG_BETA,BIG_BETA);

					g_basic->add_edge(indexmask[i][j],indexmask[i+1][j],BIG_BETA,BIG_BETA);
					
			}
		}
		double flow = g -> maxflow();
		double flow_basic = g_basic -> maxflow();
		delete g;
		delete g_basic;

		return (flow-flow_basic)/((upper_w+2*lower_w)*L);
		#endif
		return 0;
}
#endif

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

site drawRect_forDataterm(double **patch, int patch_len, int y, int x, int L, double theta, int upper_w)
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



	int x_mid = int((x_upper + x_upper_down2)/2.0);
	int y_mid = int((y_upper + y_upper_down2)/2.0);

	offset.x = (int)(x_mid - patch_len/2.0-1);
	offset.y = (int)(y_mid - patch_len/2.0-1);
	int offsetX = offset.x;
	int offsetY = offset.y;

	line(x_upper_down-offsetX,y_upper_down-offsetY,x_upper-offsetX,y_upper-offsetY,128,patch);
	line(x_upper_down-offsetX,y_upper_down-offsetY,x_upper_down2-offsetX,y_upper_down2-offsetY,128,patch);
	line(x_upper-offsetX,y_upper-offsetY,x_upper2-offsetX,y_upper2-offsetY,128,patch);
	line(x_upper_down2-offsetX,y_upper_down2-offsetY,x_upper2-offsetX,y_upper2-offsetY,128,patch);

	subfill(patch,int((x_upper_down2+x_upper)/2.0)-offsetX,int((y_upper_down2+y_upper)/2.0)-offsetY,128);

	return offset;
}

#if 0 //KDW
double dataterm_rec(int **img_seg, lineObj *line,double **patch, int patch_len, int line_w,double ***Matrix, double prior_pen1, double prior_pen2, int height,int width)
{
	site enda = line->enda;
	site endb = line->endb;
	double theta = line->theta;
	int len = line->len;
	site start_end;
	
//	if((enda.x == 313 && enda.y == 90 && endb.x == 324 && endb.y == 80)||(endb.x== 313 && endb.y == 90&& enda.x == 324 && enda.y == 80))
//		int test = 0;
	int xx = enda.x + (int)((double)len*cos(theta));
	int yy = enda.y - (int)((double)len*sin(theta));

	if(DIST(xx,yy,endb.x,endb.y)< 3*3)
		start_end = enda;
	else
	{
		xx = endb.x + (int)((double)len*cos(theta));
		yy = endb.y - (int)((double)len*sin(theta));
		if(DIST(xx,yy,enda.x,enda.y)< 3*3)
			start_end = endb;
		else
			{
				printf("find start_end wrong\n");
				return 0;
		}
	}
	for (int i =0;i<patch_len;i++)
	for (int j = 0;j<patch_len;j++)
		patch[i][j] = 0;

	site offset = drawRect(patch,patch_len, start_end.y, start_end.x, len, theta, line_w, THIN_W);

	double sum1,sum2;
	int N1,N2,N3;
//KDW	double likely;

	sum1 = sum2 =0;
	N1=N2=N3=0;
	double penalty = 0;
	double penalty1,penalty2,penalty3; 
	penalty1= penalty2 = penalty3 = 0;
#if 1
	for (int i = 0;i < patch_len;i++)
		for (int j = 0; j< patch_len;j++)
	{
		if (patch[i][j]== 128)
		{ 
			if ((i + offset.y > 0 )&& (i + offset.y < 728) && (j+offset.x > 0) &&(j+offset.x< width))
			{
				int yy =i+offset.y;
				int xx = j+offset.x;
				sum+= 0.4*img_mpp_l[yy][xx]+0.6*img_seg_l[yy][xx];
			//	if(img_mpp_l[yy][xx] > 128 && img_seg_l[yy][xx] > 128)
			//		sum++;
				N++;
			}
		}
		}
#endif
	for (int i = 0;i < patch_len;i++)
		for (int j = 0; j< patch_len;j++)
	{
		int yy = i+offset.y;
		int xx = j+offset.x;

		if (patch[i][j]== 128)
		{
			if ((yy > 0 )&& (yy < height) && (xx > 0) &&(xx< width))
			{
				if(img_seg[yy][xx] == 2)
				{
					N1++;
					penalty1-=PEN_1;
				}
				else if (img_seg[yy][xx] == 1)
				{
					N2++;
					
					penalty2 += MIN(PEN_2,Matrix[2][yy][xx]-Matrix[1][yy][xx]);
				}
				else if (img_seg[yy][xx] == 0)
				{
					N3++;
					penalty2 += MIN(PEN_3,Matrix[2][yy][xx]-Matrix[0][yy][xx]);
				}
			}
		}
		if (patch[i][j]== 100 || patch[i][j] == 255)
		{

			if ((yy > 0 )&& (yy < height) && (xx > 0) &&(xx< width))
			{
				if(img_seg[yy][xx] == 2)
				{
					N3++;
					penalty3 += MIN(PEN_1,Matrix[1][yy][xx]-Matrix[2][yy][xx]);
					//sum2+= Matrix[0][yy][xx]-Matrix[2][yy][xx];
				}
				else if (img_seg[yy][xx] == 1)
				{
				//	N2--;
				//	sum2+= Matrix[0][yy][xx]-Matrix[1][yy][xx];
				}
			}
		}

	}
//	double penalty = N2*prior_pen1+N3*prior_pen2;
//	double total = sum1+sum2+penalty;
		penalty = (penalty1/**1.5*/+penalty2+penalty3)/2.0;
	if(penalty > - 4)
		penalty = 0;
	return (penalty);


}
#else
#define STEP_DX					0.5
#define STEP_DY					0.5
#define SYM_TH	0.11//0.09
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
	double likely;
	double t_test, t_test1;//, t_test2; 
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

	num = 0;

	//if((line->enda.x==157)&&(line->enda.y==125)&&(line->endb.x==151)&&(line->endb.y==125))
	//	printf("''");
	//if((line->enda.x==160)&&(line->enda.y==122)&&(line->endb.x==149)&&(line->endb.y==123))
	//	printf("''");

	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++){
		dj_cos_t[i] = dj*cos_t;
		dj_sin_t[i] = dj*sin_t;
	}
	di = -w_5;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++){
		flag = 0;
		// non channel 1
		pt.x = cx + dj_cos_t[i] + di_sin_t;
		pt.y = cy + dj_sin_t[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
//			y[(int)pt.y][(int)pt.x] = 0;
			sum_n1ch += dtmp;
			sum_n1ch2 += dtmp*dtmp;
			num_n1ch++;
			num++;
			flag = 1;
		}
		else
			num_out++;
		// non channel 2
		pt.x = cx + dj_cos_t[i] - di_sin_t;
		pt.y = cy + dj_sin_t[i] + di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp2 = real_coord(y,pt.y,pt.x);
//			y[(int)pt.y][(int)pt.x] = 0;
			sum_n2ch += dtmp2;
			sum_n2ch2 += dtmp2*dtmp2;
			num_n2ch++;
			num++;
			if(flag){
				sum_s += fabs(dtmp-dtmp2);
				num_s++;
			}
		}
		else
			num_out++;

	}
#if 0
	// channel 1
	di = -STEP_DY;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++){
		pt.x = cx + dj_cos_t[i] + di_sin_t;
		pt.y = cy + dj_sin_t[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_1ch += dtmp;
			sum_1ch2 += dtmp*dtmp;
			num_1ch++;
			num++;
		}
		else
			num_out++;
	}
#endif
	// channel 2
	di = 0;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++){
		pt.x = cx + dj_cos_t[i] + di_sin_t;
		pt.y = cy + dj_sin_t[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_2ch += dtmp;
			sum_2ch2 += dtmp*dtmp;
			num_2ch++;
			num++;
		}
		else
			num_out++;
	}
#if 0
	// channel 3
	di = STEP_DY;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5, i = 0; dj<=l_5; dj += STEP_DX, i++){
		pt.x = cx + dj_cos_t[i] + di_sin_t;
		pt.y = cy + dj_sin_t[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_3ch += dtmp;
			sum_3ch2 += dtmp*dtmp;
			num_3ch++;
			num++;
		}
		else
			num_out++;
	}
#endif
	// Calculate Distance
	if(((double)num)/(double)(num_out+num)<OBJECT_IN_RATIO){
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
//	num_ch += num_1ch;
//	sum_ch += sum_1ch;
//	sum_ch2 += sum_1ch2;
	num_ch += num_2ch;
	sum_ch += sum_2ch;
	sum_ch2 += sum_2ch2;
//	num_ch += num_3ch;
//	sum_ch += sum_3ch;
//	sum_ch2 += sum_3ch2;
	sum_all += sum_nch;
	sum_all2 += sum_nch2;
	sum_all += sum_ch;
	sum_all2 += sum_ch2;
#if 1
	if(num_n1ch){
		sum_n1ch = sum_n1ch/(double)num_n1ch;
		sum_n1ch2 = sum_n1ch2/(double)num_n1ch;
	}
	if(num_n2ch){
		sum_n2ch = sum_n2ch/(double)num_n2ch;
		sum_n2ch2 = sum_n2ch2/(double)num_n2ch;
	}
	if(num_1ch){
		sum_1ch = sum_1ch/(double)num_1ch;
		sum_1ch2 = sum_1ch2/(double)num_1ch;
	}
	if(num_2ch){
		sum_2ch = sum_2ch/(double)num_2ch;
		sum_2ch2 = sum_2ch2/(double)num_2ch;
	}
	if(num_3ch){
		sum_3ch = sum_3ch/(double)num_3ch;
		sum_3ch2 = sum_3ch2/(double)num_3ch;
	}
	if(num_s){
		sum_s = sum_s/(double)num_s;
	}
	sum_1ch2 = sum_1ch2 - sum_1ch*sum_1ch; // variance of a channel
	sum_2ch2 = sum_2ch2 - sum_2ch*sum_2ch; // variance of a channel
	sum_3ch2 = sum_3ch2 - sum_3ch*sum_3ch; // variance of a channel
	sum_n1ch2 = sum_n1ch2 - sum_n1ch*sum_n1ch; // variance of the outside of a chanel
	sum_n2ch2 = sum_n2ch2 - sum_n2ch*sum_n2ch; // variance of the outside of a chanel
#endif
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
	if(num !=0){
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		if((sum_ch2==0)&&(sum_nch2==0))
			t_test1 = fabs(sum_ch-sum_nch);
		else
			t_test1 = fabs(sum_ch-sum_nch)/(sqrt(sum_ch2/(double)num_ch+sum_nch2/(double)num_nch));
		//t_test2 = fabs(sum_n1ch-sum_n2ch)/(sqrt(sum_n1ch2/(double)num_n1ch+sum_n2ch2/(double)num_n2ch));
		t_test = t_test1/max(SYM_TH,sum_s/sum_all)*SYM_TH;//*sqrt(sum_all2)
/* Huixi
	if (mu1>mu2)
	{
		if (B<thredh)
			likely = 1-B/thredh;
		else
			likely = exp(-(B-thredh)/(3*B))-1;
	}
	else
		likely = 1;

	if (likely >= 0)
		return(pow(likely,1.0/3.0));
	else
		return(-1*pow(-1*likely,1.0/3.0));

*/
		if ((sum_nch>sum_ch)&&(sum_all2>100))//150
		{
			if(t_test < error_th)
					likely = 1-t_test/error_th;
				else
					likely = exp(-(t_test-error_th)/(3.*t_test))-1;
		}
		else
			likely = 1;

		
		if (likely >= 0)
			//likely = pow(likely,1.0/4.0);
			likely = pow(likely,1.0/2.0);
			//likely = 0;
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/2.0);
			

	}
	else{
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

	return likely;

}
#endif

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
#if 0
double dataterm(double **input_img,site enda, site endb, double theta, int len,int theta_num,double **patch, int patch_len)
{
	int x,y;
	if (enda.y > endb.y)
	{
		x = enda.x;
		y = enda.y;
	}
	else
	{
		x = endb.x;
		y = endb.y;
	}
	
	double likely;
	if (theta_num >= 15)
		likely= GetLikelyhood(y,x,input_img,len,theta,THICK_W,THIN_W,patch, patch_len);
	else
		likely= GetLikelyhood(y,x,input_img,len,theta,THIN_W,THICK_W,patch, patch_len);

	if (likely > 80)
		return exp(likely/10.0);
	else
		return exp(-20.0);
}
#endif
#if 1
	double dataterm(double **input_img,site enda, site endb)
{
	double T = 100;
	double r = dataterm_line(enda.x, enda.y, endb.x, endb.y,input_img,T);
	if (r > 0.7)
		return exp(7.0*2*r);
	else
		return exp(-20.0);
}

	double dataterm_double(double **input_img,site enda, site endb)
{
	double T = 100;
	double r = dataterm_line(enda.x, enda.y, endb.x, endb.y,input_img,T);
	if (r > 0.7)
		return exp(7.0*2*r);
	else
		return exp(-20.0);
}
#endif

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
