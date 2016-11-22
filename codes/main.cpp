#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

//#include "cv.h"
//#include "cxcore.h"
//#include "highgui.h"
#include "allocate.h"
#include "random.h"
#include "tiff.h"
#include "neck.h"
#include "ndmulti_e.h"
#include "em.h"
//#include "synthetic.h"
#include "QualityCandy.h"




double calculate_PMP(unsigned char **xt, unsigned char **gt, int classes, int rows, int cols)
{
	unsigned char xt255;
	int i, j;
	double misclassed = 0;

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			xt255 = (unsigned char)xt[i][j] * 255 / (classes - 1);
			if(xt255!=gt[i][j]) misclassed += 1.0;
		}
	misclassed = misclassed*100/(cols*rows); // percentage of misclassified pixels
	return misclassed;
}


int compare_fn0(const void *a, const void *b)
{
	if(*(unsigned char *)a > *(unsigned char *)b) 
		return -1;
	else if(*(unsigned char *)a == *(unsigned char *)b) 
		return 0;
	else 
		return 1;
}


unsigned char median_filter(unsigned char *median, int n)
{
	qsort((void *)median, n, sizeof(unsigned char), compare_fn0);
	return median[(n-1)/2];
}



void median_filter2D(unsigned char **y, unsigned char **yf, int filter_size, int rows, int cols, int enable)
{
	unsigned char **y2,**y3, median[100];
	int n,i,j,ii,jj,fs_5;

	if(enable){
		y2 = (unsigned char **)get_img(cols+filter_size, rows+filter_size, sizeof(unsigned char));
		y3 = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));

		fs_5 = (filter_size-1)/2;
		for (i = 0; i < rows+filter_size; i++)
			for (j = 0; j < cols+filter_size; j++){
				if (i<fs_5)
					if(j<fs_5)
						y2[i][j] = y[0][0];
					else if(j>=cols+fs_5)
						y2[i][j] = y[0][cols-1];
					else
						y2[i][j] = y[0][j-fs_5];
				else if(i>=rows+fs_5) 
					if(j<fs_5)
						y2[i][j] = y[rows-1][0];
					else if(j>=cols+fs_5)
						y2[i][j] = y[rows-1][cols-1];
					else
						y2[i][j] = y[rows-1][j-fs_5];
				else
					if(j<fs_5)
						y2[i][j] = y[i-fs_5][0];
					else if(j>=cols+fs_5)
						y2[i][j] = y[i-fs_5][cols-1];
					else
						y2[i][j] = y[i-fs_5][j-fs_5];
			}
		for (i = fs_5; i < rows+fs_5; i++)
			for (j = fs_5; j < cols+fs_5; j++){
				n = 0;
				for (ii = -fs_5; ii <= fs_5; ii++)
					for (jj = -fs_5; jj <= fs_5; jj++){
						median[n] = y2[i+ii][j+jj];
						n++;
					}
				y3[i-fs_5][j-fs_5] = median_filter(median,n);
			}

		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++)
				yf[i][j] = y3[i][j];

		free_img((void **)y2);
		free_img((void **)y3);
	}
	else{
		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++)
				yf[i][j] = y[i][j];
	}
}

int QuilityCandyInterface(unsigned char **yimg, double **lm, double *mean, double *vari, 
			double variance, NeckDent *mp, MPP_Parameters mpp, double *total_e, int *mp_num, 
			int cols, int rows)
{
  FILE *fp;
  struct TIFF_img input_img, MBD_img;
  double **img1,****img2,****img3,**img_out,**patch;
  int **img_seg;
  int32_t i,j,pixel,pixel2;

  input_img.height = rows;
  input_img.width = cols;
  int patch_len = 100;
  
  /* Allocate image of double precision floats */
  img1 = (double **)get_img(input_img.width,input_img.height,sizeof(double));
  img_out = (double **)get_img(input_img.width,input_img.height,sizeof(double));
  img_seg = (int **)get_img(input_img.width,input_img.height,sizeof(int));

  img2 = (double ****)malloc(5*sizeof(double ***));
  img3 = (double ****)malloc(5*sizeof(double ***));

  patch = (double **)get_img(patch_len,patch_len,sizeof(double));

  /* copy green component to double array */
  for ( i = 0; i < input_img.height; i++ )
  {  
	for ( j = 0; j < input_img.width; j++ ) 
	{
  		img1[i][j] = (double)yimg[i][j];
		img_out[i][j] = 0;
		img_seg[i][j]=0;
	 }
  }

  srandom2(1);
  int iter = mpp.iter_num; //1000000;//100000000;
  double T = mpp.T0;//KDW 1/10.0;

   char a[30];

   int  k = 3;


	


   double *mu = (double *)malloc(MAX_K*sizeof(double));
   double *cov = (double *)malloc(MAX_K*sizeof(double));
   //K_means(img1,mu,cov,input_img.height,input_img.width,3);

   get_TIFF ( &MBD_img, input_img.height, input_img.width, 'g' );

  (*mp_num) = Candy_Model(img1,lm,img2,img3,img_out,img_seg,input_img.height,input_img.width,
		T,iter,patch,patch_len,mu,cov,3,mp,mpp);
	

  for ( i = 0; i < input_img.height; i++ )
  {
 	 for ( j = 0; j < input_img.width; j++ ) 
 	 {
    	pixel = (int32_t)img_out[i][j];
		/*mask*/
    	if(pixel>255)
    	{
      		MBD_img.mono[i][j] = 255;
    	}
    	else
    	{
      		if(pixel<0) 
      			MBD_img.mono[i][j] = 0;
      		else 
      			MBD_img.mono[i][j] = pixel;
    	}
  	}
   }
	if(1)
	{
	  char aa[30];
	  sprintf(aa,"channels_out.tiff");
	  /* open image file */
	  if ( ( fp = fopen ( aa, "wb" ) ) == NULL ) {
	    fprintf ( stderr, "cannot open file MBD_img.tif\n");
	    exit ( 1 );
	  }

	  /* write image */
	  if ( write_TIFF ( fp, &MBD_img ) ) {
	    fprintf ( stderr, "error writing TIFF file" );
	    exit ( 1 );
	  }

	  /* close image file */
	  fclose ( fp );
	}
  /* de-allocate space which was used for the images */
  free_TIFF(&(MBD_img));

  free_img( (void**)img1 );
  free(img2);
  free(img3);
  free_img( (void**)img_out ); 
  free_img((void**)img_seg);
  free_img((void**)patch);
  free((void*)mu);
  free((void*)cov);

  return  (*mp_num);
}

void difference_image(unsigned char **img1, unsigned char  **img2,int height, int width, char* name);


/*
MAIN
*/
int main( int argc , char** argv)
{
	printf(">>>>>>>>>>>>>>>>Start of the Program<<<<<<<<<<<<<\n");
	struct TIFF_img input_img, input_gt_img, output_img, output_color_img;
	unsigned char **yimg = NULL, **yfiltered, **laplacian, **lm_img;
	double **lm, **beta_dimg[2];
	char infileName[1024], outfileName[1024];
	char segfileName[1024], outfilePrefix[1024];
	char gtfileName[1024];
	char win_name[256] = "input data";
	char win_name2[256] = "output data";
	FILE *fp;
	//IplImage *image = 0, *image2 = 0; //NOTE
	int height, width, channel=3;
	int i, j, ii, jj, rows, cols;
	double mu[CLASSES], mean, vari[CLASSES], variance, min, max;
	int cnt[CLASSES];
	// for EM/MPM
	double beta[MAX_CLASSES], gamma[MAX_CLASSES];
	int mpmiter = 10, emiter = 30, classes = 2;
	unsigned char **xt, **gt;  /* output : entropy image */
	double **blur, sigma = 0., dsum, sum[CLASSES], di, dj, misclassed;
	int blur_size = 5, enable_blur = 0;
	int run_emmpm = 0;
	NeckDent mp[MAX_MKPNT_NUM];
	int np_num;
	char filename[1024];
	double stdev;
	//int syn_ch_gen = 0, gen_syn_img_en = 0, syn_num = SYN_CH_NUM;
	//NeckDent_Syn syn_ch[SYN_CH_NUM];
	double total_e[MAX_MPP_ITER_NUM+1];
	int	   mp_num[MAX_MPP_ITER_NUM+1];
	MPP_Parameters mpp;
	double avg_width, avg_amp, avg_error, avg_se, avg_U, dtmp;

	int k, m, found;
	double running_time;
	clock_t start_time, end_time;


	char  PARAM_IMAGE_NAME[50];							//1  Name
	int   PARAM_TESTTYPE = 0;		     				//2  0:Test single object,  1:Run normal code 
	
	int   PARAM_OPTIMIZATION_TYPE = 2;	    			//3  0:RJMCMC, 1:Multiple Birth and Death	2:RJMCMC Quality Candy model
	
	int   PARAM_FIXED_PARAMETERS = -1;      			//4  0:synthetic,	1:BSEImage free angle, 2:BLEImage fixed angle 3:slice_crop
	float PARAM_LAMBDA_RJMCMC = 30000;	     			//5  Lambda(RJMCMC)  (intensity = lambda*vk), alm(MBND)0.5
	float PARAM_T0 = 0.08;			      				//6  T0 (big more death)
	float PARAM__BETA_MPP = 14;		      				//7  Betampp(1)
	float PARAM_LAMBDA_l = 0;	            			//8	 Lambda_l		
	float PARAM_LAMBDA_int = 0.005;		      			//9  Lambda_int
	
	int   PARAM_CHANNEL_TYPES = 0;               		//10 0:simple channel	1:necking only	2: denting only 3:necking and denting
	
	float PARAM_ERROR_THRESHOLD = 6;	      			//11 error_th	13
	float PARAM_AMPLITUDE_THRESHOLD = 21.0;	      		//12 Amp_th 	
	float PARAM_GAUSSIAN_TAU = 17;		      			//13 Gaussina_tau 	(RJMCMC)"20" (tau,lambda_a,error_th) = (25,27,14)
	float PARAM_DECREASE_COEFFICIENT = 0.99999;   		//14 Decrease coeff = 0.999999
	int   PARAM_ITERATIONS = 3000000;		      			//15 Iterations
	float PARAM_BETTA_FOR_MPP = 2.7;	      			//16 Beta for MPP
	
	/* Parameter Parsing */
	if(1)
	{ 
		sprintf(infileName, "%s.tiff",argv[1]);
		sprintf(segfileName, "%s_seg.tiff",argv[1]);
		sprintf(gtfileName, "%s_gt.tiff",argv[1]);
		sprintf(outfilePrefix, "%s",argv[1]);
		
		mpp.test = PARAM_TESTTYPE;							//2
		mpp.optimization_type = PARAM_OPTIMIZATION_TYPE;	//3
		mpp.fixed_param_num = PARAM_FIXED_PARAMETERS;		//4
		mpp.alm = PARAM_LAMBDA_RJMCMC;
		mpp.b_zero = mpp.alm;
		mpp.T0 = PARAM_T0;
		mpp.betampp = PARAM__BETA_MPP;
		mpp.lambda_l = PARAM_LAMBDA_l;	// 0.12
		mpp.lambda_int = PARAM_LAMBDA_int;
		mpp.nd_type_num = PARAM_CHANNEL_TYPES;	// necking only = 1, 
		mpp.error_th = PARAM_ERROR_THRESHOLD;
		mpp.amp_th = PARAM_AMPLITUDE_THRESHOLD;
		mpp.gaussian_tau = PARAM_GAUSSIAN_TAU;
		mpp.de_coeff = PARAM_DECREASE_COEFFICIENT;
		mpp.iter_num = PARAM_ITERATIONS; // override 
		beta[0] = PARAM_BETTA_FOR_MPP;
		beta[1] = PARAM_BETTA_FOR_MPP;
	}

	
	
	mpp.hard_repulsion	= 4;	// 4(RJMCMC) 4
	mpp.gaussian_tau	= 15;	// 6.(MBND) 20.(RJMCMC)  
	
	
	mpp.widthmin		= 4;	// 6 6 (MBND) 9(RJMCMC)	8	9	9
	mpp.widthmax		= 11;	// 10 12 (MBND) 14(RJMCMC)	14	18	14
	mpp.lengthmin		= 4;	// 8 10 (MBND) 11(RJMCMC)	10	5	5 // should be smaller than 10
	mpp.lengthmax		= 8;	// 20 30 (MBND) 34(RJMCMC)	32	25	32
	
	//Denting and Necking Channels
	mpp.symmetry_th		= 0;	// 8
	mpp.lambda_s		= 0;	//Symmetry Potential
	mpp.lambda_nc		= 0;
	
		
	//Quality Candy
	mpp.gamma_d			= 34;//34; 
	mpp.w_eo			= 8;//1
	mpp.w_f				= 6; //1
	mpp.w_s				= 30;//1
	mpp.w_d				= -7;//1
	mpp.w_io			= 1;//1
	
	if(mpp.w_f + mpp.w_s < mpp.gamma_d)
	{
		printf("Error in Parameters. Contraint 1 is not met \n");
		exit(1);
	}
	/*
	if(2*mmp.w_f + 2*mpp.w_s < mpp.gamma_d)
	{
		printf("Error in Parameters. Contraint 1 is not met \n")
		exit(1);
	}
	if(mmp.w_f + mpp.w_s < mpp.gamma_d)
	{
		printf("Error in Parameters. Contraint 1 is not met \n")
		exit(1);
	}
	if(mmp.w_f + mpp.w_s < mpp.gamma_d)
	{
		printf("Error in Parameters. Contraint 1 is not met \n")
		exit(1);
	}
	if(mmp.w_f + mpp.w_s < mpp.gamma_d)
	{
		printf("Error in Parameters. Contraint 1 is not met \n")
		exit(1);
	}
	*/
	
	//Denting
	mpp.dent_l_w_ratio	= 1.3; // 1.25 (MBND) 1.3(RJMCMC)1.3
	
	//Multiple Birth and Death
	mpp.p_birth			= 0.3;
	mpp.p_death			= 0.3;
	mpp.p_translation	= 0.2;
	mpp.p_dilation		= 0.2;
	mpp.p_rotation		= 0.08;
	mpp.p_switching		= 0.;
 
 
	mpp.alpha			= 1;
	mpp.sigma			= 0.9;
	mpp.blur_size		= 7;
	
	//Very Important
	mpp.lambda_e = 1;			//Object Potential
	mpp.lambda_dc= 3;			//Discontinuity Potential 

	
	mpp.lambda_a = 0.5;
	mpp.length_th = mpp.lengthmin;

	readseed();
	gamma[0] = 0.;
	gamma[1] = 0.;
	sigma = 0.58;
	
	
	start_time=clock();
	/*Read Image */	
	/*NOTE: IMAGE IS CURRENTLY BEING INVERTED */
	if(1)
	{
	
		if ((fp = fopen(infileName, "rb")) == NULL) 
		{
			printf("Cannot open file %s\n", infileName);
			exit(1);
		}
		else 
		{
			if (read_TIFF(fp, &input_img))
			{
				printf("Error reading file %s\n", infileName);
				exit(1);
			}
			cols = input_img.width;
			rows = input_img.height;
			/* close image file */
			fclose(fp);
			yimg = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
			/* check the type of image data */
			if (input_img.TIFF_type != 'g') 
			{
				printf("Error:  Image must be grayscale.\n");
				exit(1);
			}
			for (i = 0; i < rows; i++)
			{
				for (j = 0; j < cols; j++)
				{
					yimg[i][j] = input_img.mono[i][j];
				}
			}
		}
	}

	/*Vectorize Image*/
	if(1)
	{
		
		mpp.vk = (double)rows*cols;
		mpp.delta = mpp.vk;
		printf ("cols = %d, rows = %d\n",cols, rows);
	}


	/* Array allocation */
	if(1)
	{
		
		laplacian = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		yfiltered = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		gt = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		xt = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		lm = (double **)get_img(cols, rows, sizeof(double));
		beta_dimg[0] = (double **)get_img(cols, rows, sizeof(double));
		beta_dimg[1] = (double **)get_img(cols, rows, sizeof(double));
		lm_img = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		get_TIFF(&output_img, rows, cols, 'g');
		get_TIFF(&output_color_img, rows, cols, 'c');
	}


	/* Copy and Paste the values of yimg to yfiltered */
	median_filter2D(yimg, yfiltered, 3,rows, cols, 0);
	
	
	/* Read original EMMPM result iamge */
	if(1)
	{
		if ((fp = fopen(segfileName, "rb")) == NULL)
		{
			printf("Cannot open file %s\n", segfileName);
			printf("Execute EM/MPM Segmentation..\n");
			run_emmpm = 1;
		}
		else
		{
			if (read_TIFF(fp, &input_gt_img))
			{/* read image */
				printf("Error reading file %s\n", segfileName);
				printf("Execute EM/MPM Segmentation..\n");
				run_emmpm = 1;
			}
			else
			{
				for (i = 0; i < rows; i++)
					for (j = 0; j < cols; j++)
					{
						// segmentation result file to use as lebesgue measure
						xt[i][j] = (int)(input_gt_img.mono[i][j]*(classes - 1)/255);
						
					}
				run_emmpm = 0;
			}
		}
	} 

	/* Calculate Blur Matrix. Innecesary */
	if(1)
	{
		blur = (double **)get_img(blur_size, blur_size, sizeof(double));
		if(sigma == 0.)
			enable_blur = 0;
		else
		{
			enable_blur = 1;
			dsum = 0.;
			for (i = 0; i < blur_size; i++)
			{
				for (j = 0; j < blur_size; j++)
				{
					di = i-(blur_size-1)/2.;
					dj = j-(blur_size-1)/2.;
					blur[i][j] = exp(-sqrt(di*di+dj*dj)/sigma);
					dsum += blur[i][j];
				}
			}
			for (i = 0; i < blur_size; i++)
			{
				for (j = 0; j < blur_size; j++)
				{
					blur[i][j] = blur[i][j]/dsum;
				}

			}
		}
	}

	/*if run em/mpm */
	if(1)
	{
		//*****************************************************************************
		//		EM/MPM for bith map and learn mu0 vari0 mu1 vari1
		//*****************************************************************************
		// initialize for EM/MPM
		//max_entropy = log10(classes)/log(2);
		if (run_emmpm)
		{
			// calculate blurring matrix 
			/*	Aexp(-x/sigma)											*/
			/*	A(amplitude), x(distance from center), sigma(width)		*/
			

			start_time=clock();
			printf("start time is:%1.0f ms\n",(double)(start_time)/CLOCKS_PER_SEC*100);
			
			#ifdef EMMPM_ONLY
			
				beta[0] = 2.9;
				beta[1] = beta[0];
				beta[2] = beta[1];
				gamma[0] = 0;
				gamma[1] = gamma[0];
				gamma[2] = gamma[1];
				emiter = 10;
				mpmiter = 10;
				classes = 2;
				enable_blur = 0;
				sprintf(segfileName, "%s_seg8nn_beta%1.1f.tiff",argv[1],beta[0]);
			
			#endif
			enable_blur = 0;
			emmpm(yimg, xt, beta, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);
			end_time=clock();
			printf("end time is:%1.0f ms\n",(double)(end_time)/CLOCKS_PER_SEC*100);
			printf("Running time is:%1.0f ms\n",(double)(end_time-start_time)/CLOCKS_PER_SEC*100);
			for (i=0; i<rows; i++)
				for (j=0; j<cols; j++)
					output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
			
			if ((fp = fopen(segfileName, "wb")) == NULL ) {
				printf("Cannot open file %s\n", outfileName);
				exit(1);
			}
			if (write_TIFF(fp, &output_img)) {
				printf("Error writing TIFF file %s\n", outfileName);
				exit(1);
			}
			fclose(fp);
			
			// PMP
			//misclassed = calculate_PMP(xt, gt, classes, rows, cols);
			//sprintf(outfileName, "%s_s%1.2f_p%1.2f.txt",outfilePrefix, sigma, misclassed);
			//if ((fp = fopen(outfileName, "wb")) == NULL ) {
			//	printf("Cannot open file %s\n", outfileName);
			//	exit(1);
			//}
			//fprintf(fp,"parameters\n image = %s\n beta0 = %1.2f\n beta1 = %1.2f\n sigma = %1.2f\n pmp = %1.2f\n ",
			//		segfileName, beta[0], beta[1], sigma, misclassed);
			//fclose(fp);
		} 
	}
	
	for(i=0; i < input_img.height; i++)
	{
		for(j=0; j< input_img.width; j++)
		{
			gt[i][j] = xt[i][j];
		}
	}
	
	/* Calculate mean and variance for foreground and background*/
	if(1)
	{
		dsum = 0.;
		for (i=0; i<CLASSES; i++)
		{
			sum[i] = 0.;
			cnt[i] = 0;
		}
		
		/*Calculate the sum and count for each class*/	
		for (i=0; i<rows; i++)
		{
			for (j=0; j<cols; j++)
			{
				sum[xt[i][j]] += yimg[i][j];
				cnt[xt[i][j]] ++;
				dsum += yimg[i][j];
			}
		}
		
		
		mean = dsum/(double)(rows*cols);
		dsum = 0.;

		/*Calculate mu for each class */
		for (i=0; i<CLASSES; i++)
		{
			mu[i] = sum[i]/(double)cnt[i];
			mpp.mean[i] = mu[i];
			sum[i] = 0.;
		}

		/* Sum to calculate variance */
		for (i=0; i<rows; i++)
		{
			for (j=0; j<cols; j++)
			{
				sum[xt[i][j]] += ((double)yimg[i][j]-mu[xt[i][j]])*((double)yimg[i][j]-mu[xt[i][j]]);
				dsum += ((double)yimg[i][j]-mean)*((double)yimg[i][j]-mean);
			}
		}
		variance = dsum/(double)(rows*cols);
		
		/*Variance for each class*/
		for (i=0; i<CLASSES; i++)
		{
			vari[i] = sum[i]/(double)cnt[i];
			mpp.vari[i] = vari[i];
			printf ("mu[%d] = %1.3f, vari[%d] = %1.3f, variance = %1.3f\n",i,mu[i],i,vari[i],variance);
		}
	}

	#define BIRTH_MAP_R		3

	/*Create Birthmap*/
	if(1)
	{
		min = 10000;
		max = 0;
		for (i=0; i<rows; i++)
		{
			for (j=0; j<cols; j++)
			{
				lm[i][j] = 0;
			}
		}
		
		for (i=BIRTH_MAP_R; i<rows-BIRTH_MAP_R; i++)
		{
			for (j=BIRTH_MAP_R; j<cols-BIRTH_MAP_R; j++)
			{
				if((xt[i][j]==0)&&((xt[i+1][j]==1)||
				   (xt[i-1][j]==1)||(xt[i][j+1]==1)||
				   (xt[i][j-1]==1)))
				{
					for (ii=-BIRTH_MAP_R; ii<=BIRTH_MAP_R; ii++)
					{
						m = (int)(sqrt((double)(BIRTH_MAP_R*BIRTH_MAP_R-ii*ii))+0.5);
						for (jj=-m; jj<=m; jj++)
						{
							lm[i+ii][j+jj] = 1;
						}
					}
				}
				else if(xt[i][j]==1)
				{
					lm[i][j] = 1;
				}
			}
		}
		
		// for output image
		for (i=0; i<rows; i++)
		{
			for (j=0; j<cols; j++)
			{
				output_img.mono[i][j] = (int)(lm[i][j] * 255.);
			}
		}

		sprintf(segfileName, "%s_birthmap.tiff",outfilePrefix);

		if ((fp = fopen(segfileName, "wb")) == NULL ) 
		{
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) 
		{
			printf("Error writing TIFF file %s\n", outfileName);
			exit(1);
		}
		fclose(fp);
	}

	//*****************************************************************************
	//		Necking Denting MPP
	//*****************************************************************************
	
	/*Create Blur Matrix*/
	if(1)
	{
		mpp.blur = (double **)get_img(mpp.blur_size, mpp.blur_size, sizeof(double));
		dsum = 0.;
		for (i = 0; i < mpp.blur_size; i++)
		{
			for (j = 0; j < mpp.blur_size; j++)
			{
				di = i-(mpp.blur_size-1)/2.;
				dj = j-(mpp.blur_size-1)/2.;
				mpp.blur[i][j] = exp(-sqrt(di*di+dj*dj)/mpp.sigma);
				dsum += mpp.blur[i][j];
			}
		}
		for (i = 0; i < mpp.blur_size; i++)
		{
			for (j = 0; j < mpp.blur_size; j++)
			{
				mpp.blur[i][j] = mpp.blur[i][j]/dsum;
			}
		}
	}	


	if(mpp.optimization_type == 0)
	{ // RJMCMC
		
		np_num = nd_mpp_rjmcmc(yfiltered, lm, mu, vari, variance, mp, mpp, total_e, mp_num, cols, rows);
	}
	else if(mpp.optimization_type == 1)
	{ // Multiple Birth and Death
		np_num = nd_mpp_multiple_birth_n_death(yfiltered, lm, mu, vari, variance, mp, mpp, total_e, mp_num, cols, rows);
	}
	else
	{// if(mpp.optimization_type == 2){ // RJMCMC Quility Candy model
		np_num = QuilityCandyInterface(yfiltered, lm, mu, vari, variance, mp, mpp, total_e, mp_num, cols, rows);
	}

	
	
	
	printf("\nCalculating New Betas \n\n");

	beta[0] = 0.0;

	calculate_betaimg(beta_dimg, beta, mp, np_num, mpp, cols, rows); // if mpp.gaussian_tau big, beta image channel become narrow
	double_to_uchar(beta_dimg[0], lm_img, cols, rows);
	

	printf("EM/MPM with Adaptive Betas\n");

	enable_blur = 0;
	emmpm_betaimg(yimg, xt, beta_dimg, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);

	
	
	end_time=clock();
	/* Print Statistics */
	if(1)
	{
		end_time=clock();	
		running_time = (double)(end_time-start_time)/CLOCKS_PER_SEC;

		printf("Running time is:%1.0f s\n",running_time);
		avg_width = 0;
		avg_amp = 0;
		avg_error = 0;
		avg_se = 0;
		avg_U = 0;
		for(i=0;i<np_num;i++)
		{
			avg_width	+=mp[i].width;
			avg_amp		+=mp[i].e[0];
			avg_error	+=mp[i].e[1];
			avg_se		+=mp[i].single_E;
			dtmp		=mpp.alpha*mp[i].single_E+(1-mpp.alpha)*mp[i].multiple_E;
			avg_U		+=dtmp;
		}
		avg_width	= avg_width	/np_num;
		avg_amp		= avg_amp	/np_num;
		avg_error	= avg_error	/np_num;
		avg_se		= avg_se	/np_num;
		avg_U		= avg_U		/np_num;

		printf("gaussian tau = %1.4f\n",mpp.gaussian_tau);
		printf("average width = %1.4f\n",avg_width);
		printf("average amp = %1.4f\n",avg_amp);
		printf("average error = %1.4f\n",avg_error);
		printf("average se = %1.4f\n",avg_se	);
		printf("average U = %1.4f\n",avg_U	);
	}

	/*write objects to a text file*/
	if(1)
	{
		sprintf(filename, "%s_all_mp.txt",outfilePrefix);
		if ((fp = fopen(filename, "wb")) == NULL ) {
			printf("Cannot open file %s\n", filename);
			exit(1);
		}
		fprintf(fp,"\r\nObject num , %d\r\n",np_num);
		fprintf(fp,"\r\n");
		for(m = 0; m<np_num; m++){
			fprintf(fp,"idx ,%d\r\n",mp[m].num);
			fprintf(fp,"center position (x,y) ,%1.2f,%1.2f\r\n",mp[m].center.x,mp[m].center.y);
			fprintf(fp,"width ,%1.2f\r\n",mp[m].width);
			fprintf(fp,"length ,%1.2f\r\n",mp[m].length);
			fprintf(fp,"theta ,%1.2f\r\n",mp[m].theta);
			fprintf(fp,"single_E ,%1.2f,\r\n",mp[m].single_E);
			fprintf(fp,"e[0], e[1], e[2],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[0],mp[m].e[1],mp[m].e[2]);
			fprintf(fp,"e[3], e[4], e[5],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[3],mp[m].e[4],mp[m].e[5]);
			fprintf(fp,"e[6], e[7], e[8],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[6],mp[m].e[7],mp[m].e[8]);
			fprintf(fp,"e[9], e[10], e[11],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[9],mp[m].e[10],mp[m].e[11]);
			fprintf(fp,"e[12], e[13], e[14],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[12],mp[m].e[13],mp[m].e[14]);
			fprintf(fp,"e[15], e[16], e[17],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[15],mp[m].e[16],mp[m].e[17]);
			fprintf(fp,"e[18], e[19], e[20],  %d, %1.2f, %1.2f\r\n",(int)mp[m].e[18],mp[m].e[19],mp[m].e[20]);
			fprintf(fp,"e[21], e[22], e[23],  %1.2f, %1.2f, %1.2f\r\n",mp[m].e[21],mp[m].e[22],mp[m].e[23]);
			fprintf(fp,"e[24], e[25], e[26],  %1.2f, %1.2f, %1.2f\r\n",mp[m].e[24],mp[m].e[25],mp[m].e[26]);
			fprintf(fp,"e[27], , ,  %1.2f,,\r\n",mp[m].e[27]);
			fprintf(fp,"\r\n");
		}
		fclose(fp);
	}

	/*Write Detected Denting/Objects in a colored image*/
	if(0)
	{
		draw_all_nds_2D(mp, mpp, np_num, yimg, laplacian, cols, rows);
		//get_TIFF(&output_img, rows, cols, 'g');
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
			{
				output_img.mono[i][j] = laplacian[i][j];
				output_color_img.color[0][i][j] = 255- (int)yimg[i][j];// + laplacian[i][j];
				output_color_img.color[1][i][j] = 255- (int)yimg[i][j];
				output_color_img.color[2][i][j] = 255- (int)yimg[i][j];
			}
			
		if ((fp = fopen("Denting_MPP_colored.tiff", "wb")) == NULL ) {
			printf("Cannot open file %s\n", "Denting_MPP.tiff");
			exit(1);
		}
		if (write_TIFF(fp, &output_color_img)) {
			printf("Error writing TIFF file %s\n", "Denting_MPP.tiff");
			exit(1);
		}
			fclose(fp);
			

		if ((fp = fopen("Denting_MPP_objects.tiff", "wb")) == NULL ) {
			printf("Cannot open file %s\n", "Denting_MPP_objects.tiff");
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) {
			printf("Error writing TIFF file %s\n", "Denting_MPP.tiff");
			exit(1);
		}
		fclose(fp);

	}
		
	

	// PMP
	//misclassed = calculate_PMP(xt, gt, classes, rows, cols);
	//printf("misclassed = %1.4f\n", misclassed);
	
	
	/* Save output Image and Write Parameters to a text file*/
	if(1)
	{
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		sprintf(segfileName, "%s_mpp_segmentation.tiff",outfilePrefix);
	

		if ((fp = fopen(segfileName, "wb")) == NULL ) {
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) {
			printf("Error writing TIFF file %s\n", outfileName);
			exit(1);
		}
		fclose(fp);

		sprintf(outfileName, "%s_mpp_segmentation_parameters.txt",outfilePrefix);

		if ((fp = fopen(outfileName, "wb")) == NULL ) {
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}

		fprintf(fp,"Summary\r\n");
		fprintf(fp,"Image Name = %s\r\n",
				argv[1]);
		
		fprintf(fp,"Iterations = %d\r\nRunning Time = %1.0f minutes\r\n",
				   mpp.iter_num, running_time/60.0);

				
		fprintf(fp,"\r\nChannel Parameters\r\n");
		fprintf(fp,"widthmin = %d\r\nwidthmax = %d\r\nlengthmin = %d\r\nlengthmax = %d\r\n",
				mpp.widthmin, mpp.widthmax, mpp.lengthmin, mpp.lengthmax);
				
		fprintf(fp,"\r\nDenting/Necking Channels Parameters\r\n");
		fprintf(fp,"dent_l_w_ratio = %1.2f\r\nlambda_a = %1.2f\r\nlambda_l = %1.2f\r\n",
				mpp.dent_l_w_ratio, mpp.amp_th, mpp.lambda_l);
				
		fprintf(fp,"\r\nOptimization Parameters\r\n");
		fprintf(fp,"p_birth = %1.2f\r\np_death = %1.2f\r\np_translation = %1.2f\r\np_dilation = %1.2f\r\np_rotation = %1.2f\r\np_switching = %1.2f\r\n",
				mpp.p_birth, mpp.p_death, mpp.p_translation, mpp.p_dilation, mpp.p_rotation, mpp.p_switching);
		
		//fprintf(fp,"single_e_test = %d\r\noptimization_type = %d\r\nfixed_param_num = %d\r\nalm(b_zero) = %1.2f\r\nT0 = %1.2f\r\n",
		//		mpp.test, mpp.optimization_type, mpp.fixed_param_num, mpp.alm, mpp.T0);
		
		fprintf(fp,"\r\nMPP Interaction Parameters\r\n");
		fprintf(fp,"lambda_e = %1.2f\r\nlambda_a = %1.2f\r\nlambda_l = %1.2f\r\nlambda_s = %1.2f\r\nlambda_nc = %1.2f\r\nlambda_dc = %1.2f\r\n",
				mpp.lambda_e, mpp.lambda_a, mpp.lambda_l, mpp.lambda_s, mpp.lambda_nc,mpp.lambda_dc);
				
		fprintf(fp,"betampp = %1.2f\r\nvk(delta) = %1.2f\r\nalpha = %1.2f\r\nde_coeff = %1.8f\r\nnd_type_num = %d\r\n",
				mpp.betampp, mpp.vk, mpp.alpha, mpp.de_coeff, mpp.nd_type_num);

		fprintf(fp,"DELTA_TRANSLATION,=,%1.2f \r\nDELTA_DILATION,=,%1.2f \r\nDELTA_ROTATION,=,%1.2f \r\n",
					DELTA_TRANSLATION,			DELTA_DILATION,			DELTA_ROTATION);

				fprintf(fp,"length_th = %1.2f\r\nsymmetry_th = %1.2f\r\nWEIGHT_AMP = %1.2f\r\nWEIGHT_OFFSET = %1.2f\r\nNECK_DISCON_TH = %1.2f\r\nDENT_DISCON_TH = %1.2f\r\n",
				mpp.length_th, mpp.symmetry_th, WEIGHT_AMP, WEIGHT_OFFSET, NECK_DISCON_TH, DENT_DISCON_TH);
				

		fprintf(fp,"\r\nEM/PMP Parameters\r\n");
		fprintf(fp,"\r\nbeta0,=,%1.2f \r\nbeta1,=,%1.2f \r\nsigma,=,%1.2f \r\nenable_blur,=,%d \r\npmp,=,%1.2f \r\n \r\n ",
				beta[0], beta[1], sigma, enable_blur, misclassed);

			//	fprintf(fp,"parameters\r\nImage Name = %s \r\n Hard_repulsion = %d     Gaussian_tau = %1.2f\r\n",
			//	argv[1], mpp.hard_repulsion, mpp.gaussian_tau);
			//fprintf(fp,"error_th = %1.2f\r\namp_th = %1.2f\r\niter_num = %d\r\nrunning time = %1.0f minutes\r\n",
			//	mpp.error_th, mpp.amp_th, mpp.iter_num, running_time/60.0);
			

		fprintf(fp,"Iteration information, ");
		#if 0
			for (i=0;i<mpp.iter_num/RJMCMC_ITER_DIV;i++)
				fprintf(fp,"%d, ", i*RJMCMC_ITER_DIV);
			fprintf(fp,"\r\n");
			fprintf(fp,"energy, ");
			for (i=0;i<mpp.iter_num/RJMCMC_ITER_DIV;i++)
				fprintf(fp,"%f, ", total_e[i]);
			fprintf(fp,"\r\n");
			fprintf(fp,"object number, ");
			for (i=0;i<mpp.iter_num/RJMCMC_ITER_DIV;i++)
				fprintf(fp,"%d, ", mp_num[i]);
			fprintf(fp,"\r\n");
		#else
		//	for (i=1;i<=mpp.iter_num;i++)
		//		fprintf(fp,"%d, ", i);
		//	fprintf(fp,"\r\n");
		//	fprintf(fp,"energy, ");
		//	for (i=1;i<=mpp.iter_num;i++)
		//		fprintf(fp,"%f, ", total_e[i]);
		//	fprintf(fp,"\r\n");
		//	fprintf(fp,"object number, ");
		//	for (i=1;i<=mpp.iter_num;i++)
		//		fprintf(fp,"%d, ", mp_num[i]);
		//	fprintf(fp,"\r\n");
		#endif
		fclose(fp);
	}

//	for(i=0; i < input_img.height; i++)
//	{
//		for(j=0; j< input_img.width; j++)
//		{
//			gt[i][j] = (int)(input_gt_img.mono[i][j];
//		}
//	}
	printf("Finding Difference Image \n");
	difference_image(xt, gt ,input_img.height, input_img.width,outfilePrefix);

	free_TIFF(&output_color_img);
	free_TIFF(&output_img);
	free_TIFF(&input_img);
	free_TIFF(&input_gt_img);
	free_img((void **)yimg);
	free_img((void **)yfiltered);
	free_img((void **)laplacian);
	free_img((void **)gt);
	free_img((void **)xt);
	free_img((void **)blur);
	free_img((void **)mpp.blur);
	free_img((void **)lm);
	free_img((void **)beta_dimg[0]);
	free_img((void **)beta_dimg[1]);
	free_img((void **)lm_img);

	

	printf(">>>>>>>>>>>>>>>>End of Program<<<<<<<<<<<<<\n");
	return 0;
}

void difference_image(unsigned char **img1, unsigned char  **img2,int height, int width, char* name)
{
	int i=0,j=0;
	struct TIFF_img output_tiff;
	FILE *fp=NULL;
	char aa[40];
	
	get_TIFF ( &output_tiff, height, width, 'c' );
	
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if(img1[i][j] == img2[i][j])
			{
				output_tiff.color[0][i][j] = img1[i][j] * 255;
				output_tiff.color[1][i][j] = img1[i][j] * 255;
				output_tiff.color[2][i][j] = img1[i][j] * 255;
			}
			else if(img1[i][j] < img2[i][j])
			{
				//Make Green Where img2 detected object but img1 didnt 
				output_tiff.color[0][i][j] = 0;
				output_tiff.color[1][i][j] = 255;
				output_tiff.color[2][i][j] = 0;
			}
			else
			{
				//Make Red Where img1 detected object but img2 didnt 
				output_tiff.color[0][i][j] = 255;
				output_tiff.color[1][i][j] = 0;
				output_tiff.color[2][i][j] = 0;
			}

		}
	}
	
	strcat(name,"_diff.tiff");
	/* open image file */
	if ( ( fp = fopen ( name, "wb" ) ) == NULL ) {
	    fprintf ( stderr, "cannot open file tif file\n");
	    exit ( 1 );
	  }

	  /* write image */
	  if ( write_TIFF ( fp, &output_tiff ) ) {
	    fprintf ( stderr, "error writing TIFF file" );
	    exit ( 1 );
	  }

	  free_TIFF(&output_tiff);
	  /* close image file */
	  fclose ( fp );
}

