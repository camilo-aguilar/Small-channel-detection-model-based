#include <stdio.h>
#define _USE_MATH_DEFINES	
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

#include "allocate.h"
#include "random.h"
#include "tiff.h"
#include "ndmulti_e.h"
#include "em.h"
#include "QualityCandy.h"
#include "Parameters.h"



#if INCLUDE_OPENCV
	#include "gui_functions.h"
#endif

#if 0

				for(i=0; i<rows;i++)
					for(j=0;j<cols;j++)
								yimg_int[i][j] = yimg[i][j];

				double my_theta;
				for(my_theta= THETA_MIN; my_theta <= THETA_MAX; my_theta = my_theta + (THETA_MAX-THETA_MIN)/20)
				{
					double super_min = 100000;
					double super_max = -100000;
					for(i=0; i<rows;i++)
						for(j=0; j<cols; j++)
							{
								channel_img[i][j] = im_filter(yimg_int, rows, cols, my_theta, i, j);
								if(channel_img[i][j] < super_min)
									super_min = channel_img[i][j];
								if(channel_img[i][j] > super_max)
									super_max = channel_img[i][j];
							}

					printf("Super Min: %f \n", super_min);
					printf("Super Max: %f \n", super_max);

					for(i=0; i<rows;i++)
						for(j=0; j<cols; j++)
						xt[i][j] = 255 -(((channel_img[i][j]-super_min)*255.00)/(super_max-super_min));
			 

				//#endif
				sprintf(segfileName, "%s_derivative_at_%.2f.tiff",outfilePrefix,my_theta*180/_PI);
				/* TO REMOVE!!! =================================================================================================*/
				

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
			}
		}
		return 0;

#endif



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

int QuilityCandyInterface(unsigned char **yimg, double **channel_img, double **lm, double *mean, double *vari, 
			double variance, NeckDent **mp, MPP_Parameters mpp, double *total_e, int *mp_num, 
			int cols, int rows)
{
    struct TIFF_img input_img;
  double **img1,****img2,****img3,**patch;
  int **img_seg;
    int32_t i,j;//,pixel,pixel2;

  input_img.height = rows;
  input_img.width = cols;
  int patch_len = 100;
  
  /* Allocate image of double precision floats */
  img1 = (double **)get_img(input_img.width,input_img.height,sizeof(double));
  //img_out = (double **)get_img(input_img.width,input_img.height,sizeof(double));
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
		channel_img[i][j] = 0;
		img_seg[i][j]=0;
	 }
  }

  srandom2(1);
  int iter = mpp.iter_num; //1000000;//100000000;
  double T = mpp.T0;//KDW 1/10.0;


   double *mu = (double *)malloc(MAX_K*sizeof(double));
   double *cov = (double *)malloc(MAX_K*sizeof(double));
   //K_means(img1,mu,cov,input_img.height,input_img.width,3);



  (*mp_num) = Candy_Model(img1,lm,img2,img3,channel_img,img_seg,input_img.height,input_img.width,
		T,iter,patch,patch_len,mu,cov,3,mp,mpp);

	
  free_img( (void**)img1 );
  free(img2);
  free(img3);
  //free_img( (void**)img_out ); 
  free_img((void**)img_seg);
  free_img((void**)patch);
  free((void*)mu);
  free((void*)cov);

  return  (*mp_num);
}

/* Finds the difference between img1 and img2. It will color green if img1 > img2, white if img1 = img2, red if img2 > img1 */
void difference_image(unsigned char **img1, unsigned char  **img2,int height, int width, char* name);

/* Saves image + detected channels */
void save_channel_image(double **channel_img, unsigned char  **img,int height, int width, char* name);

/*
MAIN
*/
int main( int argc , char** argv)
{
	printf(">>>>>>>>>>>>>>>>Start of the Program<<<<<<<<<<<<<\n");
	
	if(argc < 2)
	{
		printf("Error, only 2 arguemnts\n");
		exit(1);
	} 	
	/* Image Structures */
	struct TIFF_img input_img, input_gt_img, output_img, output_color_img;
	unsigned char **yimg = NULL, **yfiltered, **laplacian, **lm_img;
	double **lm, **beta_dimg[2], **channel_img;
	
	/*File names* arrays*/
	char infileName[1024], outfileName[1024];
	char segfileName[1024], outfilePrefix[1024];
	char gtfileName[1024], filename[1024];
	sprintf(infileName, "%s.tiff",argv[1]);
	sprintf(segfileName, "%s_seg.tiff",argv[1]);
	sprintf(gtfileName, "%s_gt.tiff",argv[1]);
	sprintf(outfilePrefix, "%s",argv[1]);

	/* Helper Variables */
	FILE *fp;
	int i, j, ii, jj, rows, cols;
	int m;
	double running_time;
	clock_t start_time, end_time;
	int np_num;
	
	// for EM/MPM
	double mu[CLASSES], mean, vari[CLASSES], variance, min, max;
	int cnt[CLASSES];
	double beta[MAX_CLASSES], gamma[MAX_CLASSES];
	int classes = 2;
	unsigned char **xt, **gt;  /* output : entropy image */
	double **blur, sigma = 0., dsum, sum[CLASSES], di, dj, misclassed;
	int blur_size = 5, enable_blur = 0;
	int run_emmpm = 0;

	
	//Variables for MPP
	double total_e[MAX_MPP_ITER_NUM+1];
	int	   mp_num[MAX_MPP_ITER_NUM+1];
	double avg_width, avg_amp, avg_error, avg_se, avg_U, dtmp;
	beta[0] = BETTA_FOR_MPP;
	beta[1] = BETTA_FOR_MPP;
	gamma[0] = 0.;
	gamma[1] = 0.;
	sigma = 0.58;
	NeckDent *mp = NULL;

	readseed();
		
	
	/* Parameter Parsing */
	MPP_Parameters mpp = parse_input_parameters(argc,argv);
	
	_print_current_parameters(mpp);

	start_time=clock();
	/*Read Image */	
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
					#if IMAGE_IVERTED
						yimg[i][j] = 255 - input_img.mono[i][j];
					#else
						yimg[i][j] = input_img.mono[i][j];
					#endif
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
		channel_img = (double **)get_img(input_img.width,input_img.height,sizeof(double));
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

	/* Calculate Blur Matrix. If necesary */
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

			start_time=clock();
			printf("start time is:%1.0f ms\n",(double)(start_time)/CLOCKS_PER_SEC*100);
			
			enable_blur = 0;
			//#if PRE_SEG_EM_MPM
			//	emmpm(yimg, xt, beta, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);
			//#else
			for(i=0; i<rows;i++)
				for(j=0; j<cols; j++)
					if(yimg[i][j]> 150)
						xt[i][j] = 1;
					else
						xt[i][j] = 0;
			


				sprintf(segfileName, "%s_seg.tiff",outfilePrefix);
				/* TO REMOVE!!! =================================================================================================*/
				

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
		double test_th = fabs((mu[0]-mu[1]))/sqrt(vari[0]/cnt[0] + vari[1]/cnt[1]);
		printf("Recommended error_th: %f\n", test_th);
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
				   (xt[i][j-1]==1) ))
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
	{ // if(mpp.optimization_type == 2){ // RJMCMC Quility Candy model
		np_num = QuilityCandyInterface(yfiltered, channel_img, lm, mu, vari, variance, &mp, mpp, total_e, mp_num, cols, rows);
	}

	
	
	//printf("\nCalculating New Betas \n\n");
	beta[0] = 0.0;

	//calculate_betaimg(beta_dimg, beta, mp, np_num, mpp, cols, rows); // if mpp.gaussian_tau big, beta image channel become narrow

	
	//double_to_uchar(beta_dimg[0], lm_img, cols, rows);
	

	//printf("EM/MPM with Adaptive Betas\n");

	//enable_blur = 0;
	//emmpm_betaimg(yimg, xt, beta_dimg, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);

	
	
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


	}

	/*write objects to a text file*/
	if(0)
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
		
	
	
	


	//	for(i=0; i < input_img.height; i++)
	//	{
	//		for(j=0; j< input_img.width; j++)
	//		{
	//			gt[i][j] = (int)(input_gt_img.mono[i][j];
	//		}
	//	}
	//printf("Finding Difference Image \n");
	//difference_image(xt, gt ,input_img.height, input_img.width,outfilePrefix);
		
	printf("Saving Channel Image\n");


		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				yimg[i][j] = input_img.mono[i][j];
			}
		}
		


		char parameter_values[100];
		float my_width = W_MAX;
		float my_length = L_MAX;
		sprintf(parameter_values, "gamma_%.2f_wf_%.2f_ws_%.2f_wd_%.2f_weo_%.2f_wio_%.2f_wid_%.2f_len_%.2f_error_th_%.2f_iter_%d_",mpp.gamma_d, mpp.w_f , mpp.w_s, mpp.w_d, mpp.w_eo, mpp.w_io,my_width, my_length, mpp.error_th,mpp.iter_num);


		//save_channel_image(channel_img, yimg,input_img.height, input_img.width, parameter_values);
		char output_name_temp[30] = "output";
		save_channel_image(channel_img, yimg,input_img.height, input_img.width, output_name_temp);

		free(mp);
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
		free_img((void**)channel_img);


    	
		printf(">>>>>>>>>>>>>>>>End of Program<<<<<<<<<<<<<\n");
		return 0;
}



void difference_image(unsigned char **img1, unsigned char  **img2,int height, int width, char* name)
{
	int i=0,j=0;
	struct TIFF_img output_tiff;
	FILE *fp=NULL;
	char aa[40];
	strcpy(aa,name);
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
	
	strcat(aa,"_diff.tiff");
	/* open image file */
	if ( ( fp = fopen ( aa, "wb" ) ) == NULL ) {
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

	/* Very Important, it saves the image INVERTED*/
void save_channel_image(double **channel_img, unsigned char  **img,int height, int width, char* name)
{

	int i=0,j=0;
	struct TIFF_img output_tiff;
	FILE *fp=NULL;
	char aa[200];
	strcpy(aa,name);
	get_TIFF ( &output_tiff, height, width, 'c' );
	
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if(channel_img[i][j]  == 0)
			{
				#if IMAGE_IVERTED
					output_tiff.color[0][i][j] = img[i][j];
					output_tiff.color[1][i][j] = img[i][j];
					output_tiff.color[2][i][j] = img[i][j];
				#else 
					output_tiff.color[0][i][j] = 255 - img[i][j];
					output_tiff.color[1][i][j] = 255 - img[i][j];
					output_tiff.color[2][i][j] = 255 - img[i][j];
				#endif
			}
			else
			{
				//Make Green Where img2 detected object but img1 didnt 
				output_tiff.color[0][i][j] = 255;
				output_tiff.color[1][i][j] = 0;
				output_tiff.color[2][i][j] = 0;
			}


		}
	}
	
	strcat(aa,"_results.tiff");
	/* open image file */
	if ( ( fp = fopen ( aa, "wb" ) ) == NULL ) {
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


