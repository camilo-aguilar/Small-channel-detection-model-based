/* mpm.c */
#include "mpm.h"
#include "em.h"

double gaussianRandom(double average, double stdev) {
  double v1, v2, s, temp;

  do {
    v1 =  2 * ((double) rand() / RAND_MAX) - 1;      // -1.0 ~ 1.0
    v2 =  2 * ((double) rand() / RAND_MAX) - 1;      // -1.0 ~ 1.0
    s = v1 * v1 + v2 * v2;
  } while (s >= 1 || s == 0);

  s = sqrt( (-2 * log(s)) / s );

  temp = v1 * s;
  temp = (stdev * temp) + average;

  return temp;
}

#ifdef	DEF_COVARIANCE
void blurring_term(unsigned char **yimg, double ***bt, double ***btkv, unsigned char **xt, double *m, double *v,
				double ****kmat, double *vmat, int rows, int cols, int classes, double **blur,
				int blur_size, int k)
#else
void blurring_term(unsigned char **yimg, double ***bt, double ***btkv, unsigned char **xt, double *m, double *v,
				int rows, int cols, int classes, double **blur,	int blur_size, int k)
#endif
{
	int i, j, ii, jj, xx, yy, l, gaussian = DEF_GAUSSIAN;
	double **dimg, sum, y_tilda_qt;
	int d = (blur_size-1)/2;

	dimg = (double **)get_img(cols, rows, sizeof(double));

	/* ys */
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++){
			if(gaussian)
				dimg[i][j] = gaussianRandom(m[xt[i][j]], sqrt(v[xt[i][j]]));
			else
				dimg[i][j] = m[xt[i][j]];
		}
	/* yt */
	for (l=0; l<classes; l++){
		if(gaussian)
			y_tilda_qt = gaussianRandom(m[l], sqrt(v[l]));
		else
			y_tilda_qt = m[l];
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++){
				sum = 0.;
				for(yy=-d; yy<=d; yy++)
					for(xx=-d; xx<=d; xx++){
						ii = i+yy;
						if (ii<0) ii = 0;
						else if (ii>=rows) ii = rows-1;
						jj = j+xx;
						if (jj<0) jj = 0;
						else if (jj>=cols) jj = cols-1;
						if((xx == 0)&&(yy == 0))
							sum += blur[yy+d][xx+d]*y_tilda_qt;
						else
							sum += blur[yy+d][xx+d]*dimg[ii][jj];
					}
#ifndef	DEF_COVARIANCE
				btkv[l][i][j] = ((double)yimg[i][j]-sum)*((double)yimg[i][j]-sum);
				bt[l][i][j] = btkv[l][i][j]/(2.0*v[l]);
#else
				for(x=0; x<=num_vpar; x++)
					km[x] = kmat[k-1][x][ii][jj];

				if (xt[i][j]!=l){
					km[xt[i][j]] -= blur[d][d]*blur[d][d];
					km[l] += blur[d][d]*blur[d][d];
					if(j+1>=cols){
						km[xt[i][j]*classes+xt[i][cols]+classes] -= blur[d][d]*blur[d][d+1];
						km[l*classes+xt[i][cols]+classes] += blur[d][d]*blur[d][d+1];
					}
					else{
						km[xt[i][j]*classes+xt[i][j+1]+classes] -= blur[d][d]*blur[d][d+1];
						km[l*classes+xt[i][j+1]+classes] += blur[d][d]*blur[d][d+1];
					}
					if(j-1<0){
						km[xt[i][j]*classes+xt[i][0]+classes] -= blur[d][d]*blur[d][d-1];
						km[l*classes+xt[i][0]+classes] += blur[d][d]*blur[d][d-1];
					}
					else{
						km[xt[i][j]*classes+xt[i][j-1]+classes] -= blur[d][d]*blur[d][d-1];
						km[l*classes+xt[i][j-1]+classes] += blur[d][d]*blur[d][d-1];
					}
					if(i+1>=rows){
						km[xt[i][j]*classes+xt[rows][j]+classes*classes+classes] -= blur[d][d]*blur[d+1][d];
						km[l*classes+xt[rows][j]+classes*classes+classes] += blur[d][d]*blur[d+1][d];
					}
					else{
						km[xt[i][j]*classes+xt[i+1][j]+classes*classes+classes] -= blur[d][d]*blur[d+1][d];
						km[l*classes+xt[i+1][j]+classes*classes+classes] += blur[d][d]*blur[d+1][d];
					}
					if(i-1<0){
						km[xt[i][j]*classes+xt[0][j]+classes*classes+classes] -= blur[d][d]*blur[d-1][d];
						km[l*classes+xt[0][j]+classes*classes+classes] += blur[d][d]*blur[d-1][d];
					}
					else{
						km[xt[i][j]*classes+xt[i-1][j]+classes*classes+classes] -= blur[d][d]*blur[d-1][d];
						km[l*classes+xt[i-1][j]+classes*classes+classes] += blur[d][d]*blur[d-1][d];
					}

				}
				tmp = 0;
				for(x=0; x<=num_vpar; x++)
					tmp += km[x]*vmat[x];

				btkv[l][i][j] = ((double)yimg[i][j]-sum)*((double)yimg[i][j]-sum);
				bt[l][i][j] = btkv[l][i][j]/(2.0*tmp);
#endif // DEF_COVARIANCE
			}
	}
	free_img((void **)dimg);
	return;
}


#ifdef	DEF_COVARIANCE
void mpm(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double *beta, double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double delta, double **diff, double *s_p, double ****kmat,
		double *vmat)
#else
void mpm(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double *beta, double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double **diff, double *s_p)
#endif
{
	//printf("Beta, %f, %f %f \n", beta[0],beta[1],beta[2]);
	double **yk[MAX_CLASSES], sqrt2pi, current, con[MAX_CLASSES], denom[MAX_CLASSES];
	double x, post[MAX_CLASSES], sum, coeff[MAX_CLASSES];
	int i, j, ii, jj, k, l, mm, prior[MAX_CLASSES], pr_sum[MAX_CLASSES], d = (blur_size-1)/2;
	double c_sum, numer_sum[MAX_CLASSES];
	double  **bt[MAX_CLASSES], **btkv[MAX_CLASSES];

#ifdef MPM_IMG_OUTPUT
	struct TIFF_img output_img;
	char  outfileName[1024];
	FILE* fp;

	get_TIFF(&output_img, rows, cols, 'g');
#endif

	for (l = 0; l < classes; l++){
		bt[l] = (double **)get_img(cols, rows, sizeof(double));
		btkv[l] = (double **)get_img(cols, rows, sizeof(double));
	}

	for(k=0;k<classes;k++){ // diff
		numer_sum[k] = 0.;
		pr_sum[k] = 0;
		for(l=0;l<classes+1;l++){ // coeff
			diff[k][l] = 0.;
		}
	}

	sqrt2pi = sqrt(2.0 * M_PI);

	for (l = 0; l < classes; l++) {
		con[l] = -log(sqrt2pi * sqrt(v[l]));
		denom[l] = -2.0 * v[l];
	}

//	printf("kmat[%d][%d][%d][%d] = %1.2f\n",n-1,num_vpar-1,rows-1,cols-1,kmat[n-1][num_vpar-1][rows-1][cols-1]);

	/* Allocate space for yk[][][] */
	for (l = 0; l < classes; l++)
		yk[l] = (double **)get_img(cols, rows, sizeof(double));

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++) {
			mm = y[i][j];
			for (l = 0; l < classes; l++) {
				pr[l][i][j] = 0;  // reset content of (16)
				yk[l][i][j] = con[l] + ((mm - m[l]) * (mm - m[l]) / denom[l]);
			}
		}

	for (k = 0; k < n; k++){
//		printf("mpmiter = %d\n",k);
		//	printf("generate yt... ");
#ifdef	DEF_COVARIANCE
			if(k==0)
				blurring_term(y, bt, btkv, xt, m, v, kmat, vmat, rows, cols, classes, blur, blur_size, n);
			else
				blurring_term(y, bt, btkv, xt, m, v, kmat, vmat, rows, cols, classes, blur, blur_size, k);
#else
			if(k==0)
				blurring_term(y, bt, btkv, xt, m, v, rows, cols, classes, blur, blur_size, n);
			else
				blurring_term(y, bt, btkv, xt, m, v, rows, cols, classes, blur, blur_size, k);
#endif
			//	printf("end\n");
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++) {
				sum = 0;
				for (l = 0; l < classes; l++) {
					bt[l][i][j] = con[l] - bt[l][i][j];
					prior[l] = 0;
					if (i - 1 >= 0)
						if (xt[i - 1][j] != l)
							(prior[l])++;
					if (i + 1 <= rows - 1)
						if (xt[i + 1][j] != l)
							(prior[l])++;
					if (j - 1 >= 0)
						if (xt[i][j - 1] != l)
							(prior[l])++;
					if (j + 1 <= cols - 1)
						if (xt[i][j + 1] != l)
							(prior[l])++;
					// New Clique
//#define _8NN
#ifdef _8NN
					if ((i - 1 >= 0)&&(j - 1 >= 0))
						if (xt[i - 1][j - 1] != l)
							(prior[l])++;
					if ((i - 1 >= 0)&&(j + 1 <= cols - 1))
						if (xt[i - 1][j + 1] != l)
							(prior[l])++;
					if ((i + 1 <= rows - 1)&&(j - 1 >= 0))
						if (xt[i + 1][j - 1] != l)
							(prior[l])++;
					if ((i + 1 <= rows - 1)&&(j + 1 <= cols - 1))
						if (xt[i + 1][j + 1] != l)
							(prior[l])++;
#endif
					if(enable_blur){
						post[l] = exp(bt[l][i][j] - beta[l] * (double)(prior[l]) - gamma[l]);
					}
					else
						post[l] = exp(yk[l][i][j] - beta[l] * (double)(prior[l]) - gamma[l]);
					sum += post[l];
				}
				x = random2();
				current = 0;

				for (l = 0; l < classes; l++) {
					if ((x >= current) && (x <= (current + post[l] / sum))) {
						xt[i][j] = l;
						pr[l][i][j] += 1.0;
					}
					current += post[l] / sum;
				}

				// to calculate blurring term EM update KDW20120222
				//
				// -------- Mean --------------
				// diff[k][l] = df(one mpm iter,one point)/dM[m]
				// = (y[i][j]-(coeff[0]*M[0]+coeff[1]*M[1]+...+coeff[classes-1]*M[classes-1])
				// *(blur[0][0]*(xt[i+d][j+d]=m)+...+blur[2d][2d]*(xt[i-d][j-d]=m))
				//
				// diff[k][l] = df/dM[m] = sum(all mpm iter) of sum(all points) of df(one mpm iter,one point)/dM[m]
				//
				// coeff[m] = (blur[0][0]*(xt[i+d][j+d]=m)+...+blur[2d][2d]*xt[i-d][j-d]=m))
				//
				// -------- Variance -----------
				// numer = (y[i][j]-(coeff[0]*M[0]+coeff[1]*M[1]+...+coeff[classes-1]*M[classes-1])^2
				// *(xt[i-d][j-d]=m)
				//
				// numer_sum[k] = sum(all point) sum (amm mpm iter) of numer
				//
				if(enable_blur){
					if((i>=d)&&(i<rows-d)&&(j>=d)&&(j<cols-d)){
						for(ii=0;ii<classes;ii++)
							coeff[ii] = 0.;
						c_sum = 0.;
						for(ii=-d;ii<=d;ii++)
							for(jj=-d;jj<=d;jj++)
								coeff[xt[i-ii][j-jj]] += blur[ii+d][jj+d];
						for(ii=0;ii<classes;ii++){ // diff
							for(l=0;l<classes;l++) // coeff
								diff[ii][l] += coeff[ii]*coeff[l]/v[xt[i][j]];
							diff[ii][classes] += coeff[ii]*(double)y[i][j]/v[xt[i][j]];
							c_sum += coeff[ii]*m[ii];
						}
						c_sum = (double)y[i][j] - c_sum;
						numer_sum[xt[i][j]] += c_sum*c_sum;
						pr_sum[xt[i][j]]++;
					}
				}

#ifdef DEF_COVARIANCE
				// to calculate blurring term EM update using covariance KDW20120522
				// -------- covariance --------------

				btk[k][i][j] = btkv[xt[i][j]][i][j];

				// Initialize kmat[][][][]
				for(l = 0; l < num_vpar; l++)
					kmat[k][l][i][j] = 0.;

				if((i>d)&&(i<rows-d-1)&&(j>d)&&(j<cols-d-1)){
					for(ii=0;ii<classes;ii++)
						for(jj=0;jj<classes;jj++){
							covh[ii][jj] = 0.;
							covv[ii][jj] = 0.;
						}
					for(ii=-d;ii<=d;ii++)
						for(jj=-d;jj<=d;jj++){
							kmat[k][xt[i-ii][j-jj]][i][j] += blur[ii+d][jj+d]*blur[ii+d][jj+d];
							covh[xt[i-ii][j-jj]][xt[i-ii-1][j-jj]] += blur[ii+d][jj+d]*blur[ii+d+1][jj+d];
							covh[xt[i-ii][j-jj]][xt[i-ii+1][j-jj]] += blur[ii+d][jj+d]*blur[ii+d-1][jj+d];
							covv[xt[i-ii][j-jj]][xt[i-ii][j-jj-1]] += blur[ii+d][jj+d]*blur[ii+d][jj+d+1];
							covv[xt[i-ii][j-jj]][xt[i-ii][j-jj+1]] += blur[ii+d][jj+d]*blur[ii+d][jj+d-1];
						}
					idx = classes;
					for(ii=0;ii<classes;ii++)
						for(jj=0;jj<classes;jj++){
							kmat[k][idx][i][j] = covh[ii][jj];
							kmat[k][idx+classes*classes][i][j] = covv[ii][jj];
							idx++;
						}
				}
#endif
			}

#ifdef MPM_IMG_OUTPUT
		// writing images
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		sprintf(outfileName, "emiter_%2d_mpmiter_%2d.tiff", emiter, k);
		if ((fp = fopen(outfileName, "wb")) == NULL ) {
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) {
			printf("Error writing TIFF file %s\n", outfileName);
			exit(1);
		}
	}
	fclose(fp);
	free_TIFF(&output_img);
#else
	}
#endif
	// to calculate blurring term EM update KDW20120222
	if(enable_blur){
		for(l=0;l<classes;l++)
			s_p[l] = numer_sum[l]/(double)pr_sum[l];
	}

	/* Normalize probabilities */
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			for (l = 0; l < classes; l++)
				pr[l][i][j] = pr[l][i][j] / (double)n;

	/* Clean Up */
	for (l = 0; l < classes; l++){
		free_img((void **)bt[l]);
		free_img((void **)btkv[l]);
		free_img((void **)yk[l]);
	}
}

#ifdef	DEF_COVARIANCE
void mpm_betaimg(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double **beta[], double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double delta, double **diff, double *s_p, double ****kmat,
		double *vmat)
#else
void mpm_betaimg(unsigned char **y, unsigned char **xt, double **btk[], double **pr[], double **beta[], double *gamma,
		double *m, double *v, int rows, int cols, int n, int classes, int emiter, double **blur,
		int blur_size, int enable_blur, double **diff, double *s_p)
#endif
{
	double **yk[MAX_CLASSES], sqrt2pi, current, con[MAX_CLASSES], denom[MAX_CLASSES];
	double x, post[MAX_CLASSES], sum, coeff[MAX_CLASSES];
	int i, j, ii, jj, k, l, mm, prior[MAX_CLASSES], prior8[MAX_CLASSES], pr_sum[MAX_CLASSES], d = (blur_size-1)/2;
	double c_sum, numer_sum[MAX_CLASSES];
	double  **bt[MAX_CLASSES], **btkv[MAX_CLASSES];

#ifdef MPM_IMG_OUTPUT
	struct TIFF_img output_img;
	char  outfileName[1024];
	FILE* fp;

	get_TIFF(&output_img, rows, cols, 'g');
#endif

	for (l = 0; l < classes; l++){
		bt[l] = (double **)get_img(cols, rows, sizeof(double));
		btkv[l] = (double **)get_img(cols, rows, sizeof(double));
	}

	for(k=0;k<classes;k++){ // diff
		numer_sum[k] = 0.;
		pr_sum[k] = 0;
		for(l=0;l<classes+1;l++){ // coeff
			diff[k][l] = 0.;
		}
	}

	sqrt2pi = sqrt(2.0 * M_PI);

	for (l = 0; l < classes; l++) {
		con[l] = -log(sqrt2pi * sqrt(v[l]));
		denom[l] = -2.0 * v[l];
	}

//	printf("kmat[%d][%d][%d][%d] = %1.2f\n",n-1,num_vpar-1,rows-1,cols-1,kmat[n-1][num_vpar-1][rows-1][cols-1]);

	/* Allocate space for yk[][][] */
	for (l = 0; l < classes; l++)
		yk[l] = (double **)get_img(cols, rows, sizeof(double));

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++) {
			mm = y[i][j];
			for (l = 0; l < classes; l++) {
				pr[l][i][j] = 0;  // reset content of (16)
				yk[l][i][j] = con[l] + ((mm - m[l]) * (mm - m[l]) / denom[l]);
			}
		}

	for (k = 0; k < n; k++){
//		printf("mpmiter = %d\n",k);
		//	printf("generate yt... ");
#ifdef	DEF_COVARIANCE
			if(k==0)
				blurring_term(y, bt, btkv, xt, m, v, kmat, vmat, rows, cols, classes, blur, blur_size, n);
			else
				blurring_term(y, bt, btkv, xt, m, v, kmat, vmat, rows, cols, classes, blur, blur_size, k);
#else
			if(k==0)
				blurring_term(y, bt, btkv, xt, m, v, rows, cols, classes, blur, blur_size, n);
			else
				blurring_term(y, bt, btkv, xt, m, v, rows, cols, classes, blur, blur_size, k);
#endif
			//	printf("end\n");
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++) {
				sum = 0;
				for (l = 0; l < classes; l++) {
					bt[l][i][j] = con[l] - bt[l][i][j];
					prior[l] = 0;
					prior8[l] = 0;
					if (i - 1 >= 0)
						if (xt[i - 1][j] != l)
							(prior[l])++;
					if (i + 1 <= rows - 1)
						if (xt[i + 1][j] != l)
							(prior[l])++;
					if (j - 1 >= 0)
						if (xt[i][j - 1] != l)
							(prior[l])++;
					if (j + 1 <= cols - 1)
						if (xt[i][j + 1] != l)
							(prior[l])++;
//#define  EIGHT_NEIGHBOR
#ifdef EIGHT_NEIGHBOR
					// New Clique
					if ((i - 1 >= 0)&&(j - 1 >= 0))
						if (xt[i - 1][j - 1] != l)
							(prior8[l])++;
					if ((i - 1 >= 0)&&(j + 1 <= cols - 1))
						if (xt[i - 1][j + 1] != l)
							(prior8[l])++;
					if ((i + 1 <= rows - 1)&&(j - 1 >= 0))
						if (xt[i + 1][j - 1] != l)
							(prior8[l])++;
					if ((i + 1 <= rows - 1)&&(j + 1 <= cols - 1))
						if (xt[i + 1][j + 1] != l)
							(prior8[l])++;
					if(enable_blur){
						post[l] = exp(bt[l][i][j] - beta[l][i][j] * ((double)prior[l]+0.7071*(double)prior8[l]) - gamma[l]);
					}
					else
						post[l] = exp(yk[l][i][j] - beta[l][i][j] * ((double)prior[l]+0.7071*(double)prior8[l]) - gamma[l]);
#else
					if(enable_blur){
						post[l] = exp(bt[l][i][j] - beta[l][i][j] * (double)(prior[l]) - gamma[l]);
					}
					else
						post[l] = exp(yk[l][i][j] - beta[l][i][j] * (double)(prior[l]) - gamma[l]);
#endif
					sum += post[l];
				}
				x = random2();
				current = 0;

				for (l = 0; l < classes; l++) {
					if ((x >= current) && (x <= (current + post[l] / sum))) {
						xt[i][j] = l;
						pr[l][i][j] += 1.0;
					}
					current += post[l] / sum;
				}

				// to calculate blurring term EM update KDW20120222
				//
				// -------- Mean --------------
				// diff[k][l] = df(one mpm iter,one point)/dM[m]
				// = (y[i][j]-(coeff[0]*M[0]+coeff[1]*M[1]+...+coeff[classes-1]*M[classes-1])
				// *(blur[0][0]*(xt[i+d][j+d]=m)+...+blur[2d][2d]*(xt[i-d][j-d]=m))
				//
				// diff[k][l] = df/dM[m] = sum(all mpm iter) of sum(all points) of df(one mpm iter,one point)/dM[m]
				//
				// coeff[m] = (blur[0][0]*(xt[i+d][j+d]=m)+...+blur[2d][2d]*xt[i-d][j-d]=m))
				//
				// -------- Variance -----------
				// numer = (y[i][j]-(coeff[0]*M[0]+coeff[1]*M[1]+...+coeff[classes-1]*M[classes-1])^2
				// *(xt[i-d][j-d]=m)
				//
				// numer_sum[k] = sum(all point) sum (amm mpm iter) of numer
				//
				if(enable_blur){
					if((i>=d)&&(i<rows-d)&&(j>=d)&&(j<cols-d)){
						for(ii=0;ii<classes;ii++)
							coeff[ii] = 0.;
						c_sum = 0.;
						for(ii=-d;ii<=d;ii++)
							for(jj=-d;jj<=d;jj++)
								coeff[xt[i-ii][j-jj]] += blur[ii+d][jj+d];
						for(ii=0;ii<classes;ii++){ // diff
							for(l=0;l<classes;l++) // coeff
								diff[ii][l] += coeff[ii]*coeff[l]/v[xt[i][j]];
							diff[ii][classes] += coeff[ii]*(double)y[i][j]/v[xt[i][j]];
							c_sum += coeff[ii]*m[ii];
						}
						c_sum = (double)y[i][j] - c_sum;
						numer_sum[xt[i][j]] += c_sum*c_sum;
						pr_sum[xt[i][j]]++;
					}
				}

#ifdef DEF_COVARIANCE
				// to calculate blurring term EM update using covariance KDW20120522
				// -------- covariance --------------

				btk[k][i][j] = btkv[xt[i][j]][i][j];

				// Initialize kmat[][][][]
				for(l = 0; l < num_vpar; l++)
					kmat[k][l][i][j] = 0.;

				if((i>d)&&(i<rows-d-1)&&(j>d)&&(j<cols-d-1)){
					for(ii=0;ii<classes;ii++)
						for(jj=0;jj<classes;jj++){
							covh[ii][jj] = 0.;
							covv[ii][jj] = 0.;
						}
					for(ii=-d;ii<=d;ii++)
						for(jj=-d;jj<=d;jj++){
							kmat[k][xt[i-ii][j-jj]][i][j] += blur[ii+d][jj+d]*blur[ii+d][jj+d];
							covh[xt[i-ii][j-jj]][xt[i-ii-1][j-jj]] += blur[ii+d][jj+d]*blur[ii+d+1][jj+d];
							covh[xt[i-ii][j-jj]][xt[i-ii+1][j-jj]] += blur[ii+d][jj+d]*blur[ii+d-1][jj+d];
							covv[xt[i-ii][j-jj]][xt[i-ii][j-jj-1]] += blur[ii+d][jj+d]*blur[ii+d][jj+d+1];
							covv[xt[i-ii][j-jj]][xt[i-ii][j-jj+1]] += blur[ii+d][jj+d]*blur[ii+d][jj+d-1];
						}
					idx = classes;
					for(ii=0;ii<classes;ii++)
						for(jj=0;jj<classes;jj++){
							kmat[k][idx][i][j] = covh[ii][jj];
							kmat[k][idx+classes*classes][i][j] = covv[ii][jj];
							idx++;
						}
				}
#endif
			}

#ifdef MPM_IMG_OUTPUT
		// writing images
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		sprintf(outfileName, "emiter_%2d_mpmiter_%2d.tiff", emiter, k);
		if ((fp = fopen(outfileName, "wb")) == NULL ) {
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) {
			printf("Error writing TIFF file %s\n", outfileName);
			exit(1);
		}
	}
	fclose(fp);
	free_TIFF(&output_img);
#else
	}
#endif
	// to calculate blurring term EM update KDW20120222
	if(enable_blur){
		for(l=0;l<classes;l++)
			s_p[l] = numer_sum[l]/(double)pr_sum[l];
	}

	/* Normalize probabilities */
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			for (l = 0; l < classes; l++)
				pr[l][i][j] = pr[l][i][j] / (double)n;

	/* Clean Up */
	for (l = 0; l < classes; l++){
		free_img((void **)bt[l]);
		free_img((void **)btkv[l]);
		free_img((void **)yk[l]);
	}
}
