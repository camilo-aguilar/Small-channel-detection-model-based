/* em.c */
/** usage:  em beta infile outfile mpm-iter em-iter #_of_classes **/
#include "em.h"
#include "mpm.h"

void Importance_Map(unsigned char **gt, unsigned char **im, unsigned int rows, unsigned int cols, int l);
void Importance_Map2(unsigned char **gt, unsigned char **im, unsigned int rows, unsigned int cols, int l);
int MRF_Term(unsigned char **xt, unsigned int rows, unsigned int cols);
void entropy(double ***probs, unsigned char **output, unsigned int rows, unsigned int cols, int classes);
void blur(double **h, unsigned char **in, unsigned char **out, int rows, int cols, int mask_size);
double weighted_pmp(unsigned char **xt, unsigned char **gt, unsigned char **im,
		unsigned int rows, unsigned int cols, unsigned int classes, unsigned int level);
void update_parameter(unsigned char **y, double ***probs, double *m, double *v, double *N, int rows,
		int cols, int classes);

#ifdef DEF_COVARIANCE
int update_parameter2(unsigned char **y, unsigned char **xt, double ***btk, double ***probs, double **diff,
		double *s_p, double *m, double *v, double *N, double ****kmat, double *vmat,
		int mpmiter, int blur_size, int rows, int cols, int classes);
#else
int update_parameter2(unsigned char **y, unsigned char **xt, double ***btk, double ***probs, double **diff,
		double *s_p, double *m, double *v, double *N, int mpmiter, int blur_size,
		int rows, int cols, int classes);
#endif

/********************************************************************************************/
/*                                                                                          */
/*                                                                                          */
/********************************************************************************************/
int old_main(int argc,char *argv[]) {
	struct TIFF_img input_nblur_img, input_gt_img, input_sharp_img, input_blur_img, output_img;
	FILE *fp;
	int i, j, l, rows, cols, level, sum;
	double beta[MAX_CLASSES], gamma[MAX_CLASSES], ga;
	int mpmiter, emiter, classes, run_type;
	unsigned char **y, **xt, **xtr[REPET], **gt, **im_old, **im_new;  /* output : entropy image */

	double max_entropy = 0;
	char  infileName[1024] = "circles.tiff";
	char  outfilePrefix[1024] = "circles_";
	char  outfileName[1024] = "circles_out.tiff";
	char  gtfileName[1024] = "circles_gt.tiff";

	double misclassed;
	unsigned char xt255;
	double new_wpmp, old_wpmp;

	double **blur, di, dj, dsum, delta, sigma;
	int blur_size = 5, enable_blur;

	/* Verify correct number of arguments */
	if (argc != 8) {
		printf("usage:  %s run_type input_filename beta0 beta1 delta sigma classes\n", argv[0]); // 4
		printf("run_type:  0(make synthetic image), 1(segmentation)\n"); // 4
		exit(1);
	}
	run_type = atoi(argv[1]);
	sprintf(infileName, "%s.tiff",argv[2]);
	sprintf(gtfileName, "%s_gt.tiff",argv[2]);
	sprintf(outfilePrefix, "%s",argv[2]);
	beta[0] = atof(argv[3]);
	beta[1] = atof(argv[4]);
	delta = atof(argv[5]);
	sigma = atof(argv[6]);
	classes = atoi(argv[7]);
	ga = 0;
	mpmiter = 10;
	emiter = 30;

	max_entropy = log10((double)classes)/log10(2.0);

//	readseed();

	for(i = 0; i < MAX_CLASSES; i++)
		gamma[i] = 0;

	gamma[0] = ga;

	if ((fp = fopen(infileName, "rb")) == NULL) {
		printf("Cannot open file %s\n", infileName);

	}
	else {
		if (read_TIFF(fp, &input_blur_img)) { 	/* read image */
			printf("Error reading file %s\n", infileName);
			exit(1);
		}
		/* close image file */
		fclose(fp);

		/* check the type of image data */
		if (input_blur_img.TIFF_type != 'g') {
			printf("Error:  Image must be grayscale.\n");
			exit(1);
		}

		/* open ground truth image file */
		if ((fp = fopen(gtfileName, "rb")) == NULL) {
			printf("Cannot open file %s\n", gtfileName);
			if ((fp = fopen(infileName, "rb")) == NULL) {
				printf("Cannot open file %s\n", infileName);
				exit(1);
			}
		}

		/* read ground truth image */
		if (read_TIFF(fp, &input_gt_img)) {
			printf("Error reading file %s\n", gtfileName);
			exit(1);
		}

		/* close image file */
		fclose(fp);

		/* check the type of image data */
		if (input_gt_img.TIFF_type != 'g') {
			printf("Error:  Image must be grayscale.\n");
			exit(1);
		}
	}
	printf("File reading end..\n");
	cols = input_blur_img.width;
	rows = input_blur_img.height;
	printf ("cols = %d, rows = %d\n",cols, rows);

	/* Copy ground truth image to gt[][] */
	gt = (unsigned char **)get_img(cols, rows, sizeof(char));
	im_old = (unsigned char **)get_img(cols, rows, sizeof(char));
	im_new = (unsigned char **)get_img(cols, rows, sizeof(char));
	xt = (unsigned char **)get_img(cols, rows, sizeof(char));
	blur = (double **)get_img(blur_size, blur_size, sizeof(double));
	for (l = 0; l < REPET; l++)
		xtr[l] = (unsigned char **)get_img(cols, rows, sizeof(char));
	y = (unsigned char **)get_img(cols, rows, sizeof(char));
	/* Allocate space for the output image, and copy a scaled xt */
	get_TIFF(&output_img, rows, cols, 'g');

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			gt[i][j] = input_gt_img.mono[i][j];

	/* calculate blurring matrix */

	/*	Aexp(-x/sigma)											*/
	/*	A(amplitude), x(distance from center), sigma(width)		*/
	if(sigma == 0.)
		enable_blur = 0;
	else{
		enable_blur = 1;
		dsum = 0.;
		for (i = 0; i < blur_size; i++){
			for (j = 0; j < blur_size; j++){
				di = i-(blur_size-1)/2.;
				dj = j-(blur_size-1)/2.;
				blur[i][j] = exp(-sqrt(di*di+dj*dj)/sigma);
				dsum += blur[i][j];
			}
		}
		for (i = 0; i < blur_size; i++){
			for (j = 0; j < blur_size; j++){
				blur[i][j] = blur[i][j]/dsum;
				printf("%1.4f  ", blur[i][j]);
			}
			printf("\n");
		}
	}

	level = 7;
	// Old Importance Map
	sprintf(outfileName, "%s_im_old.tiff",outfilePrefix);
	Importance_Map(gt, im_old, rows, cols, level);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			output_img.mono[i][j] = (int)im_old[i][j] * 255 / (level - 1);
	if ((fp = fopen(outfileName, "wb")) == NULL ) {
		printf("Cannot open file %s\n", outfileName);
		exit(1);
	}
	if (write_TIFF(fp, &output_img)) {
		printf("Error writing TIFF file %s\n", outfileName);
		exit(1);
	}
	fclose(fp);
	// New Importance Map
	sprintf(outfileName, "%s_im_new.tiff",outfilePrefix);
	Importance_Map2(gt, im_new, rows, cols, level);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			output_img.mono[i][j] = (int)im_new[i][j] * 255 / (level - 1);
	if ((fp = fopen(outfileName, "wb")) == NULL ) {
		printf("Cannot open file %s\n", outfileName);
		exit(1);
	}
	if (write_TIFF(fp, &output_img)) {
		printf("Error writing TIFF file %s\n", outfileName);
		exit(1);
	}
	fclose(fp);

	if(run_type == 0){
		/* no blur */
		printf("emmpm no blur image..\n");
		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++)
				y[i][j] = input_nblur_img.mono[i][j];

		for (l = 0; l < REPET; l++)
			emmpm(y, xtr[l], beta, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);
		// This works only for 2 classes
		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++){
				sum = 0;
				for (l = 0; l < REPET; l++){
					sum += xtr[l][i][j];
				}
				if (sum>REPET/2)
					xt[i][j] = 1;
				else
					xt[i][j] = 0;
			}

		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		sprintf(outfileName, "%s_b%1.2f_b%1.2f_d%1.3f_s%1.2f_nblur.tiff",outfilePrefix, beta[0], beta[1], delta, sigma);
		if ((fp = fopen(outfileName, "wb")) == NULL ) {
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) {
			printf("Error writing TIFF file %s\n", outfileName);
			exit(1);
		}
		fclose(fp);
	}

	/* blur */
	printf("emmpm blur image..\n");
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			y[i][j] = input_blur_img.mono[i][j];

	for (l = 0; l < REPET; l++)
		emmpm(y, xtr[l], beta, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);
	// This works only for 2 classes
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			sum = 0;
			for (l = 0; l < REPET; l++){
				sum += xtr[l][i][j];
			}
			if (sum>REPET/2)
				xt[i][j] = 1;
			else
				xt[i][j] = 0;
		}

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
	sprintf(outfileName, "%s_b%1.2f_b%1.2f_d%1.3f_s%1.2f.tiff",outfilePrefix, beta[0], beta[1], delta, sigma);
	if ((fp = fopen(outfileName, "wb")) == NULL ) {
		printf("Cannot open file %s\n", outfileName);
		exit(1);
	}
	if (write_TIFF(fp, &output_img)) {
		printf("Error writing TIFF file %s\n", outfileName);
		exit(1);
	}
	fclose(fp);

	/* Clean up */
	free_TIFF(&input_nblur_img);
	free_TIFF(&input_blur_img);
	free_TIFF(&input_sharp_img);
	free_TIFF(&output_img);

	// PMP
	misclassed = 0;
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			xt255 = (unsigned char)xt[i][j] * 255 / (classes - 1);
			if(xt255!=gt[i][j]) misclassed += 1.0;
		}
	misclassed = misclassed*100/(cols*rows); // percentage of misclassified pixels

	// Weighted PMP
	old_wpmp = weighted_pmp(xt, gt, im_old, rows, cols, classes, level);
	new_wpmp = weighted_pmp(xt, gt, im_new, rows, cols, classes, level);

	printf("misclassed = %1.4f, wpmp = %1.4f, new wpmp = %1.4f\n", misclassed, old_wpmp, new_wpmp);
	sprintf(outfileName, "%s_b%1.2f_b%1.2f_d%1.3f_s%1.2f_p%1.2f_wp%1.2f_nwp%1.2f.txt",outfilePrefix, beta[0], beta[1], delta, sigma, misclassed, old_wpmp, new_wpmp);
	if ((fp = fopen(outfileName, "wb")) == NULL ) {
		printf("Cannot open file %s\n", outfileName);
		exit(1);
	}
	fprintf(fp,"misclassed = %1.4f, wpmp = %1.4f, new wpmp = %1.4f\r\n", misclassed, old_wpmp, new_wpmp);
	fclose(fp);

	free_img((void **)gt);
	free_img((void **)im_old);
	free_img((void **)im_new);
	free_img((void **)xt);
	free_img((void **)y);
	for (l = 0; l < REPET; l++)
		free_img((void **)xtr[l]);

//	writeseed();

	return 0;
}

void emmpm(unsigned char **y, unsigned char **xt, double *beta, double *gamma, int emiter, int mpmiter,
		int rows, int cols, int classes, double **blur, int blur_size, int enable_blur)
{
	int i, j;
	double x, m[MAX_CLASSES], v[MAX_CLASSES], N[MAX_CLASSES], s_p[MAX_CLASSES], **btk[MAX_MPMITER];
	/*	m[l] - estimate of mean for class l
		v[l] - estimate of variance for class l
		N[l] - Used for normalization
	*/
	double mu, sigma;
	int k, kk, l;
	int num_vpar = 2*classes*classes+classes;
	double **probs[MAX_CLASSES], **diff;
	double tmp;
#ifdef EM_IMG_OUTPUT
	struct TIFF_img output_img;
	char  outfileName[1024];
	FILE* fp;

	get_TIFF(&output_img, rows, cols, 'g');
#endif

	for (l = 0; l < mpmiter; l++)
		btk[l] = (double **)get_img(cols, rows, sizeof(double));
	for (l = 0; l < classes; l++) {
		v[l] = 20;
		probs[l] = (double **)get_img(cols, rows, sizeof(double));
	}
#ifdef DEF_COVARIANCE
	double **kmat[mpmiter][num_vpar], vmat[num_vpar];
	for (k = 0; k < mpmiter; k++)
		for(l = 0; l < num_vpar; l++)
			kmat[k][l] = (double **)get_img(cols, rows,  sizeof(double));
	for (l = 0; l < classes; l++) {
		vmat[l] = v[l];
	}
	for (l = classes; l < num_vpar; l++) {
		vmat[l] = 0;
	}
#endif
	tmp = 0.;
	for (i=0; i<blur_size; i++)
		for (j=0; j<blur_size; j++)
			tmp += blur[i][j]*blur[i][j];
	tmp *= v[0];
//	for (i = 0; i < rows; i++)
//		for (j = 0; j < cols; j++)
//			kv[mpmiter-1][i][j] = tmp;

	diff = (double **)get_img(cols, rows, sizeof(double));

	/* Initialize classification of each pixel randomly with a uniform distribution */
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++) {
			x = random2();
			l = 0;
			while ((double)(l + 1) / classes <= x)  // may incur l = classes when x = 1
				l++;
			xt[i][j] = l;
		}
	/* Initialization of parameter estimation */
	mu = 0;
	sigma = 0;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			mu += y[i][j];
	mu /= rows*cols;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			sigma += (y[i][j]-mu)*(y[i][j]-mu);
	sigma /= rows*cols;
	sigma = sqrt((double)sigma);
	printf("mu=%f sigma=%f\n",mu,sigma);

	if (classes%2 == 0)
	{
		for (k=0; k<classes/2; k++)
		{
			m[classes/2 + k] = mu + (k+1)*sigma/2;
			m[classes/2 - 1 - k] = mu - (k+1)*sigma/2;
		}
	}
	else
	{
		m[classes/2] = mu;
		for (k=0; k<classes/2; k++)
		{
			m[classes/2 + 1 + k] = mu + (k+1)*sigma/2;
			m[classes/2 - 1 - k] = mu - (k+1)*sigma/2;
		}
	}

	/* Perform EM */
	for (k = 0; k < emiter; k++) {
//		printf("emiter = %d\n",k);
		/* Perform MPM */
#ifdef DEF_COVARIANCE
		mpm(y, xt, btk, probs, beta, gamma, m, v, rows, cols, mpmiter, classes, k,
				blur, blur_size, enable_blur, diff, s_p, kmat, vmat);

		if (enable_blur){
			if(!update_parameter2(y, xt, btk, probs, diff, s_p, m, v, N, kmat, vmat, mpmiter,
									blur_size, rows, cols, classes))
				printf("singular matrix\n");
		}
		else
			update_parameter(y, probs, m, v, N, rows, cols, classes);
#else
		mpm(y, xt, btk, probs, beta, gamma, m, v, rows, cols, mpmiter, classes, k,
				blur, blur_size, enable_blur, diff, s_p);

		if (enable_blur){
			if(!update_parameter2(y, xt, btk, probs, diff, s_p, m, v, N, mpmiter,
									blur_size, rows, cols, classes))
				printf("singular matrix\n");
		}
		else
			update_parameter(y, probs, m, v, N, rows, cols, classes);
#endif

		/* Monitor estimates of mean and variance */
//		if (emiter < 10 || (k + 1) % (emiter / 10) == 0) {
//			for (l = 0; l < classes - 1; l++)
//				printf("%.3f %.3f ", m[l], v[l]);
//			printf("%.3f %.3f\n", m[classes - 1], v[classes - 1]);
//		}

		/* Eliminate any classes that have zero probability */
		for (kk = 0; kk < classes; kk++)
			if (N[kk] == 0) {
				for (l = kk; l < classes - 1; l++) {
					/* Move other classes to fill the gap */
					N[l] = N[l + 1];
					m[l] = m[l + 1];
					v[l] = v[l + 1];
					for (i = 0; i < rows; i++)
					for (j = 0; j < cols; j++)
						if (xt[i][j] == l + 1)
							xt[i][j] = l;
				}
				classes = classes - 1;  // push the eliminated class into the last class
			}

		// writing images
#ifdef EM_IMG_OUTPUT
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		sprintf(outfileName, "emiter_%2d.tiff", k);
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

	for (l = 0; l < mpmiter; l++)
		free_img((void **)btk[l]);
#ifdef DEF_COVARIANCE
	for (k = 0; k < mpmiter; k++)
		for(l = 0; l < num_vpar; l++)
			free_img((void **)kmat[k][i]);
#endif
	for (l = 0; l < classes; l++)
		free_img((void **)probs[l]);
	free_img((void **)diff);

	return;
}

void emmpm_betaimg(unsigned char **y, unsigned char **xt, double **beta[], double *gamma, int emiter, int mpmiter,
		int rows, int cols, int classes, double **blur, int blur_size, int enable_blur)
{
	int i, j;
	double x, m[MAX_CLASSES], v[MAX_CLASSES], N[MAX_CLASSES], s_p[MAX_CLASSES], **btk[MAX_MPMITER];
	/*	m[l] - estimate of mean for class l
		v[l] - estimate of variance for class l
		N[l] - Used for normalization
	*/
	double mu, sigma;
	int k, kk, l;
	double **probs[MAX_CLASSES], **diff;
	double tmp;
#ifdef EM_IMG_OUTPUT
	struct TIFF_img output_img;
	char  outfileName[1024];
	FILE* fp;

	get_TIFF(&output_img, rows, cols, 'g');
#endif

	for (l = 0; l < mpmiter; l++)
		btk[l] = (double **)get_img(cols, rows, sizeof(double));
	for (l = 0; l < classes; l++) {
		v[l] = 20;
		probs[l] = (double **)get_img(cols, rows, sizeof(double));
	}
#ifdef DEF_COVARIANCE
	double **kmat[mpmiter][num_vpar], vmat[num_vpar];
	for (k = 0; k < mpmiter; k++)
		for(l = 0; l < num_vpar; l++)
			kmat[k][l] = (double **)get_img(cols, rows,  sizeof(double));
	for (l = 0; l < classes; l++) {
		vmat[l] = v[l];
	}
	for (l = classes; l < num_vpar; l++) {
		vmat[l] = 0;
	}
#endif
	tmp = 0.;
	for (i=0; i<blur_size; i++)
		for (j=0; j<blur_size; j++)
			tmp += blur[i][j]*blur[i][j];
	tmp *= v[0];
//	for (i = 0; i < rows; i++)
//		for (j = 0; j < cols; j++)
//			kv[mpmiter-1][i][j] = tmp;

	diff = (double **)get_img(cols, rows, sizeof(double));

	/* Initialize classification of each pixel randomly with a uniform distribution */
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++) {
			x = random2();
			l = 0;
			while ((double)(l + 1) / classes <= x)  // may incur l = classes when x = 1
				l++;
			xt[i][j] = l;
		}
	/* Initialization of parameter estimation */
	mu = 0;
	sigma = 0;
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			mu += y[i][j];
	mu /= rows*cols;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			sigma += (y[i][j]-mu)*(y[i][j]-mu);
	sigma /= rows*cols;
	sigma = sqrt((double)sigma);
	printf("mu=%f sigma=%f\n",mu,sigma);

	if (classes%2 == 0)
	{
		for (k=0; k<classes/2; k++)
		{
			m[classes/2 + k] = mu + (k+1)*sigma/2;
			m[classes/2 - 1 - k] = mu - (k+1)*sigma/2;
		}
	}
	else
	{
		m[classes/2] = mu;
		for (k=0; k<classes/2; k++)
		{
			m[classes/2 + 1 + k] = mu + (k+1)*sigma/2;
			m[classes/2 - 1 - k] = mu - (k+1)*sigma/2;
		}
	}

	/* Perform EM */
	for (k = 0; k < emiter; k++) {
//		printf("emiter = %d\n",k);
		/* Perform MPM */
#ifdef DEF_COVARIANCE
		mpm_betaimg(y, xt, btk, probs, beta, gamma, m, v, rows, cols, mpmiter, classes, k,
				blur, blur_size, enable_blur, diff, s_p, kmat, vmat);

		if (enable_blur){
			if(!update_parameter2(y, xt, btk, probs, diff, s_p, m, v, N, kmat, vmat, mpmiter,
									blur_size, rows, cols, classes))
				printf("singular matrix\n");
		}
		else
			update_parameter(y, probs, m, v, N, rows, cols, classes);
#else
		mpm_betaimg(y, xt, btk, probs, beta, gamma, m, v, rows, cols, mpmiter, classes, k,
				blur, blur_size, enable_blur, diff, s_p);

		if (enable_blur){
			if(!update_parameter2(y, xt, btk, probs, diff, s_p, m, v, N, mpmiter,
									blur_size, rows, cols, classes))
				printf("singular matrix\n");
		}
		else
			update_parameter(y, probs, m, v, N, rows, cols, classes);
#endif

		/* Monitor estimates of mean and variance */
//		if (emiter < 10 || (k + 1) % (emiter / 10) == 0) {
//			for (l = 0; l < classes - 1; l++)
//				printf("%.3f %.3f ", m[l], v[l]);
//			printf("%.3f %.3f\n", m[classes - 1], v[classes - 1]);
//		}

		/* Eliminate any classes that have zero probability */
		for (kk = 0; kk < classes; kk++)
			if (N[kk] == 0) {
				for (l = kk; l < classes - 1; l++) {
					/* Move other classes to fill the gap */
					N[l] = N[l + 1];
					m[l] = m[l + 1];
					v[l] = v[l + 1];
					for (i = 0; i < rows; i++)
					for (j = 0; j < cols; j++)
						if (xt[i][j] == l + 1)
							xt[i][j] = l;
				}
				classes = classes - 1;  // push the eliminated class into the last class
			}

		// writing images
#ifdef EM_IMG_OUTPUT
		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		sprintf(outfileName, "emiter_%2d.tiff", k);
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

	for (l = 0; l < mpmiter; l++)
		free_img((void **)btk[l]);
#ifdef DEF_COVARIANCE
	for (k = 0; k < mpmiter; k++)
		for(l = 0; l < num_vpar; l++)
			free_img((void **)kmat[k][i]);
#endif
	for (l = 0; l < classes; l++)
		free_img((void **)probs[l]);
	free_img((void **)diff);

	return;
}

void update_parameter(unsigned char **y, double ***probs, double *m, double *v, double *N, int rows, int cols, int classes)
{
	int i,j,l;

	/* Reset model parameters to zero */
	for (l = 0; l < classes; l++) {
		m[l] = 0;
		v[l] = 0;
		N[l] = 0;
	}

	/*** Some efficiency was sacrificed for readability below ***/

	/* Update estimates for mean of each class */
	for (l = 0; l < classes; l++) {
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				N[l] += probs[l][i][j];  // denominator of (20)
				m[l] += y[i][j] * probs[l][i][j];  // numerator of (20)
			}
		}
		if (N[l] != 0)
			m[l] = m[l] / N[l];  // Eq. (20)
	}

//	for (l = 0; l < classes; l++) {
//		printf("m[%d] = %2.4f \n",l,m[l]);
//	}
	/* Update estimates of variance of each class */
	for (l = 0; l < classes; l++) {
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++)
				// numerator of (21)
				v[l] += (y[i][j] - m[l]) * (y[i][j] - m[l]) * probs[l][i][j];
		}
		if (N[l] != 0)
			v[l] = v[l] / N[l];
	}

}

#ifdef DEF_COVARIANCE
int update_parameter2(unsigned char **y, unsigned char **xt, double ***btk, double ***probs, double **diff,
		double *s_p, double *m, double *v, double *N, double ****kmat, double *vmat,
		int mpmiter, int blur_size, int rows, int cols, int classes)
{
	double **kd[MAX_CLASSES], gmat[num_vpar], gnorm, dmat[num_vpar];
	double alpha, alpha_new, dQ, ddQ;
	double **kv[MAX_CLASSES];
#else
int update_parameter2(unsigned char **y, unsigned char **xt, double ***btk, double ***probs, double **diff,
		double *s_p, double *m, double *v, double *N, int mpmiter, int blur_size,
		int rows, int cols, int classes)
{
#endif
	int i,j,l,success = 0;
	double tmp;
	int num_vpar = 2*classes*classes+classes, d = (blur_size-1)/2;

//	printf("enable blur EM update\n");
//	for(i=0;i<classes;i++){
//		for(j=0;j<classes+1;j++){
//			printf("%1.2f  ",diff[i][j] );
//		}
//		printf("\n" );
//	}

	/* Reset model parameters to zero */
	for (l = 0; l < classes; l++) {
//		m[l] = 0;
//		v[l] = 0;
		N[l] = 0;
	}

	/* Update estimates of means of each class */
	// Gauss Elimination for diff matrix
	// top-down
	for (i=0;i<classes;i++){
		if(diff[i][i]==0)
			for (j=i+1;j<classes;j++){
				if(diff[j][i]!=0.){
					for(l=i;l<classes+1;l++){
						tmp = diff[i][l];
						diff[i][l] = diff[j][l];
						diff[j][l] = tmp;
					}
					break;
				}
				if(j==classes) return 0;	// singular matrix
			}
		tmp = diff[i][i];
		for(l=i;l<classes+1;l++)
			diff[i][l] = diff[i][l]/tmp;
		if(i<classes-1){
			for (j=i+1;j<classes;j++){
				if (diff[j][i]!=0.){
					tmp = diff[j][i];
					for(l=i;l<classes+1;l++){
						diff[j][l] -= tmp*diff[i][l];
					}
				}
			}
		}
	}
	// bottom-up
	for (i=classes-1;i>0;i--){
		if(diff[i][i]==0)
			return 0;						// singular matrix
		for (j=i-1;j>=0;j--){
			if (diff[j][i]!=0.){
				tmp = diff[j][i];
				for(l=i;l<classes+1;l++){
					diff[j][l] -= tmp*diff[i][l];
				}
			}
		}
	}
	// solution
	for (l = 0; l < classes; l++) {
		m[l] = diff[l][classes];
//		printf("m[%d] = %2.4f \n",l,m[l]);
	}
	for (l = 0; l < classes; l++) {
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {
				N[l] += probs[l][i][j];  // denominator of (20)
			}
		}
	}
	/* Update estimates of variance of each class */
	// Single variance
	for (l = 0; l < classes; l++) {
		v[l] = s_p[l];
	}

#ifdef DEF_COVARIANCE
	for (k = 0; k < mpmiter; k++)
		kv[k] = (double **)get_img(cols, rows, sizeof(double));

	for (k = 0; k < mpmiter; k++)
		kd[k] = (double **)get_img(cols, rows, sizeof(double));

	// Covariance : using CGA(conjugate gradient algorithm)
	// 1.
	for (l = 0; l < mpmiter; l++){
		gmat[l] = 0;
	}
	for(k=0;k<mpmiter;k++)
		for(i=d;i<rows-d-1;i++)
			for(j=d;j<cols-d-1;j++){
				kv[k][i][j] = 0.;
				for(l=d;l<num_vpar;l++)
					kv[k][i][j] += kmat[k][l][i][j]*vmat[l];
			}
	for(i=d;i<rows-d-1;i++)
		for(j=d;j<cols-d-1;j++){
			for(k=0;k<mpmiter;k++){
				for(l=0;l<num_vpar;l++)
					gmat[l] += (btk[k][i][j] - kv[k][i][j])/(kv[k][i][j]*kv[k][i][j])*kmat[k][l][i][j];
			}
		}

	gnorm = 0;
	for(l=0;l<num_vpar;l++){
		gmat[l] = gmat[l]/(2*mpmiter);
		dmat[l] = -gmat[l];
		gnorm += gmat[l]*gmat[l];
	}
	if (gnorm < 0.01) return 1;
	for (iter0 = 0; iter0<10 ;iter0++){
		alpha_new = 0;
		dQ = 0;
		ddQ = 0;
		for(iter1 = 0; iter1<10 ;iter1++){
			alpha = alpha_new;
			for(l=0;l<num_vpar;l++)
				vmat[l] += alpha*dmat[l];
			for(k=0;k<mpmiter;k++)
				for(i=d;i<rows-d-1;i++)
					for(j=d;j<cols-d-1;j++){
						kv[k][i][j] = 0.;
						for(l=d;l<num_vpar;l++)
							kv[k][i][j] += kmat[k][l][i][j]*vmat[l];
					}
		}

	}

	for (k = 0; k < mpmiter; k++)
		free_img((void **)kd[k]);
	for (k = 0; k < mpmiter; k++)
		free_img((void **)kv[k]);
#endif

	return 1;
}

double weighted_pmp(unsigned char **xt, unsigned char **gt, unsigned char **im,
		unsigned int rows, unsigned int cols, unsigned int classes, unsigned int level)
{
	unsigned char xt255;
	int i, j;
	int k, l;
	unsigned int nlkk[MAX_LEVEL][MAX_CLASSES], nl[MAX_LEVEL], nlkk_sum;
	double wis, sum, dlevel;

	dlevel = (double)level;
	for (l = 0; l < (int)level; l++){
		for (k = 0; k < (int)classes; k++){
			nlkk[l][k] = 0;
		}
		nl[l] = 0;
	}

	for (i = 0; i < (int)rows; i++)
		for (j = 0; j < (int)cols; j++){
			nl[(int)im[i][j]]++;
			xt255 = (unsigned char)xt[i][j] * 255 / (classes - 1);
			if(xt255 == gt[i][j]) nlkk[(int)im[i][j]][xt[i][j]]++;
			//else printf("(%d,%d),",xt255, gt[i][j]);
		}

	sum = 0.;
	for (l = 0; l < (int)level; l++){
		wis = 0.5/dlevel+(double)l/(dlevel*(dlevel-1.0));
		nlkk_sum = 0;
		for (k = 0; k < (int)classes; k++){
			nlkk_sum += nlkk[l][k];
		}
		//printf("wis = %1.3f, nl[%d] = %d, nlkk_sum = %d\n", wis, l, nl[l], nlkk_sum);
		sum += wis/(double)nl[l]*(double)nlkk_sum;
	}
	return (1. - sum)*100.;
}


void Importance_Map(unsigned char **gt, unsigned char **im,
		unsigned int rows, unsigned int cols, int l)
{
	int i, j, ii, jj, tmp;

	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			im[i][j] = 0;
		}
	}
	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			if(i==0){
				if(j==0){
					if((gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
				else if(j==cols-1){
					if((gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
				else{
					if((gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j+1])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
			}
			else if(i==rows-1){
				if(j==0){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
				else if(j==cols-1){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
				else{
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i][j+1])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
			}
			else{
				if(j==0){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
				else if(j==cols-1){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
				else{
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j-1])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
			}
		}
	}

	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			if(im[i][j]==l-1){
				for (ii=1-l; ii<l; ii++){
					for (jj=1-l; jj<l; jj++) {
						if((i+ii>=0)&&(i+ii<(int)rows)&&(j+jj>=0)&&(j+jj<(int)cols)){
							tmp = l-1-abs(ii)-abs(jj);
							if(tmp>im[i+ii][j+jj])
								im[i+ii][j+jj] = tmp;
						}
					}
				}
			}
		}
	}
}


void Importance_Map2(unsigned char **gt, unsigned char **im,
		unsigned int rows, unsigned int cols, int l)
{
	int i, j, ii, jj;
	double **im2, dtmp;

	im2 = (double **)get_img(cols, rows, sizeof(double));

	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			im[i][j] = 0;
			im2[i][j] = 0;
		}
	}
	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			if(i==0){
				if(j==0){
					if((gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
				else if(j==cols-1){
					if((gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
				else{
					if((gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j+1])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
			}
			else if(i==rows-1){
				if(j==0){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
				else if(j==cols-1){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
				else{
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i][j+1])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
			}
			else{
				if(j==0){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
				else if(j==cols-1){
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j-1])){
							im[i][j] = l-1;
					}
				}
				else{
					if((gt[i][j]!=gt[i-1][j])||
						(gt[i][j]!=gt[i+1][j])||
						(gt[i][j]!=gt[i][j-1])||
						(gt[i][j]!=gt[i][j+1])){
							im[i][j] = l-1;
					}
				}
			}
		}
	}

	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			if(im[i][j]==l-1){
				for (ii=1-l; ii<l; ii++){
					for (jj=1-l; jj<l; jj++) {
						if((i+ii>=0)&&(i+ii<(int)rows)&&(j+jj>=0)&&(j+jj<(int)cols)){
							dtmp = (double)(l-1-abs(ii)-abs(jj));
							if(dtmp>0) im2[i+ii][j+jj] += dtmp;
						}
					}
				}
			}
		}
	}

	// find max value
	dtmp = 0;
	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			if(im2[i][j]>dtmp) dtmp = im2[i][j];
		}
	}
	// normalize
	for (i=0; i<(int)rows; i++){
		for (j=0; j<(int)cols; j++) {
			im[i][j] = (unsigned char)((im2[i][j]*(double)(l-1))/dtmp);
//			if((i%4==0)&&(j%4==0)) printf("%d, ",im[i][j]);
		}
	}
	free_img((void **)im2);
}
