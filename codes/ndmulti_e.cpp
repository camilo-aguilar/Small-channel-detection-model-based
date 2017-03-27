#include "neck.h"
#include "ndmulti_e.h"
#include "random.h"
#include "QualityCandy.h"


void Qdelete(NeckDent *mp, int k, int *np_num, double **inter_e)
{
	int i,j;

	for (i = k; i < (*np_num)-1 ; i++){
			mp[i] = mp[i+1];	// mp[i].multiple_E = mp[i+1].multiple_E; included
			for (j = 0; j < *np_num ; j++){
				inter_e[i][j] = inter_e[i+1][j];
			}
	}
	for (j = 0; j < (*np_num) ; j++){
		inter_e[(*np_num)-1][j] = 0; // i = (*np_num)-1
	}
	mp[(*np_num)-1].multiple_E = 0;

	for (i = 0; i < (*np_num)-1 ; i++){
		mp[i].multiple_E -= inter_e[i][k];
		for (j = k; j < (*np_num)-1 ; j++){
			inter_e[i][j] = inter_e[i][j+1];
		}
		inter_e[i][j] = 0; // j = (*np_num)-1
	}
	(*np_num)--;
}


int compare_Single_E (const void * a, const void * b)
{
	if(((NeckDent*)a)->single_E > ((NeckDent*)b)->single_E) 
		return -1;
	else if(((NeckDent*)a)->single_E == ((NeckDent*)b)->single_E) 
		return 0;
	else //if(((NeckDent*)a)->single_E < ((NeckDent*)b)->single_E) 
		return 1;
}


int check_inside (DPoint a, DPoint b, DPoint c, DPoint d, DPoint p)
{
	if(( ( (p.y-a.y)*(a.x-b.x) - (p.x-a.x)*(a.y-b.y) ) *((p.y-c.y)*(c.x-d.x)-(p.x-c.x)*(c.y-d.y))<=0)
		&&(((p.y-a.y)*(a.x-c.x)-(p.x-a.x)*(a.y-c.y))*((p.y-b.y)*(b.x-d.x)-(p.x-b.x)*(b.y-d.y))<=0))
		return 1;
	else
		return 0;
}


double avg_Gaussian_save(unsigned char **y, int cols, int rows, NeckDent *mp, MPP_Parameters mpp, int *num)
//						 double mean0, double vari0, double mean1, double vari1, 
//					double variance, double gaussian_tau, double t,  
//					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, , double *e0, double *e1, double *e2, double *e3
{
	int num2;
	double amp, moffset;
	double di, dj, adi, adj;
	double w_5, l_5, mt_ww;
	double y2, dtmp;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, error2;
	DPoint pt;
	FILE *fp1, *fp2;
	char filename[100];

	sprintf(filename, "channel%d_intensity_fn.txt",mp->type);
	if ((fp1 = fopen(filename, "wb")) == NULL ) {
		printf("Cannot open file %s\n", filename);
		exit(1);
	}
	sprintf(filename, "channel%d_data.txt",mp->type);
	if ((fp2 = fopen(filename, "wb")) == NULL ) {
		printf("Cannot open file %s\n", filename);
		exit(1);
	}

	num2 = 0;
	w_5 = mp->width/2.;
	l_5 = mp->length/2.;
	mt_ww = -mpp.gaussian_tau/(mp->width*mp->width);

	(*num) = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			adj = fabs(dj);
			if(mp->type==NDTYPE_BOTH){
				y2 = (di)*(di)*mt_ww;
			}
			else if(mp->type==NDTYPE_NECKING){
				if(adj<l_5-w_5){
					y2 = di*di*mt_ww;
				}
				else{
					dtmp = distance1(l_5-w_5,w_5,adj,adi);
					if(dtmp<w_5){
						y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
					}
					else{//4
						y2 = 0;
					}
				}
			}
			else {
				if(dj>0){
					if(dj<l_5-w_5){
						y2 = di*di*mt_ww;
					}
					else{
						dtmp = distance1(l_5-w_5,w_5,dj,adi);
						if(dtmp<w_5){
							y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
						}
						else{//4
							y2 = 0;
						}
					}
				}
				else{
					if(dj>-l_5+w_5){//1
						y2 = (di)*(di)*mt_ww;
					}
					else{//3
						y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
					}
				}
			}
			y2 = -exp(y2);
			pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
			pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				fprintf(fp1,"%1.2f, ", y2);
				fprintf(fp2,"%1.2f, ", dtmp);
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp*dtmp;
				(*num)++;
			}
			num2++;
		}
		fprintf(fp1,"\r\n");
		fprintf(fp2,"\r\n");
	}
	fprintf(fp1,"\r\n");
	fprintf(fp2,"\r\n");
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		mp->e[0] = 0;
		mp->e[1] = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	fprintf(fp1,"amp, %1.2f\r\n", amp);
	fprintf(fp1,"offset, %1.2f\r\n", moffset);
	fprintf(fp1,"error, %1.2f\r\n", error2);
	mp->e[0] = amp;
	mp->e[1] = error2;
	mp->e[2] = moffset;
	fclose(fp1);
	fclose(fp2);
	return error2;
}


double avg_Gaussian_save2(unsigned char **y, int cols, int rows, NeckDent *mp, MPP_Parameters mpp, int *num)
//						 double mean0, double vari0, double mean1, double vari1, 
//					double variance, double gaussian_tau, double t,  
//					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, , double *e0, double *e1, double *e2, double *e3
{
	int num2;
	double amp, moffset;
	double di, dj;
	double w_5, l_5, mt_ww;
	double y2, dtmp;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, error2;
	DPoint pt;
	FILE *fp1, *fp2;
	char filename[100];

	sprintf(filename, "channel%d_intensity_fn.txt",mp->type);
	if ((fp1 = fopen(filename, "wb")) == NULL ) {
		printf("Cannot open file %s\n", filename);
		exit(1);
	}
	sprintf(filename, "channel%d_data.txt",mp->type);
	if ((fp2 = fopen(filename, "wb")) == NULL ) {
		printf("Cannot open file %s\n", filename);
		exit(1);
	}

	num2 = 0;
	w_5 = mp->width/2.;
	l_5 = mp->length/2.;
	mt_ww = -mpp.gaussian_tau/(mp->width*mp->width);

	(*num) = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			if(mp->type==NDTYPE_BOTH){
				y2 = (di)*(di)*mt_ww;
			}
			else if(mp->type==NDTYPE_NECKING){
				if(dj>0){
					if(l_5-dj>fabs(di)){//2
						y2 = (di)*(di)*mt_ww;
					}
					else{//4
						y2 = (l_5-dj)*(l_5-dj)*mt_ww;
					}
				}
				else{
					if(l_5+dj>fabs(di)){//1
						y2 = (di)*(di)*mt_ww;
					}
					else{//3
						y2 = (l_5+dj)*(l_5+dj)*mt_ww;
					}
				}
			}
			else {
				if(dj>0){
					if(l_5-dj>fabs(di)){//2
						y2 = (di)*(di)*mt_ww;
					}
					else{//4
						y2 = (l_5-dj)*(l_5-dj)*mt_ww;
					}
				}
				else{
					if(dj>-l_5+w_5){//1
						y2 = (di)*(di)*mt_ww;
					}
					else{//3
						y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
					}
				}
			}
			y2 = -exp(y2);
			pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
			pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				fprintf(fp1,"%1.2f, ", y2);
				fprintf(fp2,"%1.2f, ", dtmp);
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp*dtmp;
				(*num)++;
			}
			num2++;
		}
		fprintf(fp1,"\r\n");
		fprintf(fp2,"\r\n");
	}
	fprintf(fp1,"\r\n");
	fprintf(fp2,"\r\n");
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		mp->e[0] = 0;
		mp->e[1] = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	fprintf(fp1,"amp, %1.2f\r\n", amp);
	fprintf(fp1,"offset, %1.2f\r\n", moffset);
	fprintf(fp1,"error, %1.2f\r\n", error2);
	mp->e[0] = amp;
	mp->e[1] = error2;
	mp->e[2] = moffset;
	fclose(fp1);
	fclose(fp2);
	return error2;
}

// Channel Object adaptive mean map (calculate min and max to calculate mean map)
// return single_energy
double avg_Gaussian(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3)
{
	int n1, num2;
	double amp, moffset;
	double di, dj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	num2 = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			y2 = (di)*(di)*mt_ww;
			y2 = (-exp(y2));
			//y2 = -exp(y2);
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				dtmp2 = dtmp*dtmp;
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp2;
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;
	*e3 = error_s;

	return error2;
}


//#define DELTA_AMP		1.0
//#define DELTA_OFFSET	1.0

#define GAUSS_TH -0.2 //-0.2
/***********************************************************************/
//
//	n1n1n1n1n1n1n1n1n1n1n1
//						
//	c2c2c2c2c2c2c2c2c2c2c2
//						
//	n2n2n2n2n2n2n2n2n2n2n2
//
/***********************************************************************/
double avg_t_test(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double error_th, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e)
{
	int i, num_out = 0, num_s = 0;
	int num_ch = 0, num_nch = 0;
	int num_1ch = 0, num_2ch = 0, num_3ch = 0;
	int num_n1ch = 0, num_n2ch = 0;
	double di, dj, w, l;
	double w_5, l_5, tau = gaussian_tau;
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
	double sin_t = sin(t), cos_t = cos(t), di_sin_t, di_cos_t;
	double cx, cy, dj_cos_t[100], dj_sin_t[100];
	DPoint pt;
	int flag;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	w_5 = w/2.;
	l_5 = l/2.;
	cx = (pt1.x+pt4.x)/2.;
	cy = (pt1.y+pt4.y)/2.;

	(*num) = 0;

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
			sum_n1ch += dtmp;
			sum_n1ch2 += dtmp*dtmp;
			num_n1ch++;
			(*num)++;
			flag = 1;
		}
		else
			num_out++;
		// non channel 2
		pt.x = cx + dj_cos_t[i] - di_sin_t;
		pt.y = cy + dj_sin_t[i] + di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp2 = real_coord(y,pt.y,pt.x);
			sum_n2ch += dtmp2;
			sum_n2ch2 += dtmp2*dtmp2;
			num_n2ch++;
			(*num)++;
			if(flag){
				sum_s += fabs(dtmp-dtmp2);
				num_s++;
			}
		}
		else
			num_out++;

	}
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
			(*num)++;
		}
		else
			num_out++;
	}

	// Calculate distance1
	if(((double)(*num))/(double)(num_out+*num)<OBJECT_IN_RATIO){
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

	if(num_ch){
		sum_ch = sum_ch/(double)num_ch;
		sum_ch2 = sum_ch2/(double)num_ch;
	}
	if(num_nch){
		sum_nch = sum_nch/(double)num_nch;
		sum_nch2 = sum_nch2/(double)num_nch;
	}
	if(*num){
		sum_all = sum_all/(double)(num_ch+num_nch);
		sum_all2 = sum_all2/(double)(num_ch+num_nch);
	}
	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel
	sum_all2 = sum_all2 - sum_all*sum_all; // variance of the outside of a chanel
	if(*num !=0){
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		t_test1 = fabs(sum_ch-sum_nch)/(sqrt(sum_ch2/(double)num_ch+sum_nch2/(double)num_nch));
		//t_test2 = fabs(sum_n1ch-sum_n2ch)/(sqrt(sum_n1ch2/(double)num_n1ch+sum_n2ch2/(double)num_n2ch));
		t_test = t_test1/fmax(SYM_TH,sum_s/sum_all)*SYM_TH;//*sqrt(sum_all2)

		if ((sum_nch>sum_ch)&&(sum_all2>150))
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
			likely = pow(likely,1.0/3.0);
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/3.0);

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


/***********************************************************************/
//  index pixels
//	 		n1n1n1n1n1n1n1n1n1n1n1		 
//	c2		n1n1n1n1n1n1n1n1n1n1n1		c2
//	c2									c2
//	c2c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c2
//	c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2
//	c2c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c2
//	c2									c2
//	c2		n2n2n2n2n2n2n2n2n2n2n2		c2
//	 		n2n2n2n2n2n2n2n2n2n2n2		 
//
/***********************************************************************/
#define CH_E_W 2
double avg_t_test_neck(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double error_th, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e)
{
	int i,ii,num_out = 0,flag;
	int num_ch = 0, num_nch = 0;
	int num_1ch = 0, num_2ch = 0, num_3ch = 0;
	int num_n1ch = 0, num_n2ch = 0;
	int num_s = 0;
	double di, dj, w, l;
	double w_5, l_5, tau = gaussian_tau;
	double dtmp, dtmp2;
	double sum_all = 0, sum_all2 = 0;
	double sum_ch = 0, sum_ch2 = 0;
	double sum_1ch = 0, sum_1ch2 = 0;
	double sum_2ch = 0, sum_2ch2 = 0;
	double sum_3ch = 0, sum_3ch2 = 0;
	double sum_nch = 0, sum_nch2 = 0;
	double sum_n1ch = 0, sum_n1ch2 = 0;
	double sum_n2ch = 0, sum_n2ch2 = 0;
	double sum_s = 0;
	double cx, cy, likely;
	double t_test, t_test1;//, t_test2; 
	double sin_t = sin(t), cos_t = cos(t);
	double di_sin_t, di_cos_t, dj_sin_t, dj_cos_t;
	double dj_cost[100], dj_sint[100], di_cost[100], di_sint[100];
	DPoint pt;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	w_5 = w/2.;
	l_5 = l/2.;
	cx = (pt1.x+pt4.x)/2.;
	cy = (pt1.y+pt4.y)/2.;

	(*num) = 0;

	for(dj = -l_5, i = 0; dj<=l_5+1; dj += STEP_DX, i++){
		dj_cost[i] = dj*cos_t;
		dj_sint[i] = dj*sin_t;
	}
	for(di = -w_5, i = 0; di<=w_5+1; di += STEP_DY, i++){
		di_cost[i] = di*cos_t;
		di_sint[i] = di*sin_t;
	}
	#define END_WIDTH 2//4
	// non channel 1
	for(di = -w_5;di<= -w_5+STEP_DY;di+=STEP_DY){
		di_sin_t = di*sin_t;
		di_cos_t = di*cos_t;
		for(dj = -l_5+END_WIDTH*STEP_DX, i = END_WIDTH; dj<=l_5-END_WIDTH*STEP_DX; dj += STEP_DX, i++){
			flag = 0;
			// non channel 1
			pt.x = cx + dj_cost[i] + di_sin_t;
			pt.y = cy + dj_sint[i] - di_cos_t;
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_n1ch += dtmp;
				sum_n1ch2 += dtmp*dtmp;
				num_n1ch++;
				(*num)++;
				flag = 1;
			}
			else
				num_out++;
			// non channel 2
			pt.x = cx + dj_cost[i] - di_sin_t;
			pt.y = cy + dj_sint[i] + di_cos_t;
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp2 = real_coord(y,pt.y,pt.x);
				sum_n2ch += dtmp2;
				sum_n2ch2 += dtmp2*dtmp2;
				num_n2ch++;
				(*num)++;
				if(flag){
					sum_s += fabs(dtmp-dtmp2);
					num_s++;
				}
			}
			else
				num_out++;
		}
	}
	// channel 1 and 3
	di = -STEP_DY;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5+STEP_DX, i = 1; dj<=l_5-STEP_DX; dj += STEP_DX, i++){
		flag = 0;
		// channel 1
		pt.x = cx + dj_cost[i] + di_sin_t;
		pt.y = cy + dj_sint[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_1ch += dtmp;
			sum_1ch2 += dtmp*dtmp;
			num_1ch++;
			(*num)++;
			flag = 1;
		}
		else
			num_out++;
		// channel 3
		pt.x = cx + dj_cost[i] - di_sin_t;
		pt.y = cy + dj_sint[i] + di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp2 = real_coord(y,pt.y,pt.x);
			sum_3ch += dtmp2;
			sum_3ch2 += dtmp2*dtmp2;
			num_3ch++;
			(*num)++;
			if(flag){
				sum_s += fabs(dtmp-dtmp2);
				num_s++;
			}
		}
		else
			num_out++;
	}

	for(ii = 0; ii<2; ii++){
		if (ii==0)
			dj = -l_5+STEP_DX; // left end
		else
			dj = l_5-STEP_DX; // right end
		dj_sin_t = dj*sin_t;
		dj_cos_t = dj*cos_t;
		for(di = -w_5, i = 0; di<-STEP_DY; di += STEP_DY, i++){
			flag = 0;
			// channel 1
			pt.x = cx + dj_cos_t + di_sint[i];
			pt.y = cy + dj_sin_t - di_cost[i];
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_1ch += dtmp;
				sum_1ch2 += dtmp*dtmp;
				num_1ch++;
				(*num)++;
				flag = 1;
			}
			else
				num_out++;
			// channel 3
			pt.x = cx + dj_cos_t - di_sint[i];
			pt.y = cy + dj_sin_t + di_cost[i];
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp2 = real_coord(y,pt.y,pt.x);
				sum_3ch += dtmp2;
				sum_3ch2 += dtmp2*dtmp2;
				num_3ch++;
				(*num)++;
				if(flag){
					sum_s += fabs(dtmp-dtmp2);
					num_s++;
				}
			}
			else
				num_out++;
		}
	}
	// channel 2
	di = 0;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5+STEP_DX, i = 1; dj<=l_5-STEP_DX; dj += STEP_DX, i++){
		pt.x = cx + dj_cost[i] + di_sin_t;
		pt.y = cy + dj_sint[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_2ch += dtmp;
			sum_2ch2 += dtmp*dtmp;
			num_2ch++;
			(*num)++;
		}
		else
			num_out++;
	}
	for(ii = 0; ii<2; ii++){
		if (ii==0)
			dj = -l_5; // left end
		else
			dj = l_5; // right end
		dj_sin_t = dj*sin_t;
		dj_cos_t = dj*cos_t;
		for(di = -w_5, i = 0; di<=w_5; di += STEP_DY, i++){
			if(fabs(di)<=CH_E_W){
				flag = 0;
				pt.x = cx + dj_cos_t + di_sint[i];
				pt.y = cy + dj_sin_t - di_cost[i];
				if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
					dtmp = real_coord(y,pt.y,pt.x);
					sum_2ch += dtmp;
					sum_2ch2 += dtmp*dtmp;
					num_2ch++;
					(*num)++;
					flag = 1;
				}
				else
					num_out++;
				pt.x = cx + dj_cos_t - di_sint[i];
				pt.y = cy + dj_sin_t + di_cost[i];
				if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pt.y,pt.x);
					sum_2ch += dtmp2;
					sum_2ch2 += dtmp2*dtmp2;
					num_2ch++;
					(*num)++;
					if(flag){
						sum_s += fabs(dtmp-dtmp2);
						num_s++;
					}
				}
				else
					num_out++;
			}
		}
	}
	// Calculate distance1
	if(((double)(*num))/(double)(num_out+*num)<OBJECT_IN_RATIO){
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
	num_ch += num_1ch;
	sum_ch += sum_1ch;
	sum_ch2 += sum_1ch2;
	num_ch += num_2ch;
	sum_ch += sum_2ch;
	sum_ch2 += sum_2ch2;
	num_ch += num_3ch;
	sum_ch += sum_3ch;
	sum_ch2 += sum_3ch2;
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
		sum_s = sum_s/(double)num_s; // average of fabs(left hand - right hand)
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
	if(*num){
		sum_all = sum_all/(double)(*num); // average of all index pixel
		sum_all2 = sum_all2/(double)(*num);
	}
	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel
	sum_all2 = sum_all2 - sum_all*sum_all; // variance of the outside of a chanel
	if(*num !=0){
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		t_test1 = fabs(sum_ch-sum_nch)/(sqrt(sum_ch2/(double)num_ch+sum_nch2/(double)num_nch));
		t_test = t_test1/fmax(SYM_TH,sum_s/sum_all)*SYM_TH;//*sqrt(sum_all2)
		if ((sum_nch>sum_ch)&&(sum_all2>150))
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
			likely = pow(likely,1.0/3.0);
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/3.0);

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
	e[21] = 0;
	e[22] = 0;
	e[23] = 0;
	e[24] = sum_all2;
	e[25] = t_test;
	e[26] = sum_s;
	e[27] = 0;

	return likely;
}

/***********************************************************************/
//
//	n1n1n1n1n1n1n1n1n1n1n1n1n1n1n1		
//	n1n1n1n1n1n1n1n1n1n1n1n1n1n1n1		c2
//	n3n3								c2
//	n3n3    c1c1c1c1c1c1c1c1c1c1c1c1c1c1c2
//	n3n3    c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2
//	n3n3    c3c3c3c3c3c3c3c3c3c3c3c3c3c3c2
//	n3n3								c2
//	n2n2n2n2n2n2n2n2n2n2n2n2n2n2n2		c2
//	n2n2n2n2n2n2n2n2n2n2n2n2n2n2n2		
//
/***********************************************************************/
double avg_t_test_dent(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double error_th, double dent_l_w_ratio, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e)
{
	int i,ii,num_out,n1,descon_num, flag;
	int num_ch = 0, num_nch = 0;
	int num_1ch = 0, num_2ch = 0, num_3ch = 0;
	int num_n1ch = 0 , num_n2ch = 0, num_n3ch = 0;
	int num_s = 0;
	double di, dj, w, l;
	double w_5, l_5, tau = gaussian_tau;
	double dtmp, dtmp2, dtmp3, dtmp4;
	double sum_all = 0, sum_all2 = 0;
	double sum_ch = 0, sum_ch2 = 0;
	double sum_1ch = 0, sum_1ch2 = 0;
	double sum_2ch = 0, sum_2ch2 = 0;
	double sum_3ch = 0, sum_3ch2 = 0;
	double sum_nch = 0, sum_nch2 = 0;
	double sum_n1ch = 0, sum_n1ch2 = 0;
	double sum_n2ch = 0, sum_n2ch2 = 0;
	double sum_n3ch = 0, sum_n3ch2 = 0;
	double sum_s = 0;
	double sum_cr = 0;
	double cx, cy, likely;
	double t_test, t_test1, t_test_ch2, t_test_ch1; 
	double sin_t = sin(t), cos_t = cos(t);
	double di_sin_t, di_cos_t, dj_sin_t, dj_cos_t;
	double dj_cost[100], dj_sint[100], di_cost[100], di_sint[100];
	double data[400];
	double mean, std;
	DPoint pt;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	w_5 = w/2.;
	l_5 = l/2.;
	
	/*Find Center*/
	cx = (pt1.x+pt4.x)/2.;
	cy = (pt1.y+pt4.y)/2.;
	
	//	ch_end_width = w*CH_E_W; // max = w*0.5

	(*num) = 0;
	num_out = 0;

	for(dj = -l_5, i = 0; dj<=l_5+1; dj += STEP_DX, i++)
	{
		dj_cost[i] = dj*cos_t;
		dj_sint[i] = dj*sin_t;
	}
	for(di = -w_5, i = 0; di<=w_5+1; di += STEP_DY, i++)
	{
		di_cost[i] = di*cos_t;
		di_sint[i] = di*sin_t;
	}
	#define DENT_END_WIDTH 4 // 4 STEP_DX
	
	// non channel 1
	for(di = -w_5;di<=-w_5+STEP_DY;di += STEP_DY)
	{
		di_sin_t = di*sin_t;
		di_cos_t = di*cos_t;
		for(dj = -l_5, i = 0; dj<=l_5-END_WIDTH*STEP_DX; dj += STEP_DX, i++)
		{
			flag = 0;
			// non channel 1
			pt.x = cx + dj_cost[i] + di_sin_t;
			pt.y = cy + dj_sint[i] - di_cos_t;
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			{
				dtmp = real_coord(y,pt.y,pt.x);
				sum_n1ch += dtmp;
				sum_n1ch2 += dtmp*dtmp;
				num_n1ch++;
				(*num)++;
				flag = 1;
			}
			else
			{
				num_out++;
			}

			// non channel 2
			pt.x = cx + dj_cost[i] - di_sin_t;
			pt.y = cy + dj_sint[i] + di_cos_t;
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp2 = real_coord(y,pt.y,pt.x);
				sum_n2ch += dtmp2;
				sum_n2ch2 += dtmp2*dtmp2;
				num_n2ch++;
				(*num)++;
				if(flag){
					sum_s += fabs(dtmp-dtmp2);
					num_s++;
				}
			}
			else
				num_out++;
		}
	}
	// non channel 3
	for(ii = 0;ii<2;ii++){
		if(ii==0)
			dj = -l_5+STEP_DX;
		else
			dj = -l_5;
		dj_sin_t = dj*sin_t;
		dj_cos_t = dj*cos_t;
		for(di = -w_5+2*STEP_DY, i = 2; di<= 0; di += STEP_DY, i++){
			flag = 0;
			pt.x = cx + dj_cos_t + di_sint[i];
			pt.y = cy + dj_sin_t - di_cost[i];
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_n3ch += dtmp;
				sum_n3ch2 += dtmp*dtmp;
				num_n3ch++;
				(*num)++;
				flag = 1;
			}
			else
				num_out++;
			pt.x = cx + dj_cos_t - di_sint[i];
			pt.y = cy + dj_sin_t + di_cost[i];
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp2 = real_coord(y,pt.y,pt.x);
				sum_n3ch += dtmp2;
				sum_n3ch2 += dtmp2*dtmp2;
				num_n3ch++;
				(*num)++;
				if(flag){
					sum_s += fabs(dtmp-dtmp2);
					num_s++;
				}
			}
			else
				num_out++;
		}
	}
	// channel 1 and 3
	di = -STEP_DY;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5+END_WIDTH*STEP_DX, i = END_WIDTH; dj<=l_5-STEP_DX; dj += STEP_DX, i++){
		// channel 1
		pt.x = cx + dj_cost[i] + di_sin_t;
		pt.y = cy + dj_sint[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_1ch += dtmp;
			sum_1ch2 += dtmp*dtmp;
			num_1ch++;
			(*num)++;
			flag = 1;
		}
		else
			num_out++;
		// channel 3
		pt.x = cx + dj_cost[i] - di_sin_t;
		pt.y = cy + dj_sint[i] + di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp2 = real_coord(y,pt.y,pt.x);
			sum_3ch += dtmp2;
			sum_3ch2 += dtmp2*dtmp2;
			num_3ch++;
			(*num)++;
			if(flag){
				sum_s += fabs(dtmp-dtmp2);
				num_s++;
			}
		}
		else
			num_out++;
	}

	// channel 2
	di = 0;
	di_sin_t = di*sin_t;
	di_cos_t = di*cos_t;
	for(dj = -l_5+END_WIDTH*STEP_DX, i = END_WIDTH; dj<=l_5-STEP_DX; dj += STEP_DX, i++){
		pt.x = cx + dj_cost[i] + di_sin_t;
		pt.y = cy + dj_sint[i] - di_cos_t;
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			dtmp = real_coord(y,pt.y,pt.x);
			sum_2ch += dtmp;
			sum_2ch2 += dtmp*dtmp;
			num_2ch++;
			(*num)++;
		}
		else
			num_out++;
	}
	dj = l_5;
	dj_sin_t = dj*sin_t;
	dj_cos_t = dj*cos_t;
	for(di = -w_5, i = 0; di<=w_5; di += STEP_DY, i++){
		if(fabs(di)<=CH_E_W){
			flag = 0;
			pt.x = cx + dj_cos_t + di_sint[i];
			pt.y = cy + dj_sin_t - di_cost[i];
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_2ch += dtmp;
				sum_2ch2 += dtmp*dtmp;
				num_2ch++;
				(*num)++;
				flag = 1;
			}
			else
				num_out++;
			pt.x = cx + dj_cos_t - di_sint[i];
			pt.y = cy + dj_sin_t + di_cost[i];
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp2 = real_coord(y,pt.y,pt.x);
				sum_2ch += dtmp2;
				sum_2ch2 += dtmp2*dtmp2;
				num_2ch++;
				(*num)++;
				if(flag){
					sum_s += fabs(dtmp-dtmp2);
					num_s++;
				}
			}
			else
				num_out++;
		}
	}
	// calculate distance1
	if(((double)(*num))/(double)(num_out+*num)<OBJECT_IN_RATIO){
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
	
	num_nch += num_n3ch;
	sum_nch += sum_n3ch;
	sum_nch2 += sum_n3ch2;
	
	num_ch += num_1ch;
	sum_ch += sum_1ch;
	sum_ch2 += sum_1ch2;
	
	num_ch += num_2ch;
	sum_ch += sum_2ch;
	sum_ch2 += sum_2ch2;
	
	num_ch += num_3ch;
	sum_ch += sum_3ch;
	sum_ch2 += sum_3ch2;
	
	sum_all += sum_nch;
	sum_all2 += sum_nch2;
	
	sum_all += sum_ch;
	sum_all2 += sum_ch2;

	if(num_n1ch){
		sum_n1ch = sum_n1ch/(double)num_n1ch;
		sum_n1ch2 = sum_n1ch2/(double)num_n1ch;
	}
	if(num_n2ch){
		sum_n2ch = sum_n2ch/(double)num_n2ch;
		sum_n2ch2 = sum_n2ch2/(double)num_n2ch;
	}
	if(num_n3ch){
		sum_n3ch = sum_n3ch/(double)num_n3ch;
		sum_n3ch2 = sum_n3ch2/(double)num_n3ch;
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
	
	sum_1ch2 = sum_1ch2 - sum_1ch*sum_1ch; // variance of a channel 1
	sum_2ch2 = sum_2ch2 - sum_2ch*sum_2ch; // variance of a channel 2
	sum_3ch2 = sum_3ch2 - sum_3ch*sum_3ch; // variance of a channel 3
	
	sum_n1ch2 = sum_n1ch2 - sum_n1ch*sum_n1ch; // variance of the outside of a chanel
	sum_n2ch2 = sum_n2ch2 - sum_n2ch*sum_n2ch; // variance of the outside of a chanel
	sum_n3ch2 = sum_n3ch2 - sum_n3ch*sum_n3ch; // variance of the outside of a chanel

	if(num_ch){
		sum_ch = sum_ch/(double)num_ch;
		sum_ch2 = sum_ch2/(double)num_ch;
	}
	if(num_nch){
		sum_nch = sum_nch/(double)num_nch;
		sum_nch2 = sum_nch2/(double)num_nch;
	}
	if(*num){
		sum_all = sum_all/(double)(*num);
		sum_all2 = sum_all2/(double)(*num);
	}
	
	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel
	sum_all2 = sum_all2 - sum_all*sum_all; // variance of the outside of a chanel
	if(*num !=0){
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		t_test1 = fabs(sum_ch-sum_nch)/(sqrt(sum_ch2/(double)num_ch+sum_nch2/(double)num_nch));
		t_test = t_test1/fmax(SYM_TH,sum_s/sum_all)*SYM_TH;//*sqrt(sum_all2)
		
		t_test_ch1 = fabs(sum_1ch-sum_n1ch)/(sqrt(sum_1ch2/(double)num_1ch+sum_n1ch2/(double)num_n1ch));
		t_test_ch2 = fabs(sum_3ch-sum_n2ch)/(sqrt(sum_3ch2/(double)num_3ch+sum_n2ch2/(double)num_n2ch));

		int ratio_between_both_channels_t_test1, ratio_between_both_channels_t_test2;
		int ratio_boolean = 0;
		
		ratio_between_both_channels_t_test1 = t_test_ch1 < 1.2*t_test_ch2;
		ratio_between_both_channels_t_test2 = t_test_ch2 < 1.2*t_test_ch1;
		
		ratio_boolean = ratio_between_both_channels_t_test1 * ratio_between_both_channels_t_test2;
		
		if ((sum_nch>sum_ch)&&(sum_all2>150) &&  ratio_boolean)
		{
			if(t_test < 0.1 * error_th)
					likely = 1-t_test/error_th;
				else
					likely = exp(-(t_test-error_th)/(3.*t_test))-1;
		}
		else
			likely = 1;

		if (likely >= 0)
			//likely = pow(likely,1.0/4.0);
			likely = pow(likely,1.0/3.0);
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/3.0);

	}
	else{
		likely = INF;
	}

	// check discontinuity method 2
	n1 = 1;
	di = -w_5;			//Middle Point of Width
	data[0] = mean1; 	//Mean of class 1
	sum_cr = 0;

	for(dj = l_5-w_5; dj>-l_5+w_5; dj -= STEP_DX)
	{
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		
		sum_cr += data[n1];
		n1++;
	}
	
	//Curvature
	dtmp2 = M_PI/2.;
	dtmp3 = M_PI*3./2.;
	dtmp4 = M_PI/(2.*w);
	for(dtmp = dtmp2; dtmp<dtmp3; dtmp += dtmp4){
		di = -w_5*sin(dtmp);
		dj = w_5*cos(dtmp)+w_5-l_5;
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);

		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	
	
	di = w_5;
	
	for(dj = -l_5+w_5; dj<=l_5-w_5; dj += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	
	mean = sum_cr/(double)(n1-1);
	std = 0;
	for(i = 1; i< n1; i++){
		std += (data[i]-mean)*(data[i]-mean);
	}
	std = sqrt(std/(double)(n1-1));

	for(i = 1; i< n1; i++){
		data[i] = (data[i]-mean)/std; //it was /mean; -> Create a normal distribution N(0,1)
	}
	descon_num = 0;
	for(i = 3; i< n1-2; i++){
		//dtmp = (data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2])/5.;
		dtmp = data[i];
		dtmp = dtmp < 0? -dtmp: dtmp;
		if(dtmp > 2)
		{
			
			//printf("Discontinuity Here, dtmp: %f \n", dtmp);
			descon_num++;
		}
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
	e[21] = (double)num_n3ch;
	e[22] = sum_n3ch;
	e[23] = sum_n3ch2;
	e[24] = sum_all2;
	e[25] = t_test;
	e[26] = sum_s;
	e[27] = (double)descon_num;

	return likely;
}


void draw_nds_neck_2D(NeckDent *mp, MPP_Parameters mpp, unsigned char **image, int cols, int rows)
{
	unsigned char num;
	double amp, moffset;
	double di, dj, adi, adj, w, l;
	double w_5, l_5, mt_ww, tau = mpp.gaussian_tau;
	double y2, dtmp;
//	double min = 255, max = 0;
	DPoint pt;
	
	l = mp->length;
	w = mp->width;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);
	amp = mp->e[0];
	moffset = mp->e[2];

	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			adj = fabs(dj);
			if(adj<l_5-w_5){
				y2 = di*di*mt_ww;
			}
			else{
				dtmp = distance1(l_5-w_5,w_5,adj,adi);
				if(dtmp<w_5){
					y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
				}
				else{//4
					y2 = 0;
				}
			}
			y2 = -exp(y2);
//			if(y2<min) min = y2;
//			if(y2>max) max = y2;
			//y2 = -exp(y2);
			pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
			pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				if(y2<GAUSS_TH)
					num = 0;
				else 
					num = 255;
				image[(int)pt.y][(int)pt.x] = num;
			}
		}
	}
//	printf("min = %1.4f, max = %1.4f, diff = %1.4f, amp = %1.4f\n",min, max, max-min, mp->e0);
}



void draw_nds_dent_2D(NeckDent *mp, MPP_Parameters mpp, unsigned char **image, int cols, int rows)
{
	unsigned char num;
	double amp, moffset;
	double di, dj, adi, w, l;
	double w_5, l_5, mt_ww, tau = mpp.gaussian_tau;
	double y2, dtmp;
	DPoint pt;
	
	l = mp->length;
	w = mp->width;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);
	amp = mp->e[0];
	moffset = mp->e[2];

	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			if(dj>0){
				if(dj<l_5-w_5){
					y2 = di*di*mt_ww;
				}
				else{
					dtmp = distance1(l_5-w_5,w_5,dj,adi);
					if(dtmp<w_5){
						y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
					}
					else{//4
						y2 = 0;
					}
				}
			}
			else{
				if(dj>-l_5+w_5){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
				}
			}
			y2 = -exp(y2);
			//y2 = -exp(y2);
			pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
			pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				if(y2<GAUSS_TH)
					num = 255-0;
				else 
					num = 255-255;
				image[(int)pt.y][(int)pt.x] = num;
			}
		}
	}
}



void draw_all_nds_2D(NeckDent *mp, MPP_Parameters mpp, int np_num, 
					 unsigned char **image, unsigned char **image2, int cols, int rows)
{
	int i,j,k;

	for (i = 0;i<rows;i++){
		for (j = 0;j<cols;j++){
			image2[i][j] = 0;//image[i][j];
		}
	}

	for(k = 0; k<np_num; k++){
		if(mp[k].type == NDTYPE_NECKING)
			draw_nds_neck_2D(&(mp[k]), mpp, image2, cols, rows);
		else
			draw_nds_dent_2D(&(mp[k]), mpp, image2, cols, rows);

	}
}




double avg_Bhattacharya(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double error_th, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e)
{
	int num2, n1;
	double di, dj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2;
	double sum_s, sum_all, sum_ch, sum_nch, sum_all2, sum_ch2, sum_nch2, error_s, likely;
	double Bhatta, Bhatta1, Bhatta2 = 0; 
	int num_ch, num_nch;
	DPoint pt,pts;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_ch = 0;
	sum_nch = 0;
	sum_ch2 = 0;
	sum_nch2 = 0;
	num_ch = 0;
	num_nch = 0;
	sum_s = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			y2 = (di)*(di)*mt_ww;
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di for symmetry potential
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di for symmetry potential
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				if(y2<GAUSS_TH){
					sum_ch += dtmp;
					sum_ch2 += dtmp*dtmp;
					num_ch++;
				}
				else{
					sum_nch += dtmp;
					sum_nch2 += dtmp*dtmp;
					num_nch++;
				}
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){
		e[0] = 0;
		e[1] = 0;
		e[2] = 0;
		e[3] = 0;
		e[4] = 0;
		e[5] = 0;
		return INF;
	}
	sum_all = (sum_ch+sum_nch)/(double)num2;
	sum_all2 = (sum_ch2+sum_nch2)/(double)num2;
	sum_all2 = sum_all2 - sum_all*sum_all;
	sum_ch = sum_ch/(double)num_ch; // mu of a channel
	sum_ch2 = sum_ch2/(double)num_ch; 
	sum_nch = sum_nch/(double)num_nch; // mu of the outside of a channel
	sum_nch2 = sum_nch2/(double)num_nch; 
	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel

	if(*num !=0){
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		Bhatta1 = pow(fabs(sum_ch-sum_nch),1)/(sqrt(sum_ch2+sum_nch2));
		Bhatta2 = 1*log(2.*sqrt(sum_ch2*sum_nch2)/(sum_ch2+sum_nch2));
		Bhatta = Bhatta1 + Bhatta2;
		error_s = sum_s/((double)(n1));

		if (sum_nch>sum_ch)
		{
			if(Bhatta < error_th)
					likely = 1-Bhatta/error_th;
				else
					likely = exp(-(Bhatta-error_th)/(3.*Bhatta))-1;
		}
		else
			likely = 1;

		if (likely >= 0)
			//likely = pow(likely,1.0/4.0);
			likely = pow(likely,1.0/3.0);
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/3.0);

	}
	else{
		likely = INF;
	}
	error_s = sqrt(error_s);
	e[0] = sum_ch;
	e[1] = sum_nch;
	e[2] = sum_ch2;
	e[3] = sum_nch2;
	e[4] = sum_all2;
	e[5] = 0;

	return likely;
}

double avg_Bhattacharya_neck(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double error_th, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e)
{
	int num2, n1;
	double di, dj, adi, adj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2;
	double sum_s, sum_all, sum_ch, sum_nch, sum_all2, sum_ch2, sum_nch2, error_s, likely;
	double Bhatta, Bhatta1, Bhatta2 = 0; 
	int num_ch, num_nch;
	DPoint pt,pts;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_ch = 0;
	sum_nch = 0;
	sum_ch2 = 0;
	sum_nch2 = 0;
	num_ch = 0;
	num_nch = 0;
	sum_s = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			adj = fabs(dj);
			if(adj<l_5-w_5){
				y2 = di*di*mt_ww;
			}
			else{
				dtmp = distance1(l_5-w_5,w_5,adj,adi);
				if(dtmp<w_5){
					y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
				}
				else{//4
					y2 = 0;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di for symmetry potential
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di for symmetry potential
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				if(y2<GAUSS_TH){
					sum_ch += dtmp;
					sum_ch2 += dtmp*dtmp;
					num_ch++;
				}
				else{
					sum_nch += dtmp;
					sum_nch2 += dtmp*dtmp;
					num_nch++;
				}
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		e[0] = 0;
		e[1] = 0;
		e[2] = 0;
		e[3] = 0;
		e[4] = 0;
		e[5] = 0;
		return INF;
	}
	sum_all = (sum_ch+sum_nch)/(double)num2;
	sum_all2 = (sum_ch2+sum_nch2)/(double)num2;
	sum_all2 = sum_all2 - sum_all*sum_all;
	sum_ch = sum_ch/(double)num_ch; // mu of a channel
	sum_ch2 = sum_ch2/(double)num_ch; 
	sum_nch = sum_nch/(double)num_nch; // mu of the outside of a channel
	sum_nch2 = sum_nch2/(double)num_nch; 
	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel

	if(*num !=0){
		//Bhatta = (sum_ch-sum_nch)*(sum_ch-sum_nch)/(sqrt(sum_ch2+sum_nch2));
		Bhatta1 = pow(fabs(sum_ch-sum_nch),1)/(sqrt(sum_ch2+sum_nch2));
		Bhatta2 = 1*log(2.*sqrt(sum_ch2*sum_nch2)/(sum_ch2+sum_nch2));
		Bhatta = Bhatta1 + Bhatta2;
		error_s = sum_s/((double)(n1));

		if (sum_nch>sum_ch)
		{
			if(Bhatta < error_th)
					likely = 1-Bhatta/error_th;
				else
					likely = exp(-(Bhatta-error_th)/(3.*Bhatta))-1;
		}
		else
			likely = 1;

		if (likely >= 0)
			//likely = pow(likely,1.0/4.0);
			likely = pow(likely,1.0/4.0);
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/4.0);

	}
	else{
		likely = INF;
	}
	error_s = sqrt(error_s);
	e[0] = sum_ch;
	e[1] = sum_nch;
	e[2] = sum_ch2;
	e[3] = sum_nch2;
	e[4] = sum_all2;
	e[5] = sum_all2;

	return likely;
}

double avg_Bhattacharya_dent(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double error_th, double dent_l_w_ratio, double t,
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e)
{
	int i, num2, n1, descon_num;
	double di, dj, adi, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2, dtmp3, dtmp4;
	double sum_s, sum_all, sum_ch, sum_nch, sum_all2, sum_ch2, sum_nch2, error_s, likely;
	double Bhatta, Bhatta1, Bhatta2 = 0; 
	int num_ch, num_nch;
	DPoint pt,pts;
	double sum_cr;
	double data[400];
	double mean, std;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	if (dent_l_w_ratio*w>l){
		e[0] = 0;
		e[1] = 0;
		e[2] = 0;
		e[3] = 0;
		e[4] = 0;
		e[5] = 0;
		return INF;
	}
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_ch = 0;
	sum_nch = 0;
	sum_ch2 = 0;
	sum_nch2 = 0;
	num_ch = 0;
	num_nch = 0;
	sum_s = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			if(dj>0){
				if(dj<l_5-w_5){
					y2 = di*di*mt_ww;
				}
				else{
					dtmp = distance1(l_5-w_5,w_5,dj,adi);
					if(dtmp<w_5){
						y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
					}
					else{//4
						y2 = 0;
					}
				}
			}
			else{
				if(dj>-l_5+w_5){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				if(y2<GAUSS_TH){
					sum_ch += dtmp;
					sum_ch2 += dtmp*dtmp;
					num_ch++;
				}
				else{
					sum_nch += dtmp;
					sum_nch2 += dtmp*dtmp;
					num_nch++;
				}
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		e[0] = 0;
		e[1] = 0;
		e[2] = 0;
		e[3] = 0;
		e[4] = 0;
		e[5] = 0;
		return INF;
	}
	sum_all = (sum_ch+sum_nch)/(double)num2;
	sum_all2 = (sum_ch2+sum_nch2)/(double)num2;
	sum_all2 = sum_all2 - sum_all*sum_all;
	sum_ch = sum_ch/(double)num_ch; // mu of a channel
	sum_ch2 = sum_ch2/(double)num_ch; 
	sum_nch = sum_nch/(double)num_nch; // mu of the outside of a channel
	sum_nch2 = sum_nch2/(double)num_nch; 
	sum_ch2 = sum_ch2 - sum_ch*sum_ch; // variance of a channel
	sum_nch2 = sum_nch2 - sum_nch*sum_nch; // variance of the outside of a chanel

	if(*num !=0){
		Bhatta1 = pow(fabs(sum_ch-sum_nch),1)/(sqrt(sum_ch2+sum_nch2));
		Bhatta2 = 1*log(2.*sqrt(sum_ch2*sum_nch2)/(sum_ch2+sum_nch2));
		Bhatta = Bhatta1 + Bhatta2;

		error_s = sum_s/((double)(n1));

		if (sum_nch>sum_ch)
		{
			if(Bhatta < error_th)
					likely = 1-Bhatta/error_th;
				else
					likely = exp(-(Bhatta-error_th)/(3.*Bhatta))-1;
		}
		else
			likely = 1;

		if (likely >= 0)
			//likely = pow(likely,1.0/4.0);
			likely = pow(likely,1.0/4.0);
		else
			//likely = -1*pow(-1*likely,1.0/4.0);
			likely = -1*pow(-1*likely,1.0/4.0);

	}
	else{
		likely = INF;
	}
	error_s = sqrt(error_s);

	// check discontinuity method 2
	n1 = 1;
	di = -w_5;
	data[0] = mean1;
	sum_cr = 0;
	for(dj = l_5-w_5; dj>-l_5+w_5; dj -= STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	dtmp2 = M_PI/2.;
	dtmp3 = M_PI*3./2.;
	dtmp4 = M_PI/(2.*w);
	for(dtmp = dtmp2; dtmp<dtmp3; dtmp += dtmp4){
		di = -w_5*sin(dtmp);
		dj = w_5*cos(dtmp)+w_5-l_5;
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	di = w_5;
	for(dj = -l_5+w_5; dj<=l_5-w_5; dj += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	mean = sum_cr/(double)(n1-1);
	std = 0;
	for(i = 1; i< n1; i++){
		std += (data[i]-mean)*(data[i]-mean);
	}
	std = sqrt(std/(double)(n1-1));

	for(i = 1; i< n1; i++){
		data[i] = (data[i]-mean)/mean;
	}
	descon_num = 0;
	for(i = 3; i< n1-2; i++){
		dtmp = (data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2])/5.;
		if(dtmp<DENT_DISCON_TH)
			descon_num++;
	}
	e[0] = sum_ch;
	e[1] = sum_nch;
	e[2] = sum_ch2;
	e[3] = sum_nch2;
	e[4] = sum_all2;
	e[5] = (double)descon_num;

	return likely;
}



double avg_Gaussian_neck(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3, double *e4)
{
	int i, num2, n1, descon_num;
	double amp, moffset;
	double di, dj, adi, adj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2, weight, weight_sum;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	double sum_cr;
	double data[400];
	double mean, std;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	weight_sum = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			adj = fabs(dj);
			if(adj<l_5-w_5){
				y2 = di*di*mt_ww;
				weight = y2;
			}
			else{
				dtmp = distance1(l_5-w_5,w_5,adj,adi);
				if(dtmp<w_5){
					y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
				}
				else{//4
					y2 = 0;
				}
				weight = ((l_5-w_5-adj)*(l_5-w_5-adj)+adi*adi)*mt_ww;
			}
			y2 = (-exp(y2));
			weight = WEIGHT_AMP*exp(weight)+WEIGHT_OFFSET;
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += weight*y2;
				sum_fy += weight*y2*dtmp;
				sum_y += weight*dtmp;
				sum_f2 += weight*y2*y2;
				sum_y2 += weight*dtmp*dtmp;
				(*num)++;
				weight_sum += weight;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		*e4 = 0;
		return INF;
	}
	amp = (weight_sum*sum_fy-sum_y*sum_f)/(weight_sum*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/(weight_sum);
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + moffset*moffset*weight_sum;
	if(*num !=0){
		error2 = error2/(weight_sum);
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;
	*e3 = error_s;

		// check discontinuity method 2
	n1 = 1;
	data[0] = mean1;
	sum_cr = 0;
	dj = -l_5;
	for(di = -w_5; di<=w_5; di += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	dj = l_5;
	for(di = -w_5; di<=w_5; di += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	mean = sum_cr/(double)(n1-1);
	std = 0;
	for(i = 1; i< n1; i++){
		std += (data[i]-mean)*(data[i]-mean);
	}
	std = sqrt(std/(double)(n1-1));

	for(i = 1; i< n1; i++){
		data[i] = (data[i]-mean)/mean;
	}
	descon_num = 0;
	for(i = 3; i< n1-2; i++){
		dtmp = (data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2])/5.;
		if(dtmp<NECK_DISCON_TH)
			descon_num++;
	}
	*e4 = *e3;// tmp  (double)descon_num;

	return error2;
}

double avg_Gaussian_dent(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double dent_l_w_ratio, double t,
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3, double *e4)
{
	int i, num2, n1, descon_num;
	double amp, moffset;
	double di, dj, adi, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2, dtmp3, dtmp4, weight, weight_sum;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	double sum_cr;
	double data[400];
	double mean, std;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	if (dent_l_w_ratio*w>l){
		(*e0) = 0;
		(*e1) = 0;
		(*e2) = 0;
		(*e3) = 0;
		(*e4) = 0;
		return INF;
	}
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	weight_sum = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			if(dj>0){
				if(dj<l_5-w_5){
					y2 = di*di*mt_ww;
					weight = y2;
				}
				else{
					dtmp = distance1(l_5-w_5,w_5,dj,adi);
					if(dtmp<w_5){
						y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
					}
					else{//4
						y2 = 0;
					}
					weight = ((l_5-w_5-dj)*(l_5-w_5-dj)+adi*adi)*mt_ww;
				}
			}
			else{
				if(dj>-l_5+w_5){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
				}
				weight = y2;
			}
			y2 = (-exp(y2));
			weight = WEIGHT_AMP*exp(weight)+WEIGHT_OFFSET;
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += weight*y2;
				sum_fy += weight*y2*dtmp;
				sum_y += weight*dtmp;
				sum_f2 += weight*y2*y2;
				sum_y2 += weight*dtmp*dtmp;
				(*num)++;
				weight_sum += weight;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		*e4 = 0;
		return INF;
	}
	amp = (weight_sum*sum_fy-sum_y*sum_f)/(weight_sum*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/(weight_sum);
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + moffset*moffset*weight_sum;
	if(*num !=0){
		error2 = error2/(weight_sum);
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;
	*e3 = error_s;

#if 0	
// check discontinuity method 1
	n1 = 0;
	sum_cr = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<-l_5+w_5; dj += STEP_DX){
			y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
			y2 = amp*(-exp(y2))+moffset;
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_cr += (y2-dtmp)*(y2-dtmp);
				(n1)++;
			}
		}
	}
	if(n1 !=0){
		sum_cr = sum_cr/((double)(n1));
	}
	else{
		sum_cr = INF;
	}
	sum_cr = sqrt(sum_cr);
	*e4 = sum_cr;

#else
	// check discontinuity method 2
	n1 = 1;
	di = -w_5;
	data[0] = mean1;
	sum_cr = 0;
	for(dj = l_5-w_5; dj>-l_5+w_5; dj -= STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	dtmp2 = M_PI/2.;
	dtmp3 = M_PI*3./2.;
	dtmp4 = M_PI/(2.*w);
	for(dtmp = dtmp2; dtmp<dtmp3; dtmp += dtmp4){
		di = -w_5*sin(dtmp);
		dj = w_5*cos(dtmp)+w_5-l_5;
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	di = w_5;
	for(dj = -l_5+w_5; dj<=l_5-w_5; dj += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	mean = sum_cr/(double)(n1-1);
	std = 0;
	for(i = 1; i< n1; i++){
		std += (data[i]-mean)*(data[i]-mean);
	}
	std = sqrt(std/(double)(n1-1));

	for(i = 1; i< n1; i++){
		data[i] = (data[i]-mean)/mean;
	}
	descon_num = 0;
	for(i = 3; i< n1-2; i++){
		dtmp = (data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2])/5.;
		if(dtmp<DENT_DISCON_TH)
			descon_num++;
	}
	*e4 = (double)descon_num;
#endif

	return error2;
}

// new fixed mean map (width : -w/2 ~ +w/2, length : -l/2 ~ +l/2)
double avg_Gaussian_neck2(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3)
{
	int num2, n1;
	double amp, moffset;
	double di, dj, adi, adj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			adj = fabs(dj);
			if(adj<l_5-w_5){
				y2 = di*di*mt_ww;
			}
			else{
				dtmp = distance1(l_5-w_5,w_5,adj,adi);
				if(dtmp<w_5){
					y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
				}
				else{//4
					y2 = 0;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp*dtmp;
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;
	*e3 = error_s;

	return error2;
}

double avg_Gaussian_dent2(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double dent_l_w_ratio, double t,
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3)
{
	int i, num2, n1, descon_num;
	double amp, moffset;
	double di, dj, adi, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2, dtmp3, dtmp4;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	double sum_cr;
	double data[400];
	double mean, std;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	if (dent_l_w_ratio*w>l){
		(*e0) = 0;
		(*e1) = 0;
		return INF;
	}
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			adi = fabs(di);
			if(dj>0){
				if(dj<l_5-w_5){
					y2 = di*di*mt_ww;
				}
				else{
					dtmp = distance1(l_5-w_5,w_5,dj,adi);
					if(dtmp<w_5){
						y2 = (w_5-dtmp)*(w_5-dtmp)*mt_ww;
					}
					else{//4
						y2 = 0;
					}
				}
			}
			else{
				if(dj>-l_5+w_5){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp*dtmp;
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;

#if 0
// check discontinuity method 1
	n1 = 0;
	sum_cr = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<-l_5+w_5; dj += STEP_DX){
			y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
			y2 = amp*(-exp(y2))+moffset;
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_cr += (y2-dtmp)*(y2-dtmp);
				(n1)++;
			}
		}
	}
	if(*num !=0){
		sum_cr = sum_cr/((double)(n1));
	}
	else{
		sum_cr = INF;
	}
	sum_cr = sqrt(sum_cr);
	*e3 = sum_cr;
#else
// check discontinuity method 2
	n1 = 1;
	di = -w_5;
	data[0] = mean1;
	sum_cr = 0;
	for(dj = l_5-w_5; dj>-l_5+w_5; dj -= STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	dtmp2 = M_PI/2.;
	dtmp3 = M_PI*3./2.;
	dtmp4 = M_PI/(2.*w);
	for(dtmp = dtmp2; dtmp<dtmp3; dtmp += dtmp4){
		di = -w_5*sin(dtmp);
		dj = w_5*cos(dtmp)+w_5-l_5;
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	di = w_5;
	for(dj = -l_5+w_5; dj<=l_5-w_5; dj += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		sum_cr += data[n1];
		n1++;
	}
	mean = sum_cr/(double)(n1-1);
	std = 0;
	for(i = 1; i< n1; i++){
		std += (data[i]-mean)*(data[i]-mean);
	}
	std = sqrt(std/(double)(n1-1));

	for(i = 1; i< n1; i++){
		data[i] = (data[i]-mean)/mean;
	}
	descon_num = 0;
	for(i = 3; i< n1-2; i++){
		dtmp = (data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2])/5.;
		if(dtmp<DENT_DISCON_TH)
			descon_num++;
	}
	*e3 = (double)descon_num;
#endif

	return error2;
}

double avg_Gaussian_neck3(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double t,  
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3)
{
	int num2, n1;
	double amp, moffset;
	double di, dj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			if(dj>0){
				if(l_5-dj>fabs(di)){//2
					y2 = (di)*(di)*mt_ww;
				}
				else{//4
					y2 = (l_5-dj)*(l_5-dj)*mt_ww;
				}
			}
			else{
				if(l_5+dj>fabs(di)){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = (l_5+dj)*(l_5+dj)*mt_ww;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp*dtmp;
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;
	*e3 = error_s;

	return error2;
}

#define WEIGHT_CRITICAL_REGION 0.1
// new fixed mean map (width : -w/2 ~ +w/2, length : -l/2 ~ +l/2)
double avg_Gaussian_dent3(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double dent_l_w_ratio, double t,
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3)
{
	int num2, n1;
	double amp, moffset;
	double di, dj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2, dtmp3, dtmp4;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	double sum_cr;
	double data[400];
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	if (dent_l_w_ratio*w>l){
		(*e0) = 0;
		(*e1) = 0;
		return INF;
	}
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			if(dj>0){
				if(l_5-dj>fabs(di)){//2
					y2 = (di)*(di)*mt_ww;
				}
				else{//4
					y2 = (l_5-dj)*(l_5-dj)*mt_ww;
				}
			}
			else{
				if(dj>-l_5+w_5){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += y2;
				sum_fy += y2*dtmp;
				sum_y += dtmp;
				sum_f2 += y2*y2;
				sum_y2 += dtmp*dtmp;
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		return INF;
	}
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;

	// check discontinuity method 1
	n1 = 0;
	sum_cr = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<-l_5+w_5; dj += STEP_DX){
			y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
			y2 = amp*(-exp(y2))+moffset;
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_cr += (y2-dtmp)*(y2-dtmp);
				(n1)++;
			}
		}
	}
	if(*num !=0){
		sum_cr = sum_cr/((double)(n1));
	}
	else{
		sum_cr = INF;
	}
	sum_cr = sqrt(sum_cr);
	*e3 = sum_cr;

	// check discontinuity method 2
	n1 = 1;
	di = -w_5;
	data[0] = mean1;
	for(dj = l_5-w_5; dj>-l_5+w_5; dj -= STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		n1++;
	}
	dtmp2 = M_PI/2.;
	dtmp3 = M_PI*3./2.;
	dtmp4 = M_PI/(2.*w);
	for(dtmp = dtmp2; dtmp<dtmp3; dtmp += dtmp4){
		di = -w_5*sin(dtmp);
		dj = w_5*cos(dtmp)+w_5-l_5;
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		n1++;
	}
	di = w_5;
	for(dj = -l_5+w_5; dj<=l_5-w_5; dj += STEP_DX){
		pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
		pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.))
			data[n1] = real_coord(y,pt.y,pt.x);
		else
			data[n1] = data[n1-1];
		n1++;
	}

	return error2;
}


void avg_Gaussian_dent_save(unsigned char **y, int cols, int rows, NeckDent *mp, MPP_Parameters mpp, int *num)
{
	int i, n1;
	double di, dj;
	double w_5, l_5;
	double dtmp, dtmp2, dtmp3, dtmp4;
	DPoint pt;
	double data[100],ptx[100],pty[100];
	FILE *fp1;
	char filename[400];

	sprintf(filename, "dent_continuity.txt");
	if ((fp1 = fopen(filename, "wb")) == NULL ) {
		printf("Cannot open file %s\n", filename);
		exit(1);
	}
	
	w_5 = mp->width/2.;
	l_5 = mp->length/2.;

	// check discontinuity method 2
	n1 = 1;
	di = -w_5;
	data[0] = mpp.mean[1];
	for(dj = l_5-w_5; dj>-l_5+w_5; dj -= STEP_DX){
		pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
		pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			data[n1] = real_coord(y,pt.y,pt.x);
			ptx[n1] = pt.x;
			pty[n1] = pt.y;
		}
		else{
			data[n1] = data[n1-1];
			ptx[n1] = pt.x;
			pty[n1] = pt.y;
		}
		n1++;
	}
	dtmp2 = M_PI/2.;
	dtmp3 = M_PI*3./2.;
	dtmp4 = M_PI/(2.*mp->width);
	for(dtmp = dtmp2; dtmp<dtmp3; dtmp += dtmp4){
		di = -w_5*sin(dtmp);
		dj = w_5*cos(dtmp)+w_5-l_5;
		pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
		pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			data[n1] = real_coord(y,pt.y,pt.x);
			ptx[n1] = pt.x;
			pty[n1] = pt.y;
		}
		else{
			data[n1] = data[n1-1];
			ptx[n1] = pt.x;
			pty[n1] = pt.y;
		}
		n1++;
	}
	di = w_5;
	for(dj = -l_5+w_5; dj<=l_5-w_5; dj += STEP_DX){
		pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
		pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
		if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
			data[n1] = real_coord(y,pt.y,pt.x);
			ptx[n1] = pt.x;
			pty[n1] = pt.y;
		}
		else{
			data[n1] = data[n1-1];
			ptx[n1] = pt.x;
			pty[n1] = pt.y;
		}
		n1++;
	}
	for(i = 1; i<n1; i++)
		fprintf(fp1,"%1.2f, ", data[i]);
	fprintf(fp1,"\r\n");
	for(i = 1; i<n1; i++)
		fprintf(fp1,"%1.2f, ", ptx[i]);
	fprintf(fp1,"\r\n");
	for(i = 1; i<n1; i++)
		fprintf(fp1,"%1.2f, ", pty[i]);
	fprintf(fp1,"\r\n");
	fprintf(fp1,"m1, %1.2f\r\n", mpp.mean[1]);
	fprintf(fp1,"v1, %1.2f\r\n", mpp.vari[1]);
	fclose(fp1);
}

// new fixed mean map (width : -w/2 ~ +w/2, length : -l/2 ~ +l/2)
double avg_Gaussian_dent4(unsigned char **y, int cols, int rows, double mean0, double vari0, double mean1, double vari1, 
					double variance, double gaussian_tau, double dent_l_w_ratio, double t,
					DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, int *num, double *e0, double *e1, double *e2, double *e3)
{
	int num2, n1;
	double amp, moffset;
	double di, dj, w, l;
	double w_5, l_5, mt_ww, tau = gaussian_tau;
	double y2, dtmp, dtmp2;
	double sum_f, sum_fy, sum_y, sum_f2, sum_y2, sum_s, error_s, error2;
	DPoint pt,pts;
	double weight, sum_weight;
	
	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	if (dent_l_w_ratio*w>l){
		(*e0) = 0;
		(*e1) = 0;
		return INF;
	}
	num2 = 0;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);

	(*num) = 0;
	n1 = 0;
	sum_weight = 0;
	sum_s = 0;
	sum_f = 0;
	sum_fy = 0;
	sum_y = 0;
	sum_f2 = 0;
	sum_y2 = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			weight = 1;
			if(dj>0){
				if(l_5-dj>fabs(di)){//2
					y2 = (di)*(di)*mt_ww;
				}
				else{//4
					y2 = (l_5-dj)*(l_5-dj)*mt_ww;
				}
			}
			else{
				if(dj>-l_5+w_5){//1
					y2 = (di)*(di)*mt_ww;
				}
				else{//3
					y2 = ((w_5-l_5-dj)*(w_5-l_5-dj)+di*di)*mt_ww;
					weight = WEIGHT_CRITICAL_REGION;
				}
			}
			y2 = (-exp(y2));
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			pts.x = (pt1.x+pt4.x)/2. + dj*cos(t) - di*sin(t);//-di
			pts.y = (pt1.y+pt4.y)/2. + dj*sin(t) + di*cos(t);//-di
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				dtmp = real_coord(y,pt.y,pt.x);
				sum_f += weight*y2;
				sum_fy += weight*y2*dtmp;
				sum_y += weight*dtmp;
				sum_f2 += weight*y2*y2;
				sum_y2 += weight*dtmp*dtmp;
				sum_weight += weight;
				(*num)++;
				if((pts.x >= 0.)&&(pts.x < (double)cols-1.)&&(pts.y >= 0.)&&(pts.y < (double)rows-1.)){
					dtmp2 = real_coord(y,pts.y,pts.x);
					sum_s += (dtmp-dtmp2)*(dtmp-dtmp2);
					n1++;
				}
			}
			num2++;
		}
	}
	if(((double)(*num))/(double)num2<OBJECT_IN_RATIO){ 
		*e0 = 0;
		*e1 = 0;
		*e2 = 0;
		*e3 = 0;
		return INF;
	}
	sum_f = sum_f/sum_weight*(double)(*num);
	sum_fy = sum_fy/sum_weight*(double)(*num);
	sum_y = sum_y/sum_weight*(double)(*num);
	sum_f2 = sum_f2/sum_weight*(double)(*num);
	sum_y2 = sum_y2/sum_weight*(double)(*num);
	amp = ((double)(*num)*sum_fy-sum_y*sum_f)/((double)(*num)*sum_f2-sum_f*sum_f);
	moffset = (sum_y-amp*sum_f)/((double)(*num));
	error2 = amp*amp*sum_f2 + sum_y2 - 2*amp*sum_fy +2*amp*moffset*sum_f - 2*moffset*sum_y + (double)(*num)*moffset*moffset;
	if(*num !=0){
		error2 = error2/((double)(*num));
		error_s = sum_s/((double)(n1));
	}
	else{
		error2 = INF;
	}
	error2 = sqrt(error2);
	error_s = sqrt(error_s);
	*e0 = amp;
	*e1 = error2;
	*e2 = moffset;
	*e3 = error_s;

	return error2;
}



int channel_shape(unsigned char **y, int cols, int rows, double t, DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, double *ch_y,  double *ch_x)
{
	double di, dj, w, w_5, l_5, sum0;
	DPoint pt;
	int i, n;

	l_5 = dist(pt2, pt1)/2.;
	w = dist(pt3, pt1);
	w_5 = w/2.;

	i = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		sum0 = 0.;
		n = 0;
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				sum0 += real_coord(y,pt.y,pt.x);
				n++;
			}
		}
		ch_y[i] = sum0/(double)n;	// mean
		ch_x[i] = di/w;
		i++;
	}
	i = 0;
	for(di = -w_5; di<=w_5; di += STEP_DY){
		sum0 = 0.;
		n = 0;
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			pt.x = (pt1.x+pt4.x)/2. + dj*cos(t) + di*sin(t);
			pt.y = (pt1.y+pt4.y)/2. + dj*sin(t) - di*cos(t);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				sum0 += (real_coord(y,pt.y,pt.x)-ch_y[i])*(real_coord(y,pt.y,pt.x)-ch_y[i]);
				n++;
			}
		}
//		ch_y[i] = sum0/(double)n;	// variance
		i++;
	}
	return i;
}

double min_val(double a, double b, double c, double d)
{
	double tmp;

	tmp = a;
	if (a>b){
		tmp = b;
		if (b>c){
			tmp = c;
			if (c>d) tmp = d;
		}
		else{ // c>=b
			if (b>d) tmp = d;
		}
	}
	else{ // b>=a
		if (a>c){
			tmp = c;
			if (c>d) tmp = d;
		}
		else{ // c>=a
			if (a>d) tmp = d;
		}
	}
	return tmp;
}

double max_val(double a, double b, double c, double d)
{
	double tmp;

	tmp = a;
	if (a<b){
		tmp = b;
		if (b<c){
			tmp = c;
			if (c<d) tmp = d;
		}
		else{ // c<=b
			if (b<d) tmp = d;
		}
	}
	else{ // b<=a
		if (a<c){
			tmp = c;
			if (c<d) tmp = d;
		}
		else{ // c<=a
			if (a<d) tmp = d;
		}
	}
	return tmp;
}

void Gaussian_beta(double **beta_img[], double *beta, int cols, int rows, NeckDent *mp, MPP_Parameters mpp)
{
	double di,dj;
	double w, l, w_5, l_5, w_n, mt_ww, tau = mpp.gaussian_tau;//20 
	double ma, y2;
	int i,j,minx, maxx, miny, maxy;
	DPoint pt, pt1, pt2, pt3, pt4, center;
	double t, dtmp;


	w = mp->width;// 5 mp->width;
	mp->length += 3;// 6 draw more zeros both side of length direction.
	

	l = mp->length;
	nd_GenCornerPoint(mp);
	pt1 = mp->r[0];
	pt2 = mp->r[1];
	pt3 = mp->r[2];
	pt4 = mp->r[3];
	t = mp->theta;
	center = mp->center;

	w_5 = w/2.;
	l_5 = l/2.; 
	w_n = w/4.;
	mt_ww = -tau/(w*w);
	ma = (beta[1]-beta[0])/(1-exp(-tau/4.));

	minx = (int)min_val(pt1.x, pt2.x, pt3.x, pt4.x);
	maxx = (int)max_val(pt1.x, pt2.x, pt3.x, pt4.x);
	miny = (int)min_val(pt1.y, pt2.y, pt3.y, pt4.y);
	maxy = (int)max_val(pt1.y, pt2.y, pt3.y, pt4.y);
	if(minx<0) minx = 0;
	if(miny<0) miny = 0;
	if(maxx>(double)cols-1) maxx = cols-1;
	if(maxy>(double)rows-1) maxy = rows-1;

	for (i = miny; i<= maxy; i++){
		for (j = minx; j<= maxx; j++){
			pt.x = (double)j;
			pt.y = (double)i;

			if(check_inside (pt1, pt2, pt3, pt4, pt)){
				pt.x = pt.x - center.x;
				pt.y = pt.y - center.y;
				di = rotatey(pt.x,pt.y,-t);
				dj = rotatex(pt.x,pt.y,-t);
				if ((-l_5 <= dj)&&(dj< -l_5 + w_n))
					y2 = ((dj-w_n+l_5)*(dj-w_n+l_5)+di*di)*mt_ww;
				else if ((-l_5 + w_n <= dj)&&(dj < l_5 - w_n))
					y2 = (di)*(di)*mt_ww;
				else
					y2 = ((dj+w_n-l_5)*(dj+w_n-l_5)+di*di)*mt_ww;
				dtmp = f_exp(ma, beta[0], y2);//*1.1-0.1; // *1.1-0.1 -> to make more zeros
				if (dtmp<beta_img[0][i][j])
					beta_img[0][i][j] = dtmp;
			}
		}
	}
}






double intersection_area(DPoint pt1, DPoint pt2, DPoint pt3, DPoint pt4, 
						 DPoint a, DPoint b, DPoint c, DPoint d, double t)
{
	int num, in_num, num2;
	double di, dj, w, l, w2, l2;
	DPoint pt;

	l = dist(pt2, pt1);
	w = dist(pt3, pt1);
	l2 = dist(b, a);
	w2 = dist(c, a);
	num2 = ((int)l2)*((int)w2);
	num = ((int)l)*((int)w);
	num = fmin(num, num2);

	in_num = 0;
	for(di = 0.; di<=w; di += 1.0){
		for(dj = 0.; dj<=l; dj += 1.0){

			pt.x = pt1.x + dj*cos(t) + di*sin(t);
			pt.y = pt1.y + dj*sin(t) - di*cos(t);
			if(check_inside (a, b, c, d, pt))
				in_num++;
		}
	}
	return (double)in_num/(double)num;
}



#define TAU_INT 10//4
#define ALPHA_THETA 0.6
#define LAMBDA_CONNECTION 1
int C_prior (NeckDent *mp, int k, int np_num, double epsilon, double **inter_e)
{	
	int i,m=0;
	double dtmp, dtmp2, a, b, norm = exp((double)TAU_INT)-1;
	double theta_diff;
	double d1=0,d2, d3, d4;
	DPoint epi1, epi2, epj1, epj2;
	
	for(i=0;i<np_num;i++)
	{
		if(i!=k)
		{
			dtmp = dist(mp[k].center,mp[i].center);
			if (dtmp<epsilon)
			{
				mp[k].multiple_E = INF;
				return 0;
			}
			a = dist(mp[k].r[0], mp[k].center);
			b = dist(mp[i].r[0], mp[i].center);
			if(a+b>dtmp)
			{ 		
				dtmp2 = intersection_area(mp[k].r[0], mp[k].r[1], mp[k].r[2], mp[k].r[3],
				mp[i].r[0], mp[i].r[1], mp[i].r[2], mp[i].r[3], mp[k].theta);
				if(mp[k].type == NDTYPE_BOTH)
				{

					epi1.x = (mp[i].r[0].x+mp[i].r[2].x)/2;
					epi1.y = (mp[i].r[0].y+mp[i].r[2].y)/2;
					epi2.x = (mp[i].r[1].x+mp[i].r[3].x)/2;
					epi2.y = (mp[i].r[1].y+mp[i].r[3].y)/2;
					epj1.x = (mp[k].r[0].x+mp[k].r[2].x)/2;
					epj1.y = (mp[k].r[0].y+mp[k].r[2].y)/2;
					epj2.x = (mp[k].r[1].x+mp[k].r[3].x)/2;
					epj2.y = (mp[k].r[1].y+mp[k].r[3].y)/2;
					d1 = dist(epi1,epj1);
					d2 = dist(epi1,epj2);
					d1 = fmin(d1,d2);
					d3 = dist(epi2,epj1);
					d4 = dist(epi2,epj2);
					d2 = fmin(d3,d4);
					d1 = fmin(d1,d2);
					theta_diff = fmin(fabs(mp[i].theta-mp[k].theta),M_PI-fabs(mp[i].theta-mp[k].theta));
					inter_e[k][i] = (exp(TAU_INT*dtmp2/(theta_diff+ALPHA_THETA))-1)+LAMBDA_CONNECTION*d1;///norm;
				}
				else
				{
					inter_e[k][i] = (exp(TAU_INT*dtmp2)-1);///norm;
				}
				
				inter_e[i][k] = inter_e[k][i];
				mp[i].multiple_E += inter_e[i][k];
				mp[k].multiple_E += inter_e[k][i];
			}
		}
	}
	return 1;
}


/*==================================================================================================*/

void draw_nds_both_2D(NeckDent *mp, MPP_Parameters mpp, unsigned char **image, int cols, int rows)
{
	unsigned char num;
	double amp, moffset;
	double di, dj, w, l;
	double w_5, l_5, mt_ww, tau = mpp.gaussian_tau;
	double y2;
	DPoint pt;
	
	l = mp->length;
	w = mp->width;
	w_5 = w/2.;
	l_5 = l/2.;
	mt_ww = -tau/(w*w);
	amp = mp->e[0];
	moffset = mp->e[2];

	for(di = -w_5; di<=w_5; di += STEP_DY){
		for(dj = -l_5; dj<=l_5; dj += STEP_DX){
			y2 = (di)*(di)*mt_ww;
			y2 = amp*(-exp(y2))+moffset;
			//y2 = -exp(y2);
			pt.x = (mp->r[0].x+mp->r[3].x)/2. + dj*cos(mp->theta) + di*sin(mp->theta);
			pt.y = (mp->r[0].y+mp->r[3].y)/2. + dj*sin(mp->theta) - di*cos(mp->theta);
			if((pt.x >= 0.)&&(pt.x < (double)cols-1.)&&(pt.y >= 0.)&&(pt.y < (double)rows-1.)){
				if(y2<0)
					num = 0;
				else if(y2>255)
					num = 255;
				else
					num = (unsigned char)y2;
				image[(int)pt.y][(int)pt.x] = num;
			}
		}
	}
}



int select_jump_type(double *interval)
{
	double a;
	int i;

	a = random2();

	for(i = 0;i<NUM_JUMP_TYPE ;i++)
	{
		if ((a>interval[i])&&(a<=interval[i+1])){
			return i;
		}
	}
	return 0;
}

#define VARI_SLOPE 0
/*Calculate T test between classes and return f(e) based off t-test */
int nd_Single_Object(unsigned char **yimg, NeckDent *mp, double *mean, double *vari, double variance,
					 MPP_Parameters mpp, int cols, int rows)
{
	
	int i,n;
	double error;

	for(i = 0; i<ERR_NUM; i++)
		mp->e[i] = 0;

	if (mp->type == NDTYPE_BOTH)
	{
		error  = avg_t_test(yimg, cols, rows, mean[0], vari[0], mean[1], vari[1], variance, mpp.gaussian_tau, 
		mpp.error_th, mp->theta, mp->r[0], mp->r[1],mp->r[2], mp->r[3], &n, mp->e);
		
		mp->single_E  = mpp.lambda_e*error;
	}
	#if 0
		else if (mp->type == NDTYPE_NECKING)
		{
			mp->single_E  = avg_Gaussian_neck3(yimg, cols, rows, mean[0], vari[0], mean[1], vari[1], variance, mpp.gaussian_tau,
							mp->theta, mp->r[0], mp->r[1],mp->r[2], mp->r[3], mpp.blur, mpp.blur_size, &n, &e0, &e1, &e2);
							
			mp->single_E  = exp(mpp.lambda_e*(mp->single_E - mpp.error_th))-(e0-mpp.amp_th);
		}
		else if (mp->type == NDTYPE_DENTING)
		{
			mp->single_E  = avg_Gaussian_dent3(yimg, cols, rows, mean[0], vari[0], mean[1], vari[1], variance, mpp.gaussian_tau,
			mpp.dent_l_w_ratio, mp->theta, mp->r[0], mp->r[1],mp->r[2], mp->r[3], mpp.blur, mpp.blur_size, &n, &e0, &e1, &e2);
			mp->single_E  = exp(mpp.lambda_e*(mp->single_E - mpp.error_th))-(e0-mpp.amp_th);
		}
	#else
		else if (mp->type == NDTYPE_NECKING)
		{
			error  = avg_t_test_neck(yimg, cols, rows, mean[0], vari[0], mean[1], vari[1], variance, mpp.gaussian_tau,
				mpp.error_th, mp->theta, mp->r[0], mp->r[1],mp->r[2], mp->r[3], &n, mp->e);
			mp->single_E  = mpp.lambda_e*error;//-mpp.lambda_l*(mp->length-mpp.length_th)+mpp.lambda_s*exp(mp->e[3]-mpp.symmetry_th)+mpp.lambda_nc*mp->e[5];
		}
		else if (mp->type == NDTYPE_DENTING)
		{
			error  = avg_t_test_dent(yimg, cols, rows, mean[0], vari[0], mean[1], vari[1], variance, mpp.gaussian_tau,
				mpp.error_th, mpp.dent_l_w_ratio, mp->theta, mp->r[0], mp->r[1],mp->r[2], mp->r[3], &n, mp->e);
			mp->single_E  = mpp.lambda_e*error + mpp.lambda_dc*mp->e[27]; //-mpp.lambda_l*(mp->length-mpp.length_th)+mpp.lambda_s*exp(mp->e[3]-mpp.symmetry_th)+mpp.lambda_dc*mp->e[5];
		}
	#endif
	

	else
	{
		return 0;
	}
	return 1;
}

/* Rotate Corner Points */
void nd_GenCornerPoint(NeckDent *mp)
{
	int k;
	DPoint r[12];

	r[0].x  = -mp->length/2;			r[0].y  =  mp->width/2;
	r[1].x  =  mp->length/2;			r[1].y  =  mp->width/2;
	r[2].x  = -mp->length/2;			r[2].y  =  -mp->width/2;
	r[3].x  =  mp->length/2;			r[3].y  =  -mp->width/2;
	for(k=0;k<4;k++){
		mp->r[k].x = rotatex(r[k].x,r[k].y,mp->theta)+mp->center.x;
		mp->r[k].y = rotatey(r[k].x,r[k].y,mp->theta)+mp->center.y;
	}
}

#if 1
//Give birth with random values of a and v depending on image.
void nd_Qbirth(int i, int j, NeckDent *mp, MPP_Parameters mpp)
{
	
	double a, b;
	if(mpp.nd_type_num == 0)
	{
		mp->type = NDTYPE_BOTH;
	}
	if(mpp.nd_type_num == 1)
	{
		mp->type = NDTYPE_NECKING;
	}
	else if(mpp.nd_type_num == 2)
	{
		mp->type = NDTYPE_DENTING;
	}
	else //if mpp.nd_type_num == 3
	{
		// rand() fn is [0,1]
		mp->type = (ND_Type)int_random(NDTYPE_NECKING,NDTYPE_DENTING);
	}
	
	mp->center.x = (double)j;
	mp->center.y = (double)i;
	do
	{
		mp->width = ubnd_random( mpp.widthmin, mpp.widthmax);
		mp->length = ubnd_random( mpp.lengthmin, mpp.lengthmax);
	}while(mp-> length <= 2*mp->width);
		
	
	if(mpp.fixed_param_num == 2)
	{
		a = random2();
		if (a < 0.5)
		{
			a = 0.20*M_PI;//-0.1;
			b = 0.20*M_PI;//+0.1;
		}
		else
		{
			a = 0.70*M_PI;//-0.1;
			b = 0.70*M_PI;//+0.1;
		}
		mp->theta = ubnd_random(a,b);
	}
	else if(mpp.fixed_param_num == 6)
	{
		a = random2();
		if (a < 0.4)
		{
			a = 0.31*M_PI;//-0.1;
			b = 0.31*M_PI;//+0.1;
		}
		else if (a < 0.8)
		{
			a = 0.77*M_PI;//-0.1;
			b = 0.77*M_PI;//+0.1;
		}
		else
		{
			a = 0.05*M_PI;//-0.1;
			b = 0.05*M_PI;//+0.1;
		}
		mp->theta = ubnd_random(a,b);
	}
	else
	{
		a = 0;
		if(mp->type == NDTYPE_DENTING)
			b = 2*M_PI;
		else
			b = M_PI;
		mp->theta = ubnd_random(a,b);
	}
	mp->single_E = 0;
	mp->multiple_E = 0;
	
}
#endif
#if 0
//Give birth with random values of a and v depending on image.
void nd_Qbirth(int i, int j, NeckDent *mp, MPP_Parameters mpp)
{
	
	double a, b;
	if(mpp.nd_type_num == 2)
	{
		// rand() fn is [0,1]
		mp->type = (ND_Type)int_random(NDTYPE_NECKING,NDTYPE_DENTING);
	}
	else if(mpp.nd_type_num == 1)
	{
		mp->type = NDTYPE_NECKING;
	}
	else
	{
		mp->type = NDTYPE_BOTH;
	}
	
	mp->center.x = (double)j;
	mp->center.y = (double)i;
	mp->width = ubnd_random( mpp.widthmin, mpp.widthmax);
	a = mpp.lengthmin;
	b = mpp.lengthmax;
	mp->length = ubnd_random( a, b);

	if(mpp.fixed_param_num == 2)
	{
		a = random2();
		if (a < 0.5)
		{
			a = 0.20*M_PI;//-0.1;
			b = 0.20*M_PI;//+0.1;
		}
		else
		{
			a = 0.70*M_PI;//-0.1;
			b = 0.70*M_PI;//+0.1;
		}
		mp->theta = ubnd_random(a,b);
	}
	else if(mpp.fixed_param_num == 6)
	{
		a = random2();
		if (a < 0.4)
		{
			a = 0.31*M_PI;//-0.1;
			b = 0.31*M_PI;//+0.1;
		}
		else if (a < 0.8)
		{
			a = 0.77*M_PI;//-0.1;
			b = 0.77*M_PI;//+0.1;
		}
		else
		{
			a = 0.05*M_PI;//-0.1;
			b = 0.05*M_PI;//+0.1;
		}
		mp->theta = ubnd_random(a,b);
	}
	else
	{
		a = 0;
		if(mp->type == NDTYPE_DENTING)
			b = 2*M_PI;
		else
			b = M_PI;
		mp->theta = ubnd_random(a,b);
	}
	mp->single_E = 0;
	mp->multiple_E = 0;
}
#endif





int jmp_birth(unsigned char **yimg, double **lm, unsigned char **MP_exist, NeckDent *mp, int *np_num, 
			   double beta, double epsilon, MPP_Parameters mpp,
			   int cols, int rows, double *mean, double *vari, double variance, double **inter_e)
{
	
	int i, j;
	double a, Qb, total_energy, e_ratio, G_ratio, bd_ratio = mpp.p_death/mpp.p_birth*(mpp.alm*mpp.vk);

	if ((*np_num)>=MAX_MKPNT_NUM-1)
		return -1;

	do{
		j = int_random(0,cols-1);
		i = int_random(0,rows-1);
	}while(MP_exist[i][j]||(lm[i][j]==0));

	nd_Qbirth(i, j, &(mp[*np_num]), mpp);
	nd_GenCornerPoint(&(mp[*np_num]));
	
	// Green Ratio
	nd_Single_Object(yimg, &(mp[*np_num]), mean, vari, variance, mpp, cols, rows);
	
	C_prior (mp, *np_num, *np_num, epsilon, inter_e);

	total_energy =  mpp.alpha*mp[*np_num].single_E+mpp.lambda_int*mp[*np_num].multiple_E;
	e_ratio = exp(-(beta*total_energy));
	Qb = bd_ratio/((*np_num)+1);
	G_ratio = fmin(1, e_ratio*Qb);

	a = random2();
	if (a < G_ratio)	// accept
	{
		MP_exist[i][j] = 1;
		mp[*np_num].num = *np_num;
		mp[*np_num].state = STATE_NEW_BORN;
		(*np_num)++;
		return 1;
	}
	else{
		(*np_num)++;
		Qdelete(mp, (*np_num)-1, np_num, inter_e);
		return 0;
	}

}

int jmp_death(unsigned char **yimg, unsigned char **MP_exist, NeckDent *mp, int *np_num, 
			   double beta, MPP_Parameters mpp,
			   int cols, int rows, double *mean, double *vari, double variance,
			   double **inter_e)
{
	int k;
	double a, Qd, total_energy, e_ratio, G_ratio, db_ratio = mpp.p_birth/mpp.p_death/(mpp.alm*mpp.vk);

	if ((*np_num) <=0)
		return -1;

	k = int_random(0,(*np_num)-1);
 



	// Green Ratio
	total_energy =  mpp.alpha*mp[k].single_E + mpp.lambda_int*mp[k].multiple_E;
	e_ratio = exp(beta*total_energy);
	Qd = db_ratio*(*np_num);
	G_ratio = fmin(1, e_ratio*Qd);

	a = random2();

	if (a < G_ratio)	// accept to delete
	{

		MP_exist[(int)mp[k].center.y][(int)mp[k].center.x] = 0;
		Qdelete(mp, k, np_num, inter_e);
		return 1;
	}
	return 0;
}

int jmp_translation(unsigned char **yimg, NeckDent *mp,int np_num, double beta, double epsilon,
					MPP_Parameters mpp, 
					 int cols, int rows, double *mean, double *vari, double variance, double **inter_e)
{
	double a, b, dx, dy, e_ratio, G_ratio, total_energy;
	int i, k;
	NeckDent tmp_mp, new_mp;
	double inter_e_tmp[MAX_MKPNT_NUM];

	if (np_num <=0)
		return -1;

	k = int_random(0,np_num-1);

	tmp_mp = mp[k];	// backup
	new_mp = tmp_mp;
	new_mp.multiple_E = 0;

	a = -DELTA_TRANSLATION;
	b = DELTA_TRANSLATION;
	dx = ubnd_random(a,b);
	dy = ubnd_random(a,b);

	new_mp.center.x = new_mp.center.x + dx;
	if ((new_mp.center.x>=0)&&(new_mp.center.x<=(cols-1))){
		new_mp.center.y = new_mp.center.y + dy;
		if ((new_mp.center.y>=0)&&(new_mp.center.y<=(rows-1))){
			// calculate Green ratio
			for(i=0;i<k;i++){
				inter_e_tmp[i] = inter_e[k][i];
			}
			for(i=k;i<np_num-1;i++){
				inter_e_tmp[i] = inter_e[k][i+1];
			}
			inter_e_tmp[np_num-1] = 0;

			Qdelete(mp, k, &np_num, inter_e); // np_num--;
			np_num++;
			nd_GenCornerPoint(&(new_mp));
			mp[np_num-1] = new_mp;
			nd_Single_Object(yimg, &(mp[np_num-1]), mean, vari, variance, mpp, cols, rows);
			C_prior (mp, np_num-1, np_num, epsilon, inter_e);
			total_energy = mpp.alpha*(mp[np_num-1].single_E-tmp_mp.single_E)+mpp.lambda_int*(mp[np_num-1].multiple_E-tmp_mp.multiple_E);
			e_ratio = exp(-beta*total_energy);
			G_ratio = fmin(1, e_ratio);

			a = random2();

			if (a >= G_ratio)	// reject
			{
				Qdelete(mp, np_num-1, &np_num, inter_e);
				np_num++;
				mp[np_num-1] = tmp_mp;	// restore
				for(i=0;i<np_num;i++){
					inter_e[np_num-1][i] = inter_e_tmp[i];
					inter_e[i][np_num-1] = inter_e_tmp[i];
					mp[i].multiple_E += inter_e_tmp[i];
				}
				return 0;
			}
			return 1;
		}
		return -1;
	}
	return -1;
}

int jmp_dilation(unsigned char **yimg, NeckDent *mp,int np_num, double beta, double epsilon,
					MPP_Parameters mpp, 
				  int cols, int rows, double *mean, double *vari, double variance, double **inter_e)
{
	double a, b, dx, dy, e_ratio, G_ratio, total_energy;
	int i, k;
	NeckDent tmp_mp, new_mp;
	double inter_e_tmp[MAX_MKPNT_NUM];

	if (np_num <=0)
		return -1;

	k = int_random(0,np_num-1);

	tmp_mp = mp[k];	// backup
	new_mp = tmp_mp;
	new_mp.multiple_E = 0;

	a = -DELTA_DILATION;
	b = DELTA_DILATION;
	dx = ubnd_random(a,b);
	dy = ubnd_random(a,b);

	new_mp.width = new_mp.width + dx;
	if ((new_mp.width>=mpp.widthmin)&&(new_mp.width<=mpp.widthmax)){
		new_mp.length = new_mp.length + dy;
	if((new_mp.length>=mpp.lengthmin)&&(new_mp.length<=mpp.lengthmax)){

			// calculate Green ratio
			for(i=0;i<k;i++){
				inter_e_tmp[i] = inter_e[k][i];
			}
			for(i=k;i<np_num-1;i++){
				inter_e_tmp[i] = inter_e[k][i+1];
			}
			inter_e_tmp[np_num-1] = 0;
			Qdelete(mp, k, &np_num, inter_e); // np_num--;
			np_num++;
			nd_GenCornerPoint(&(new_mp));
			mp[np_num-1] = new_mp;
			nd_Single_Object(yimg, &(mp[np_num-1]), mean, vari, variance, mpp, cols, rows);
			C_prior (mp, np_num-1, np_num, epsilon, inter_e);
			total_energy = mpp.alpha*(mp[np_num-1].single_E-tmp_mp.single_E)+mpp.lambda_int*(mp[np_num-1].multiple_E-tmp_mp.multiple_E);
			e_ratio = exp(-beta*total_energy);
			G_ratio = fmin(1, e_ratio);


			a = random2();

			if (a >= G_ratio)	// reject
			{
				Qdelete(mp, np_num-1, &np_num, inter_e);
				np_num++;
				mp[np_num-1] = tmp_mp;	// restore
				for(i=0;i<np_num;i++){
					inter_e[np_num-1][i] = inter_e_tmp[i];
					inter_e[i][np_num-1] = inter_e_tmp[i];
					mp[i].multiple_E += inter_e_tmp[i];
				}
				return 0;
			}
			return 1;
		}
		return -1;
	}
	return -1;
}

int jmp_rotation(unsigned char **yimg, NeckDent *mp,int np_num, double beta, double epsilon,
					MPP_Parameters mpp, 
				  int cols, int rows, double *mean, double *vari, double variance, double **inter_e)
{
	double a, b, dx, e_ratio, G_ratio, total_energy;
	int i, k;
	NeckDent tmp_mp, new_mp;
	double inter_e_tmp[MAX_MKPNT_NUM];

	if (np_num <=0)
		return -1;

	k = int_random(0,np_num-1);
	tmp_mp = mp[k];	// backup
	new_mp = tmp_mp;
	new_mp.multiple_E = 0;

	a = -DELTA_ROTATION;
	b = DELTA_ROTATION;
	dx = ubnd_random(a,b);

	new_mp.theta = new_mp.theta + dx;
	if(new_mp.type == NDTYPE_DENTING){
		if (new_mp.theta<0)
			new_mp.theta = new_mp.theta + 2*M_PI;
		else if (new_mp.theta>2*M_PI)
			new_mp.theta = new_mp.theta - 2*M_PI;
	}
	else{
		if (new_mp.theta<0)
			new_mp.theta = new_mp.theta + M_PI;
		else if (new_mp.theta>M_PI)
			new_mp.theta = new_mp.theta - M_PI;
	}

	if (1){
		// calculate Green ratio
		for(i=0;i<k;i++){
			inter_e_tmp[i] = inter_e[k][i];
		}
		for(i=k;i<np_num-1;i++){
			inter_e_tmp[i] = inter_e[k][i+1];
		}
		inter_e_tmp[np_num-1] = 0;
		Qdelete(mp, k, &np_num, inter_e); // np_num--;
		np_num++;
		nd_GenCornerPoint(&(new_mp));
		mp[np_num-1] = new_mp;
		nd_Single_Object(yimg, &(mp[np_num-1]), mean, vari, variance, mpp, cols, rows);
		C_prior (mp, np_num-1, np_num, epsilon, inter_e);
		total_energy = mpp.alpha*(mp[np_num-1].single_E-tmp_mp.single_E)+mpp.lambda_int*(mp[np_num-1].multiple_E-tmp_mp.multiple_E);
		e_ratio = exp(-beta*total_energy);
		G_ratio = fmin(1, e_ratio);
	
		a = random2();

		if (a >= G_ratio)	// reject
		{
			Qdelete(mp, np_num-1, &np_num, inter_e);
			np_num++;
			mp[np_num-1] = tmp_mp;	// restore
			for(i=0;i<np_num;i++){
				inter_e[np_num-1][i] = inter_e_tmp[i];
				inter_e[i][np_num-1] = inter_e_tmp[i];
				mp[i].multiple_E += inter_e_tmp[i];
			}
			return 0;
		}
		return 1;
	}
	return -1;
}



int switch_neckdent(NeckDent tmp_mp, NeckDent *new_mp, MPP_Parameters mpp, int cols, int rows, int switch_type)
{
	double da, db, dl, udl, udl2, fact;

	*new_mp = tmp_mp;
	new_mp->multiple_E = 0;
	fact = ubnd_random(0,2);
	udl = fact*new_mp->width/4.;
	udl2 = udl/2.;
	switch (switch_type)
	{
		case 0: // dent to neck with shorting
			da = new_mp->center.x + udl2*cos(new_mp->theta);
			if((da>=0)&&(da<=(cols-1))){
				new_mp->center.x = da;
				db = new_mp->center.y + udl2*sin(new_mp->theta);
				if((db>=0)&&(db<=(rows-1))){
					new_mp->center.y = db;
					dl = new_mp->length-udl;
					if(dl>=mpp.lengthmin){
						new_mp->length = dl;
						if(new_mp->theta>M_PI)
							new_mp->theta = new_mp->theta - M_PI;
					}
					else
						return -1;
				}
				else
					return -1;
			}
			else
				return -1;
			new_mp->type = NDTYPE_NECKING;
			break;
	
		case 1: // dent to neck with lengthening
			da = new_mp->center.x - udl2*cos(new_mp->theta);
			if((da>=0)&&(da<=(cols-1))){
				new_mp->center.x = da;
				db = new_mp->center.y - udl2*sin(new_mp->theta);
				if((db>=0)&&(db<=(rows-1))){
					new_mp->center.y = db;
					dl = new_mp->length+udl;
					if(dl<=mpp.lengthmax){
						new_mp->length = dl;
						if(new_mp->theta>M_PI)
							new_mp->theta = new_mp->theta - M_PI;
					}
					else
						return -1;
				}
				else
					return -1;
			}
			else
				return -1;
			new_mp->type = NDTYPE_NECKING;
			break;

		case 2: // neck to dent with lengthening 1
			da = new_mp->center.x - udl2*cos(new_mp->theta);
			if((da>=0)&&(da<=(cols-1))){
				new_mp->center.x = da;
				db = new_mp->center.y - udl2*sin(new_mp->theta);
				if((db>=0)&&(db<=(rows-1))){
					new_mp->center.y = db;
					dl = new_mp->length+udl;
					if(dl<=mpp.lengthmax){
						new_mp->length = dl;
					}
					else
						return -1;
				}
				else
					return -1;
			}
			else
				return -1;
			new_mp->type = NDTYPE_DENTING;
			break;

		case 3: // neck  to dent with lengthening 2
			da = new_mp->center.x + udl2*cos(new_mp->theta);
			if((da>=0)&&(da<=(cols-1))){
				new_mp->center.x = da;
				db = new_mp->center.y + udl2*sin(new_mp->theta);
				if((db>=0)&&(db<=(rows-1))){
					new_mp->center.y = db;
					dl = new_mp->length+udl;
					if(dl<=mpp.lengthmax){
						new_mp->length = dl;
						new_mp->theta = new_mp->theta + M_PI;

					}
					else
						return -1;
				}
				else
					return -1;
			}
			else
				return -1;
			new_mp->type = NDTYPE_DENTING;
			break;

		case 4: // neck to dent with shorting 1
			da = new_mp->center.x + udl2*cos(new_mp->theta);
			if((da>=0)&&(da<=(cols-1))){
				new_mp->center.x = da;
				db = new_mp->center.y + udl2*sin(new_mp->theta);
				if((db>=0)&&(db<=(rows-1))){
					new_mp->center.y = db;
					dl = new_mp->length-udl;
					if(dl>=mpp.lengthmin){
						new_mp->length = dl;
					}
					else
						return -1;
				}
				else
					return -1;
			}
			else
				return -1;
			new_mp->type = NDTYPE_DENTING;
			break;

		case 5: // neck to dent with shorting 2
			da = new_mp->center.x - udl2*cos(new_mp->theta);
			if((da>=0)&&(da<=(cols-1))){
				new_mp->center.x = da;
				db = new_mp->center.y - udl2*sin(new_mp->theta);
				if((db>=0)&&(db<=(rows-1))){
					new_mp->center.y = db;
					dl = new_mp->length-udl;
					if(dl>=mpp.lengthmin){
						new_mp->length = dl;
						new_mp->theta = new_mp->theta + M_PI;
						new_mp->theta = new_mp->theta - 2*M_PI;
					}
					else
						return -1;
				}
				else
					return -1;
			}
			else
				return -1;
			new_mp->type = NDTYPE_DENTING;
			break;

		default:
			return -1;
			break;
	}
	return 1;
}


int jmp_switch(unsigned char **yimg, NeckDent *mp,int np_num, double beta, double epsilon,
				MPP_Parameters mpp, 
				int cols, int rows, double *mean, double *vari, 
				double variance, double **inter_e)
{
	double a, e_ratio, G_ratio, occur_rate, total_energy;
	int i, k;
	NeckDent tmp_mp, new_mp;
	double inter_e_tmp[MAX_MKPNT_NUM];

	if(mpp.nd_type_num == 1){
		return -1;
	}
	if (np_num <=0)
		return -1;



	#define J_DENT2NECK		0.5
	#define J_NECK2DENT		0.5
	a = random2();
	if(a<J_DENT2NECK/2.){	// dent to neck1 (l-w/4)
		i=0;
		do{
			k = int_random(0,np_num-1);
			i++;
		}while((mp[k].type != NDTYPE_DENTING)&&(i<1000));
		if (i>=1000)
			return -1;
		tmp_mp = mp[k];	// backup
		if (switch_neckdent(tmp_mp, &new_mp, mpp, cols, rows, 0)==-1)
			return -1;
		occur_rate = J_NECK2DENT/J_DENT2NECK;
	}
	else if(a<J_DENT2NECK){  // dent to neck2 (l+w/4)
		i=0;
		do{
			k = int_random(0,np_num-1);
			i++;
		}while((mp[k].type != NDTYPE_DENTING)&&(i<1000));
		if (i>=1000)
			return -1;
		tmp_mp = mp[k];	// backup
		if(switch_neckdent(tmp_mp, &new_mp, mpp, cols, rows, 1)==-1)
			return -1;
		occur_rate = J_NECK2DENT/J_DENT2NECK;
	}
	else if(a<J_DENT2NECK+J_NECK2DENT/4.){	// neck to dent
		i=0;
		do{
			k = int_random(0,np_num-1);
			i++;
		}while((mp[k].type != NDTYPE_NECKING)&&(i<1000));
		if (i>=1000)
			return -1;
		tmp_mp = mp[k];	// backup
		if(switch_neckdent(tmp_mp, &new_mp, mpp, cols, rows, 2)==-1)
			return -1;
		occur_rate = J_DENT2NECK/J_NECK2DENT;
	}
	else if(a<J_DENT2NECK+J_NECK2DENT/2.){	// neck to dent2
		i=0;
		do{
			k = int_random(0,np_num-1);
			i++;
		}while((mp[k].type != NDTYPE_NECKING)&&(i<1000));
		if (i>=1000)
			return -1;
		tmp_mp = mp[k];	// backup
		if(switch_neckdent(tmp_mp, &new_mp, mpp, cols, rows, 3)==-1)
			return -1;
		occur_rate = J_DENT2NECK/J_NECK2DENT;
	}
	else if(a<J_DENT2NECK+J_NECK2DENT*3./4.){	// neck to dent3
		i=0;
		do{
			k = int_random(0,np_num-1);
			i++;
		}while((mp[k].type != NDTYPE_NECKING)&&(i<1000));
		if (i>=1000)
			return -1;
		tmp_mp = mp[k];	// backup
		if(switch_neckdent(tmp_mp, &new_mp, mpp, cols, rows, 4)==-1)
			return -1;
		occur_rate = J_DENT2NECK/J_NECK2DENT;
	}
	else {	// neck to dent4
		i=0;
		do{
			k = int_random(0,np_num-1);
			i++;
		}while((mp[k].type != NDTYPE_NECKING)&&(i<1000));
		if (i>=1000)
			return -1;
		tmp_mp = mp[k];	// backup
		if(switch_neckdent(tmp_mp, &new_mp, mpp, cols, rows, 5)==-1)
			return -1;
		occur_rate = J_DENT2NECK/J_NECK2DENT;
	}

	// calculate Green ratio
	for(i=0;i<k;i++){
		inter_e_tmp[i] = inter_e[k][i];
	}
	for(i=k;i<np_num-1;i++){
		inter_e_tmp[i] = inter_e[k][i+1];
	}
	inter_e_tmp[np_num-1] = 0;
	Qdelete(mp, k, &np_num, inter_e); // np_num--;
	np_num++;
	nd_GenCornerPoint(&(new_mp));
	mp[np_num-1] = new_mp;
	nd_Single_Object(yimg, &(mp[np_num-1]), mean, vari, variance, mpp, cols, rows);
	C_prior (mp, np_num-1, np_num, epsilon, inter_e);
	total_energy = mpp.alpha*(mp[np_num-1].single_E-tmp_mp.single_E)+mpp.lambda_int*(mp[np_num-1].multiple_E-tmp_mp.multiple_E);
	e_ratio = exp(-beta*total_energy);
	G_ratio = fmin(1, occur_rate*e_ratio);

	a = random2();

	if (a >= G_ratio)	// reject
	{
		Qdelete(mp, np_num-1, &np_num, inter_e);
		np_num++;
		mp[np_num-1] = tmp_mp;	// restore
		for(i=0;i<np_num;i++){
			inter_e[np_num-1][i] = inter_e_tmp[i];
			inter_e[i][np_num-1] = inter_e_tmp[i];
			mp[i].multiple_E += inter_e_tmp[i];
		}
		return 0;
	}
	return 1;
}


void qsort_with_multienergy(NeckDent *mp, int np_num, int size, double **inter_e, int (*compare_fn)(const void * ,const void * ))
{
	int i,j;
	double tmp[MAX_MKPNT_NUM][MAX_MKPNT_NUM];

	for(i=0;i<np_num;i++)
		mp[i].sort_idx = i;
	qsort (mp, np_num, size, compare_fn);
	for(i=0;i<np_num;i++)
		for(j=0;j<np_num;j++){
			tmp[i][j] = inter_e[mp[i].sort_idx][j];
		}
	for(j=0;j<np_num;j++)
		for(i=0;i<np_num;i++){
			inter_e[i][j] = tmp[i][mp[j].sort_idx];
		}
}


int nd_mpp_multiple_birth_n_death(unsigned char **yimg, double **lm, double *mean,
			double *vari, double variance, NeckDent *mp,
			MPP_Parameters mpp, double *total_e, int *mp_num, 
			int cols, int rows)
{
	
	unsigned char **MP_exist;
	int i, j, k, iter;
	double dtmp;
	double birth_rate, d_beta, d_rate;
	double epsilon = HARDCORE_REPULSION;
	int np_num;
	char filename[1024];
	double **inter_e, total_energy;
	double betampp = mpp.betampp;
	double delta = 0.9;

	MP_exist = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
	inter_e = (double **)get_img(MAX_MKPNT_NUM, MAX_MKPNT_NUM, sizeof(double));
	
	/* Initialize MP_Exist to 0*/
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			MP_exist[i][j] = 0;
		}
		

	
	for(i = 0; i < MAX_MKPNT_NUM; i++){
		mp[i].state = STATE_NON_EXIST;
	}
	
	for (i = 0; i < MAX_MKPNT_NUM; i++)
		for (j = 0; j < MAX_MKPNT_NUM; j++){
			inter_e[i][j] = 0;
		}
	np_num = 0;


	
	for (iter = 1;iter <= mpp.iter_num; iter++)
	{ //MAX_MPP_ITER_NUM
		birth_rate = delta*mpp.b_zero;
		printf("iter=%d\n",iter);

		//----------------------- Birth -------------------------//
		for (i = 5;i<rows-5;i++){
			for (j = 5;j<cols-5;j++){
				if (MP_exist[i][j]==0){
					dtmp = ((double)rand()/RAND_MAX);
					
					if (dtmp < birth_rate*lm[i][j])
					{ // birth
						if (np_num>MAX_MKPNT_NUM-2)
						{
							printf("BREAK ALERT, i: %d j: %d \n", i,j);
							break;
						}
						/*Create Semi-Random Parameters for w,l,theta */
						nd_Qbirth(i, j, &(mp[np_num]), mpp);
						
						/* Rotate Corner Points */
						nd_GenCornerPoint(&(mp[np_num]));
						
						/* Find Likelyhood */
						nd_Single_Object(yimg, &(mp[np_num]), mean, vari, variance, mpp, cols, rows);
						
						if (mp[np_num].single_E < 0.15)
						{
							MP_exist[i][j] = 1;
							mp[np_num].num = np_num;
							mp[np_num].state = STATE_NEW_BORN;
							mp[np_num].multiple_E = 0.;
							np_num ++;
						}
					}//End Birth
				}//End if it doesnt exist
			}
		}

		//----------------------- death -------------------------//
		qsort_with_multienergy (mp, np_num, sizeof(NeckDent), inter_e, compare_Single_E);
		for(k = 0; k<np_num; k++){
			if(mp[k].state == STATE_NEW_BORN)
				C_prior(mp, k, np_num, epsilon, inter_e);
		}

		for(k = 0; k<np_num; k++){
			if(mp[k].multiple_E>=INF) {
				d_rate = 1.;
			}
			else{
				total_energy =  mpp.alpha*mp[k].single_E+mpp.lambda_int*mp[k].multiple_E;
				d_beta = exp(betampp*total_energy);
				if (d_beta > INF/2) // if d_beta is too big
					d_rate = 1;
				else
					d_rate = (delta*d_beta)/(1.0+delta*d_beta);
//				printf("k=%d, dp=%1.2f\n",k,d_rate);
			}

			dtmp = ((double)rand()/RAND_MAX);
			if (dtmp < d_rate){
				MP_exist[(int)mp[k].center.y][(int)mp[k].center.x] = 0;
				Qdelete(mp, k, &np_num, inter_e);
				k--;
			}
			else{
				mp[k].num = k;
				mp[k].state = STATE_EXIST;
			}
		}
		//---------------- Set up new parameters ----------------//
		delta = delta*mpp.de_coeff;
		betampp = betampp/mpp.de_coeff;
		//----------- for energy and mp number graph ------------//
		total_e[iter] = 0;
		for(i=0;i<np_num;i++)
			total_e[iter] +=  mpp.alpha*mp[i].single_E+mpp.lambda_int*mp[i].multiple_E/2.;
		mp_num[iter] = np_num;
		printf("obj num =%d, betampp = %1.2f\n",np_num,betampp);
		//if(iter%mpp.iter_num/2 == 0)
		//{
			
		//}
		
		
		
	} 
	free_img((void **)MP_exist);
	free_img((void **)inter_e);
	return np_num;
	}


int nd_mpp_rjmcmc(unsigned char **yimg, double **lm, double *mean, double *vari, 
			double variance, NeckDent *mp, 
			MPP_Parameters mpp, double *total_e, int *mp_num, 
			int cols, int rows)
{
	unsigned char **MP_exist;
	int i, j, iter;
	double dtmp;
	double epsilon = mpp.hard_repulsion;
	int np_num;
	char filename[1024];
	int jmp_type;
//	double ch_y[MAX_MKPNT_NUM], ch_x[MAX_MKPNT_NUM];
	FILE *fp;
	double accum_jmp_prob[NUM_JUMP_TYPE+1];
	double betaT0 = mpp.betampp/mpp.T0;
	double betaT = betaT0;
	double **inter_e;
//	double avg_U, vari_U;
	int	   return_val;
	
	accum_jmp_prob[0]= 0.;
	accum_jmp_prob[1]= accum_jmp_prob[0] + mpp.p_birth;
	accum_jmp_prob[2]= accum_jmp_prob[1] + mpp.p_death;
	accum_jmp_prob[3]= accum_jmp_prob[2] + mpp.p_translation;
	accum_jmp_prob[4]= accum_jmp_prob[3] + mpp.p_dilation;
	accum_jmp_prob[5]= accum_jmp_prob[4] + mpp.p_rotation;
	accum_jmp_prob[NUM_JUMP_TYPE]= 1.;


	MP_exist = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
	inter_e = (double **)get_img(MAX_MKPNT_NUM, MAX_MKPNT_NUM, sizeof(double));
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			MP_exist[i][j] = 0;
		}
		
	for (i = 0; i < MAX_MKPNT_NUM; i++)
		for (j = 0; j < MAX_MKPNT_NUM; j++){
			inter_e[i][j] = 0;
		}
	for(i = 0; i < MAX_MKPNT_NUM; i++){
		mp[i].state = STATE_NON_EXIST;
	}
	

	np_num = 0;
	sprintf(filename, "mpp.txt");
	if ((fp = fopen(filename, "wb")) == NULL ) {
		printf("Cannot open file %s\n", filename);
		exit(1);
	}
		
	for (iter = 0;iter < mpp.iter_num; iter++){ 
		jmp_type = select_jump_type(accum_jmp_prob);
		switch(jmp_type)
		{
			case JMP_BIRTH:
				return_val = jmp_birth(yimg, lm, MP_exist,mp, &np_num, betaT, epsilon, mpp, cols, rows, mean, vari, variance,inter_e);
				break;
			case JMP_DEATH:
				return_val = jmp_death(yimg, MP_exist, mp, &np_num, betaT, mpp, cols, rows, mean, vari, variance, inter_e);
				break;
			case JMP_TRANSLATION:
				return_val = jmp_translation(yimg, mp, np_num, betaT, epsilon, mpp, cols, rows, mean, vari, variance, inter_e);
				break;
			case JMP_DILATION:
				return_val = jmp_dilation(yimg, mp,np_num, betaT, epsilon, mpp, cols, rows, mean, vari, variance, inter_e);
				break;
			case JMP_ROTATION:
				return_val = jmp_rotation(yimg, mp,np_num, betaT, epsilon, mpp, cols, rows, mean, vari, variance, inter_e);
				break;
			case JMP_SWITCHING:
				return_val = jmp_switch(yimg, mp,np_num, betaT, epsilon, mpp, cols, rows, mean, vari, variance, inter_e);
				break;
			default:
				break;
		}
		
		
		// simulated annealing
		betaT = betaT0/(pow(mpp.de_coeff, (double)iter));


		if(iter%(RJMCMC_ITER_DIV)==0){
			total_e[iter/RJMCMC_ITER_DIV] = 0;
			for(i=0;i<np_num;i++){
				dtmp = mpp.alpha*mp[i].single_E+mpp.lambda_int*mp[i].multiple_E/2.;
				total_e[iter/RJMCMC_ITER_DIV] += dtmp; 
			}
			mp_num[iter/RJMCMC_ITER_DIV] = np_num;
		}
		
		if(iter%(RJMCMC_ITER_DIV)==0){		
			fprintf(fp,"%07d\r\n ", iter);
			for (j=0;j<np_num;j++){
				fprintf(fp,"%1.4f, ", mp[j].e[1]);
			}
			fprintf(fp,"\r\n");
			for (j=0;j<np_num;j++){
				fprintf(fp,"%1.4f, ", mp[j].e[0]);
			}
			fprintf(fp,"\r\n");
			fprintf(fp,"\r\n");
		}
	} 

	fclose(fp);



	for(i=0;i<np_num;i++){
		mp[i].num = i;
	}


	free_img((void **)MP_exist);
	free_img((void **)inter_e);
	
	return np_num;
}

/******************************************************************************/
// calculate_betaimg : 
// 
/******************************************************************************/
void calculate_betaimg(double **beta_img[], double *beta, NeckDent *mp, 
					   int np_num, MPP_Parameters mpp, int cols, int rows)
{
	int i, j, k;


	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			beta_img[0][i][j] = beta[1];
			beta_img[1][i][j] = beta[1];
		}
	}

	for(k = 0; k<np_num; k++){
		Gaussian_beta(beta_img, beta, cols, rows, &(mp[k]), mpp);
	}
}
