#include "QualityCandy.h"


void UpdateItsNeighboors_born_freeSeg (lineObj* seg, Candy* M)
{
	lineObj* target;

	if (seg->enda_L_Num != 0)
	{
		for (int i = 0; i < seg->enda_L_Num;i++)
		{
			target = seg->enda_L[i];
			updateNeighboorLink_born(M,target,seg,seg->enda);
		}
	}
	if (seg->endb_L_Num != 0)
	{
		for (int i = 0; i < seg->endb_L_Num;i++)
		{
			target = seg->endb_L[i];
			updateNeighboorLink_born(M,target,seg,seg->endb);
		}
	}
}

void UpdateItsNeighboors_death_freeSeg (lineObj* seg, Candy* M)
{
	lineObj* target;

	if (seg->enda_L_Num != 0)
	{
		for (int i = 0; i < seg->enda_L_Num;i++)
		{
			target = seg->enda_L[i];
			updateNeighboorLink_death(M,target,seg,seg->enda);
		}
	}
	if (seg->endb_L_Num != 0)
	{
		for (int i = 0; i < seg->endb_L_Num;i++)
		{
			target = seg->endb_L[i];
			updateNeighboorLink_death(M,target,seg,seg->endb);
		}
	}
}




void AddFreeSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,int test,double **patch, int patch_len, double ***Matrix, double prior_pen1, double prior_pen2, int **dilated_img)
{
	
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;
	double lambda = M->lambda;

	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));

	freeSeg->width = (int)floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5);
	freeSeg->len = (int)floor(random2()*(L_MAX-L_MIN)+L_MIN+0.5);
	
	freeSeg->theta = random2()*(THETA_MAX-THETA_MIN)+THETA_MIN;
	
	//freeSeg->theta = 0.2;

//KDW	freeSeg->enda.x  = (int)floor(random2()*(img_width-4*L_MAX)+2*L_MAX+0.5);
//KDW	freeSeg->enda.y = (int)floor(random2()*(img_height-4*L_MAX)+2*L_MAX+0.5);
	freeSeg->enda.x  = (int)floor(random2()*(img_width-2)+0.5);
	freeSeg->enda.y = (int)floor(random2()*(img_height-2)+0.5);

	int img_num = 0;
	int w_num = freeSeg->width-W_MIN;

	double iterval = _PI/(RADIUS_SEGS);
	while (freeSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;


#if 0 //KDW	 
	while(dilated_img[freeSeg->enda.y][freeSeg->enda.x] != 1)
	{
	//freeSeg->enda.x  = floor(random2()*(img_width-4*L_MAX)+2*L_MAX+0.5);
	//freeSeg->enda.y = floor(random2()*(img_height-4*L_MAX)+2*L_MAX+0.5);
	freeSeg->enda.x  = (int)floor(random2()*(img_width-2*5)+5+0.5);
	freeSeg->enda.y = (int)floor(random2()*(img_height-2*5)+5+0.5);
	}
#endif
	//test
	//	freeSeg->enda.x = 344;
	 //   freeSeg->enda.y = 104;
	//	freeSeg->theta = _PI-1.2391724;
	//	freeSeg->len = 43;
	//	freeSeg->width = 6;

	freeSeg->endb = GenerateEndb(freeSeg->enda,freeSeg->len,freeSeg->theta,img_height,img_width);

	freeSeg->len = (int)sqrt(double(DistSquare(freeSeg->endb,freeSeg->enda)));

	if (freeSeg->len >= L_MIN)
	{
	freeSeg->x = (int)floor(0.5*(double(freeSeg->enda.x+freeSeg->endb.x))+0.5);
	freeSeg->y = (int)floor(0.5*(double(freeSeg->enda.y+freeSeg->endb.y))+0.5);
	if(lm[freeSeg->y][freeSeg->x]==0){
		free(freeSeg);
		return;
	}

	freeSeg->endb_L_Num = 0;
	freeSeg->enda_L_Num = 0;
	freeSeg->enda_C_Num = 0;
	freeSeg->endb_C_Num = 0;
	freeSeg->img_num = img_num;

	double searchRatio = 0.25;
	double g_Rio = 0;
	double g_Rc =0.;

	int n_io = Bad_IO(freeSeg, M,&g_Rio);
	int n_eo = Bad_EO_freeSeg(NEIGHBOORHOOD,searchRatio,freeSeg,M,&g_Rc);

	//should add data term here, assume homogeneous Poisson 

	//double dterm = dataterm(mid_img[img_num],freeSeg->enda,freeSeg->endb);
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	//double dterm = dataterm(img,freeSeg->enda,freeSeg->endb,freeSeg->theta,freeSeg->len,img_num,patch,patch_len);
	//printf("%f\n",dterm);
	//if(dterm < 0) 
	//	printf("%f\n",dterm);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*n_io-w_eo*n_eo);
#endif
	freeSeg->dataterm = dterm;
	freeSeg->engergy_for_transition = Echange;

	double R = pow(Echange,T) *(M->p_f_d)* (M->lambda) * (L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN)/((M->p_f_b)*(n_f+1));
	//printf("free_seg_r = %f\n",R);
	double r = random2();
	if (r < MIN(1,R))
	{
		if (freeSeg->enda_C_Num == 0 &&  freeSeg->endb_C_Num == 0)
		{
//	if(g_Rio>10000)
//		printf("..");
			M->n_f++;
			freeSeg->type = 0;
		    LinkedListInsert( M->link_f,M->n_f, freeSeg);  //this is a free segment	
			UpdateItsNeighboors_born_freeSeg(freeSeg,M);
		}
		else
			free(freeSeg);
	}
	else
		free(freeSeg);
	}
	else
		free(freeSeg);
}
void KillFreeSeg(Candy *M, double **img, int img_height, int img_width, double T)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;
	double lambda = M->lambda;

	if (n_f != 0)
	{
		int num = (int)floor((n_f-1)*random2()+1+0.5);
		Node *p;
		Node *pre = M->link_f;

		for (int i = 1;i<num;i++)
			pre = pre->next;

		p = pre->next;

		lineObj *l = p->index;

		double g_Rio = 0; 
		int n_io = Bad_IO(l, M,&g_Rio);
		double searchRatio = 0.25;
		double g_Rc =0.;
		int n_eo = Bad_EO_freeSeg_death(NEIGHBOORHOOD,searchRatio,l,M,&g_Rc);

#ifdef QUALITY_CANDY
		double Echange = beta*exp(-gamma_d*l->dataterm-w_f-w_io*g_Rio-w_eo*g_Rc);
#else
		double Echange = beta*exp(-gamma_d*l->dataterm-w_f-w_io*n_io-w_eo*n_eo);
#endif

		if (Echange < 0.00000001)
			Echange = 0.00000001;


		double R = (M->p_f_b)*n_f/(pow(Echange,T) *(M->p_f_d)* (M->lambda) * (L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN));
		 
		double r = random2();
		if (r < MIN(1,R))
		{
		   pre->next = p->next;
		   M->n_f--;
		   UpdateItsNeighboors_death_freeSeg(l,M);

		   free(l);
		}
	}
}
