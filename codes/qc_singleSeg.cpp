#include "QualityCandy.h"


void UpdateItsNeighboors_born_SingleSeg (lineObj* seg, Candy* M)
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

	if (seg->enda_C_Num != 0)
	{
		for (int i = 0; i < seg->enda_C_Num;i++)
		{
			target = seg->enda_C[i];
			 if (target->type == 0)   //free seg change to single seg 
			 {
				 updateNeighboorLink_born(M,target,seg,seg->enda);

				 LinkedListDelete(M->link_f,target);
				 M->n_f--;
				 M->n_s++;
				 LinkedListInsert(M->link_s,M->n_s,target);
				 target->type = 1;
			 }
			 else if (target->type == 1)
			 {
				updateNeighboorLink_born(M,target,seg,seg->enda);
				if (target->enda_C_Num != 0 && target->endb_C_Num != 0)   //single change to double seg, otherwise, remain in single
				{
					LinkedListDelete(M->link_s,target);
					M->n_s--;
					M->n_d++;
					LinkedListInsert(M->link_d,M->n_d,target);
					target->type = 2;
				}
			 }
			 else   
			 {
				 updateNeighboorLink_born(M,target,seg,seg->enda);  //double seg remain double
			 }
		}
	}
	if (seg->endb_C_Num != 0)
	{
		for (int i = 0; i < seg->endb_C_Num;i++)
		{
			target = seg->endb_C[i];
			 if (target->type == 0)   //free seg change to single seg 
			 {
				 updateNeighboorLink_born(M,target,seg,seg->endb);

				 LinkedListDelete(M->link_f,target);
				 M->n_f--;
				 M->n_s++;
				 LinkedListInsert(M->link_s,M->n_s,target);
				 target->type = 1;
			 }
			 else if (target->type == 1)
			 {
				updateNeighboorLink_born(M,target,seg,seg->endb);
				if (target->enda_C_Num != 0 && target->endb_C_Num != 0)   //single change to double seg, otherwise, remain in single
				{
					LinkedListDelete(M->link_s,target);
					M->n_s--;
					M->n_d++;
					LinkedListInsert(M->link_d,M->n_d,target);
					target->type = 2;
				}
			 }
			 else   
			 {
				 updateNeighboorLink_born(M,target,seg,seg->endb);  //double seg remain double
			 }
		}
	}
}

void UpdateItsNeighboors_death_SingleSeg (lineObj* seg, Candy* M)
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

	if (seg->enda_C_Num != 0)
	{
		for (int i = 0; i < seg->enda_C_Num;i++)
		{
			target = seg->enda_C[i];
			 if (target->type == 2)
			 {
				updateNeighboorLink_death(M,target,seg,seg->enda);
				if (target->enda_C_Num == 0 || target->endb_C_Num == 0)   //double change to single seg, otherwise, remain in double
				{
					LinkedListDelete(M->link_d,target);
					M->n_d--;
					M->n_s++;
					LinkedListInsert(M->link_s,M->n_s,target);
					target->type = 1;
				}
			 }
			 else if(target->type == 1)    
			 {
				 updateNeighboorLink_death(M,target,seg,seg->enda);  
				 if(target->enda_C_Num == 0 && target->endb_C_Num == 0)//single seg change to free seg, otherwise remain single
				 {
					LinkedListDelete(M->link_s,target);
					M->n_s--;
					M->n_f++;
					LinkedListInsert(M->link_f,M->n_f,target);
					target->type = 0;
				 }
			 }
		}
	}
	if (seg->endb_C_Num != 0)
	{
		for (int i = 0; i < seg->endb_C_Num;i++)
		{
			target = seg->endb_C[i];
			 if (target->type == 2)
			 {
				updateNeighboorLink_death(M,target,seg,seg->endb);
				if (target->enda_C_Num == 0 || target->endb_C_Num == 0)   //double change to single seg, otherwise, remain in double
				{
					LinkedListDelete(M->link_d,target);
					M->n_d--;
					M->n_s++;
					LinkedListInsert(M->link_s,M->n_s,target);
					target->type = 1;
				}
			 }
			 else if(target->type == 1)    
			 {
				 updateNeighboorLink_death(M,target,seg,seg->endb);  
				 if(target->enda_C_Num == 0 && target->endb_C_Num == 0)//single seg change to free seg, otherwise remain single
				 {
					LinkedListDelete(M->link_s,target);
					M->n_s--;
					M->n_f++;
					LinkedListInsert(M->link_f,M->n_f,target);
					target->type = 0;
				 }
			 }
		}
	}
}


void AddSingleSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len, double ***Matrix, double prior_pen1, double prior_pen2)
{
	
	double w_s = M->w_s;
	double w_f = M->w_f;
	double w_d = M->w_d;

	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_d = M->n_d;
	int n_s = M->n_s;
	int n_f = M->n_f;

	double lambda = M->lambda;

	

	lineObj *SingleSeg = (lineObj *)malloc(sizeof(lineObj));

	SingleSeg->width = (int)floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5);
	SingleSeg->len = (int)floor(random2()*(L_MAX-L_MIN)+L_MIN+0.5);
	SingleSeg->theta = random2()*(THETA_MAX-THETA_MIN)+THETA_MIN;

	int w_num = SingleSeg->width - W_MIN;
	
	//Site is a struct x,y
	site  *end = (site *)malloc(2*(n_f+n_s+n_d)*sizeof(site));
	#ifdef EN_NEW_SEG_AT_NON_FREE_END

	int num = count_all_ends(M,end);
	if (num > 0){
		int choose_num = (int)floor(random2()*(num-1)+0.5);
		SingleSeg->enda = end[choose_num];
	#else
	int num = 2*n_f+n_s;
	int choose_num = (int)floor(random2()*(num-1)+1+0.5);
	if (num > 0)
	{
		if(choose_num <= 2*n_f)
			SingleSeg->enda = SelectEndfromFreeLink(M->link_f,n_f,(int)floor(double(choose_num+1)/2));
		else
			SingleSeg->enda = SelectEndfromSingleLink(M->link_s,n_s,choose_num-2*n_f);

	if (SingleSeg->enda.x == 0 && SingleSeg->enda.y == 0)
	{
		//printf("selcting end error\n");
	}
	// original quality candy
	//SingleSeg->enda.x += floor((2*random2()-1)*NEIGHBOORHOOD+0.5);
	//SingleSeg->enda.y += floor((2*random2()-1)*NEIGHBOORHOOD+0.5);
#endif
	SingleSeg->endb = GenerateEndb(SingleSeg->enda,SingleSeg->len,SingleSeg->theta,img_height,img_width);

//KDW	SingleSeg->len = sqrt(double(DistSquare(SingleSeg->endb,SingleSeg->enda)));
	SingleSeg->len = (int)sqrt(double(DistSquare(SingleSeg->endb,SingleSeg->enda))); //KDW

	if (SingleSeg->len >= L_MIN)
	{
	SingleSeg->x = (int)floor(0.5*(double(SingleSeg->enda.x+SingleSeg->endb.x))+0.5);
	SingleSeg->y = (int)floor(0.5*(double(SingleSeg->enda.y+SingleSeg->endb.y))+0.5);
	if(lm[SingleSeg->y][SingleSeg->x]==0){
		free(SingleSeg);
		free(end);
		return;
	}

	SingleSeg->endb_L_Num = 0;
	SingleSeg->enda_L_Num = 0;
	SingleSeg->enda_C_Num = 0;
	SingleSeg->endb_C_Num = 0;

	SingleSeg->type = 1;

	double searchRatio = 0.25;
	double g_Rio = 0; 
	double g_Rc = 0.;

	int n_io = Bad_IO(SingleSeg, M,&g_Rio);
	int n_eo = Bad_EO(NEIGHBOORHOOD,searchRatio,SingleSeg,M,&g_Rc);

	//should add data term here, assume homogeneous Poisson 
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (SingleSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;

	//double dterm = dataterm(mid_img[img_num],SingleSeg->enda,SingleSeg->endb);
	double dterm = dataterm_rec(img_seg, SingleSeg,patch,patch_len,SingleSeg->width,Matrix, prior_pen1,prior_pen2, img_height,img_width);
	//double dterm = dataterm(img,SingleSeg->enda,SingleSeg->endb,SingleSeg->theta,SingleSeg->len,img_num,patch,patch_len);
	double Echange_from_neighboor = Echange_from_neighboors_born(M,SingleSeg);

#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_s-w_io*g_Rio-w_eo*g_Rc-Echange_from_neighboor);
	SingleSeg->engergy_for_transition = beta*exp(-gamma_d*dterm-w_s-w_io*g_Rio-w_eo*g_Rc);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_s-w_io*n_io-w_eo*n_eo-Echange_from_neighboor);
	SingleSeg->engergy_for_transition = dterm*beta*exp(-w_s-w_io*n_io-w_eo*n_eo);
#endif
	SingleSeg->dataterm = dterm;
	SingleSeg->img_num = img_num;

#ifdef EN_NEW_SEG_AT_NON_FREE_END
	double R = pow(Echange,T) *(M->p_s_d)* (num)* (L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN)/((M->p_s_b)*(n_s+1));
#else
	double R = pow(Echange,T) *(M->p_s_d)* (2*n_f + n_s)* (L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN)/((M->p_s_b)*(n_s+1));
#endif
	//printf("single_seg_r = %f\n",R);
	double r = random2();
	if (r < MIN(1,R))
	{
		// if ((SingleSeg->enda_C_Num != 0 && SingleSeg->endb_C_Num == 0) ||(SingleSeg->endb_C_Num != 0 && SingleSeg->enda_C_Num ==0))
		//{
//	if(g_Rio>10000)
//		printf("..");
			M->n_s++;
			SingleSeg->type = 1;
		    LinkedListInsert( M->link_s,M->n_s, SingleSeg);  //this is a single segment
			UpdateItsNeighboors_born_SingleSeg(SingleSeg,M);
		//}
		// else
		//	 free(SingleSeg);
	}
	else
		free(SingleSeg);
	}
	else
		free(SingleSeg);
	}
	else
		free(SingleSeg);
	free(end);
}

void KillSingleSeg(Candy *M, double **img, int img_height, int img_width, double T)
{
	double w_s = M->w_s;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_s = M->n_s;
	int n_f = M->n_f;
	double lambda = M->lambda;

	if (n_s != 0)
	{
		int num = (int)floor((n_s-1)*random2()+1+0.5);
		Node *p;
		Node *pre = M->link_s;

		for (int i = 1;i<num;i++)
			pre = pre->next;

		p = pre->next;

		lineObj *l = p->index;
		double g_Rio = 0; 
		int n_io = Bad_IO(l, M,&g_Rio);
		double searchRatio = 0.25;
		double g_Rc = 0.;
		int n_eo = Bad_EO_death(NEIGHBOORHOOD,searchRatio,l,M,&g_Rc);

		double Echange_from_neighboor = Echange_from_neighboors_death(M,l);
#ifdef QUALITY_CANDY
		double Echange = beta*exp(-gamma_d*l->dataterm-w_s-w_io*g_Rio-w_eo*g_Rc+Echange_from_neighboor);
#else
		double Echange = beta*exp(-gamma_d*l->dataterm-w_s-w_io*n_io-w_eo*n_eo+Echange_from_neighboor);
#endif

		if (Echange < 0.00000001)
			Echange = 0.00000001;

		double R = (M->p_s_b)*n_s/(pow(Echange,T) *(M->p_s_d) * (2*n_f + n_s)*(L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN));
		double r = random2();

		if (r < MIN(1,R))
		{
		   pre->next = p->next;
		   M->n_s--;
		   UpdateItsNeighboors_death_SingleSeg(l,M);

		   free(l);
		}
	}
}
