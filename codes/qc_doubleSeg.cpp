#include "QualityCandy.h"

void UpdateItsNeighboors_born_DoubleSeg (lineObj* seg, Candy* M)
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

void UpdateItsNeighboors_death_DoubleSeg (lineObj* seg, Candy* M)
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


site SelectNearestFreeEnd(site a,int a_num, int a_type, Candy *M,int searchR, int *num, int *end_type)
{
	int n_f = M->n_f;
	int n_s = M->n_s;

	double min_dist = searchR*searchR;
	site candidate;
	candidate.x = 0;
	candidate.y = 0;
	Node *p = M->link_f;
	if (a_type == 0) //site a is a freeSeg
	{
	for (int i = 0;i< n_f;i++)
	{
		p = p->next;
		if (DistSquare(a,p->index->enda)<min_dist)
		{
			if (i!=a_num)
			{
				candidate = p->index->enda;
				min_dist = DistSquare(a,p->index->enda);
				*end_type = 0;
				*num = i;
			}
		}
		if (DistSquare(a,p->index->endb)<min_dist)
		{
			if (i!=a_num)
			{
				candidate = p->index->endb;
				min_dist = DistSquare(a,p->index->endb);
				*end_type = 0;
				*num = i;
			}
		}
	}
	p = M->link_s;
		for (int i = 0;i< n_s;i++)
	{
		p = p->next;
		if (p->index->enda_C_Num == 0)
		{
			if (DistSquare(a,p->index->enda)<min_dist)
			{
				candidate = p->index->enda;
				min_dist = DistSquare(a,p->index->enda);
				*end_type = 1;
				*num = i;
			}
		}
		else if (p->index->endb_C_Num == 0)
		{
			if (DistSquare(a,p->index->endb)<min_dist)
			{
				candidate = p->index->endb;
				min_dist = DistSquare(a,p->index->endb);
				*end_type = 1;
				*num = i;
			}
		}
	}
	}
	else
	{	
		for (int i = 0;i< n_f;i++)
	{
		p = p->next;
		if (DistSquare(a,p->index->enda)<min_dist)
		{
				candidate = p->index->enda;
				min_dist = DistSquare(a,p->index->enda);
				*end_type = 0;
				*num = i;
		}
		if (DistSquare(a,p->index->endb)<min_dist)
		{
				candidate = p->index->endb;
				min_dist = DistSquare(a,p->index->endb);
				*end_type = 0;
				*num = i;
		}
	}
	p = M->link_s;
		for (int i = 0;i< n_s;i++)
	{
		p = p->next;
		if (p->index->enda_C_Num == 0)
		{
			if (DistSquare(a,p->index->enda)<min_dist)
			{
				if (i!=a_num)
				{
				candidate = p->index->enda;
				min_dist = DistSquare(a,p->index->enda);
				*end_type = 1;
				*num = i;
				}
			}
		}
		else if (p->index->endb_C_Num == 0)
		{
			if (DistSquare(a,p->index->endb)<min_dist)
			{
				if (i!=a_num)
				{
				candidate = p->index->endb;
				min_dist = DistSquare(a,p->index->endb);
				*end_type = 1;
				*num = i;
				}
			}
		}
		else
		{
//KDW			printf("single seg wrong\n");
		}
	}
	}

	return candidate;
}
int check_exist(site a, site b, lineObj *belong)
{
	int n;
	lineObj *p;

	if((a.x == b.x)&&(a.y == b.y))
		return 1;
	if((a.x == belong->enda.x)&&(a.y == belong->enda.y)){
		if((b.x == belong->endb.x)&&(b.y == belong->endb.y))
			return 1;
		n = belong->enda_C_Num;
		for(int i = 0; i<n; i++){
			p = belong->enda_C[i];
			if((a.x == p->enda.x)&&(a.y == p->enda.y))
			{
				if((b.x == p->endb.x)&&(b.y == p->endb.y))
					return 1;
			}
			else
			{
				if((b.x == p->enda.x)&&(b.y == p->enda.y))
					return 1;
			}
		}
	}
	if((a.x == belong->endb.x)&&(a.y == belong->endb.y)){
		if((b.x == belong->enda.x)&&(b.y == belong->enda.y))
			return 1;
		n = belong->endb_C_Num;
		for(int i = 0; i<n; i++){
			p = belong->endb_C[i];
			if((a.x == p->enda.x)&&(a.y == p->enda.y))
			{
				if((b.x == p->endb.x)&&(b.y == p->endb.y))
					return 1;
			}
			else
			{ 
				if((b.x == p->enda.x)&&(b.y == p->enda.y))
					return 1;
			}
		}
	}
	return 0;
}

int SelectNearestEnd_in_link(LinkedList link, int n, site a, lineObj *belong, double *min_dist, site *candid)
{
	Node *p = link;
	int found = 0;


	for (int i = 0;i< n;i++)
	{
		p = p->next;
		if (DistSquare(a,p->index->enda)<*min_dist)
		{
			if (!check_exist(a, p->index->enda, belong))
			{
				found = 1;
				(*candid) = p->index->enda;
				*min_dist = DistSquare(a,p->index->enda);
			}
		}
		if (DistSquare(a,p->index->endb)<*min_dist)
		{
			if (!check_exist(a, p->index->endb, belong))
			{
				found = 1;
				(*candid) = p->index->endb;
				*min_dist = DistSquare(a,p->index->endb);
			}
		}
	}
	return found;
}

int SelectNearestEnd(site a,lineObj *belong,Candy *M,int searchR, site *candid)
{
	int found=0;
	double min_dist = searchR*searchR;

	found |= SelectNearestEnd_in_link(M->link_f,M->n_f,a,belong,&min_dist,candid);
	found |= SelectNearestEnd_in_link(M->link_s,M->n_s,a,belong,&min_dist,candid);
	found |= SelectNearestEnd_in_link(M->link_d,M->n_d,a,belong,&min_dist,candid);

	return found;

}


void AdddoubleSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l, double ****img_seg_l,int img_height, int img_width, double T,double **patch, int patch_len , double ***Matrix, double prior_pen1, double prior_pen2, int test)
{
	
	double w_d = M->w_d;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_s = M->n_s;
	int n_f = M->n_f;
	int n_d = M->n_d;


	

	lineObj *DoubleSeg = (lineObj *)malloc(sizeof(lineObj));

	//DoubleSeg->width = floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5);
	//DoubleSeg->len = floor(random2()*(L_MAX-L_MIN)+L_MIN+0.5);
	//DoubleSeg->theta = random2()*(THETA_MAX-THETA_MIN)+THETA_MIN;

	site  *end = (site *)malloc(2*(n_f+n_s+n_d)*sizeof(site));
	lineObj **obj = (lineObj **)malloc(2*(n_f+n_s+n_d)*sizeof(lineObj*));

#ifdef EN_NEW_SEG_AT_NON_FREE_END
	int num = count_all_ends(M,end,obj);
	if (num > 0){
		int choose_num = (int)floor(random2()*(num-1)+0.5);
		DoubleSeg->enda = end[choose_num];
//		site tmp;
//		tmp.x = 46; tmp.y = 59;
//		if(DistSquare(DoubleSeg->enda,tmp)<5)
//			printf("..");
		lineObj *belong = obj[choose_num];
		site candidate;
		if(SelectNearestEnd(DoubleSeg->enda,belong,M,L_MAX,&candidate))
			DoubleSeg->endb = candidate;
		else{
			free(end);
			free(obj);
			return;
		}
#else
	int enda_type, endb_type; //0:freeSeg 1:singleSeg
	int enda_num,endb_num;
	int num = 2*n_f+n_s;
	if (n_f+n_s > 1)
	{
		int choose_num = (int)floor(random2()*(num-1)+1+0.5);
		if(choose_num <= 2*n_f)
			{
				DoubleSeg->enda = SelectEndfromFreeLink(M->link_f,n_f,(int)floor(double(choose_num+1)/2));
				enda_type = 0;
				enda_num = (int)floor(double(choose_num+1)/2)-1;
		}
		else
			{
				DoubleSeg->enda = SelectEndfromSingleLink(M->link_s,n_s,choose_num-2*n_f);
				enda_type = 1;
				enda_num = choose_num-2*n_f-1;
		}

		if (DoubleSeg->enda.x == 0 && DoubleSeg->enda.y == 0)
		{
			//printf("selcting end error double\n");
		}
		DoubleSeg->endb = SelectNearestFreeEnd(DoubleSeg->enda,enda_num,enda_type,M,L_MAX,&endb_num,&endb_type);
#endif
	//original quality candy
	//DoubleSeg->enda.x += floor((2*random2()-1)*NEIGHBOORHOOD+0.5);
	//DoubleSeg->enda.y += floor((2*random2()-1)*NEIGHBOORHOOD+0.5);
	//DoubleSeg->endb.x += floor((2*random2()-1)*NEIGHBOORHOOD+0.5);
	//DoubleSeg->endb.y += floor((2*random2()-1)*NEIGHBOORHOOD+0.5);

	//	int choose_num2 = floor(random2()*(num-1)+1+0.5);
	//	while (choose_num2 == choose_num)
	//	{
	//		choose_num2 = floor(random2()*(num-1)+1+0.5);
	//	}
	//	if(choose_num2 <= 2*n_f)
	//		DoubleSeg->endb = SelectEndfromFreeLink(M->link_f,n_f,floor(double(choose_num2+1)/2));
	//	else
	//		DoubleSeg->endb = SelectEndfromSingleLink(M->link_s,n_s,choose_num2-2*n_f);

	DoubleSeg->len = (int)sqrt(double(DistSquare(DoubleSeg->endb,DoubleSeg->enda)));

	if (DoubleSeg->endb.x == 0 && DoubleSeg->endb.y == 0)
		DoubleSeg->len  = (int)L_MAX+1;

	//printf("DoubleSeg->len = %d\n",DoubleSeg->len);
	if (DoubleSeg->len >= L_MIN && DoubleSeg->len <= L_MAX)
	//if ( DoubleSeg->len <= L_MAX)
	{
	DoubleSeg->x = (int)floor(0.5*(double(DoubleSeg->enda.x+DoubleSeg->endb.x))+0.5);
	DoubleSeg->y = (int)floor(0.5*(double(DoubleSeg->enda.y+DoubleSeg->endb.y))+0.5);
	if(lm[DoubleSeg->y][DoubleSeg->x]==0){
		free(DoubleSeg);
		free(end);
		free(obj);
		return;
	}

	double theta;

	if(DoubleSeg->enda.y < DoubleSeg->endb.y)
	{
		if (DoubleSeg->enda.x-DoubleSeg->endb.x != 0)
		 theta = atan(double((DoubleSeg->endb.y-DoubleSeg->enda.y)/double(DoubleSeg->enda.x-DoubleSeg->endb.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
		DoubleSeg->theta = theta;
	}
	else
	{
		if (DoubleSeg->enda.x-DoubleSeg->endb.x != 0)
		 theta = atan(double((DoubleSeg->enda.y-DoubleSeg->endb.y)/double(DoubleSeg->endb.x-DoubleSeg->enda.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
		DoubleSeg->theta = theta;
	}

//KDW	DoubleSeg->width = floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5);
	DoubleSeg->width = (int)floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5); //KDW
	DoubleSeg->endb_L_Num = 0;
	DoubleSeg->enda_L_Num = 0;
	DoubleSeg->enda_C_Num = 0;
	DoubleSeg->endb_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0;

	Bad_IO(DoubleSeg, M ,&g_Rio);



	double g_Rc =0.;
	Bad_EO(NEIGHBOORHOOD,searchRatio,DoubleSeg,M,&g_Rc);

	//should add data term here, assume homogeneous Poisson 
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (DoubleSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}
//	if (random2()<0.5)
//		img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;

	//double dterm = dataterm_double(mid_img[img_num],DoubleSeg->enda,DoubleSeg->endb);
	double dterm =  dataterm_rec(img_seg, DoubleSeg,patch,patch_len,DoubleSeg->width, Matrix, prior_pen1,prior_pen2, img_height,img_width);
	//double dterm = dataterm(img,DoubleSeg->enda,DoubleSeg->endb,DoubleSeg->theta,DoubleSeg->len,img_num,patch,patch_len);
	
	double Echange_from_neighboor = Echange_from_neighboors_born(M,DoubleSeg);
	
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_d-w_io*g_Rio-w_eo*g_Rc-Echange_from_neighboor);
	DoubleSeg->engergy_for_transition = beta*exp(-gamma_d*dterm-w_d-w_io*g_Rio-w_eo*g_Rc);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_d-w_io*n_io-w_eo*n_eo-Echange_from_neighboor);
	DoubleSeg->engergy_for_transition = beta*exp(-gamma_d*dterm-w_d-w_io*n_io-w_eo*n_eo);
#endif
	DoubleSeg->dataterm = dterm;
	DoubleSeg->img_num = img_num;
	
	double R = pow(Echange,T) *(M->p_d_d)* (2*n_f + n_s)/** (2*n_f + n_s-1)*/ /((M->p_d_b)*(n_d+1));
	
	//printf("double_seg_r = %f,n_io = %d, n_eo = %d\n",R,n_io,n_eo);
	double r = random2();
	if (r < MIN(1,R))
	{
		//give birth to such segment
		if (DoubleSeg->enda_C_Num * DoubleSeg->endb_C_Num != 100)
		{
//	if(g_Rio>10000)
//		printf("..");
			M->n_d++;
			DoubleSeg->type = 2;
		    LinkedListInsert( M->link_d,M->n_d, DoubleSeg);  //this is a double segment
			UpdateItsNeighboors_born_DoubleSeg(DoubleSeg,M);
		}
		else
		{
			printf("add double seg wrong\n");	
			free(DoubleSeg);
		}
	}
	else
		free(DoubleSeg);
	}
	else
		free(DoubleSeg);
	}
	else
		free(DoubleSeg);
	free(end);
	free(obj);
}

void KillDoubleSeg(Candy *M, double **img, int img_height, int img_width, double T)
{
	double w_d = M->w_d;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_s = M->n_s;
	int n_f = M->n_f;
	int n_d = M->n_d;

	if (n_d != 0)
	{
		int num = (int)floor((n_d-1)*random2()+1+0.5);
		Node *p;
		Node *pre = M->link_d;

		for (int i = 1;i<num;i++)
			pre = pre->next;

		p = pre->next;

		lineObj *l = p->index;
		double g_Rio = 0; 
		Bad_IO(l, M,&g_Rio);
		double searchRatio = 0.25;
		double g_Rc =0.;
		Bad_EO_death(NEIGHBOORHOOD,searchRatio,l,M,&g_Rc);

		double Echange_from_neighboor = Echange_from_neighboors_death(M,l);
#ifdef QUALITY_CANDY
		double Echange = beta*exp(-gamma_d*l->dataterm-w_d-w_io*g_Rio-w_eo*g_Rc+Echange_from_neighboor);
#else
		double Echange = beta*exp(-gamma_d*l->dataterm-w_d-w_io*n_io-w_eo*n_eo+Echange_from_neighboor);
#endif

		if (Echange < 0.00000001)
			Echange = 0.00000001;

		double R = (M->p_d_b)*n_d/(pow(Echange,T) *(M->p_d_d) * (2*n_f + n_s));//* (2*n_f + n_s-1));
		double r = random2();

		if (r < MIN(1,R))
		{
		   pre->next = p->next;
		   M->n_d--;
		   UpdateItsNeighboors_death_DoubleSeg(l,M);
		 
		   free(l);
		}
	}
}
