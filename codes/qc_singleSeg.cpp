#include "QualityCandy.h"


#if INCLUDE_OPENCV
	#include "gui_functions.h"
#endif


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



/* Adds a Signle Segment and returns the change in energy */
/*******************************************************
Inputs: 													
Candy:   		Linked List with Candy Objectss
img  :   		Original Input Image in Double
lm   :   		Birth Map
img_seg: 		Original Input Image in int
img_height: 	Input Image Heigth
img_width:		Input Image Width
T: 				System Temperature
patch:			unused
patch_len:		unused
Matrix:			unused
prior_pen1:		unused
prior_pen2:		unused
*********************************************************/
double AddSingleSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len, double ***Matrix, double prior_pen1, double prior_pen2)
{
	
	double w_s = M->w_s;


	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_d = M->n_d;
	int n_s = M->n_s;
	int n_f = M->n_f;
	double Echange = 0;

	

	lineObj *SingleSeg = (lineObj *)malloc(sizeof(lineObj));
	SingleSeg->width = (int)floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5);
	SingleSeg->len = (int)floor(random2()*(L_MAX-L_MIN)+L_MIN+0.5);
	SingleSeg->theta = random2()*(THETA_MAX-THETA_MIN)+THETA_MIN;

	
	//Site is a struct x,y
	site  *end = (site *)malloc(2*(n_f+n_s+n_d)*sizeof(site));
	
	int num = count_all_ends(M,end);
	if (num > 0)
	{
		int choose_num = (int)floor(random2()*(num-1)+0.5);
		SingleSeg->enda = end[choose_num];

		SingleSeg->endb = GenerateEndb(SingleSeg->enda,SingleSeg->len,SingleSeg->theta,img_height,img_width);
		SingleSeg->len = (int)sqrt(double(DistSquare(SingleSeg->endb,SingleSeg->enda))); 

		if (SingleSeg->len >= L_MIN)
		{
			SingleSeg->x = (int)floor(0.5*(double(SingleSeg->enda.x+SingleSeg->endb.x))+0.5);
			SingleSeg->y = (int)floor(0.5*(double(SingleSeg->enda.y+SingleSeg->endb.y))+0.5);
			if(lm[SingleSeg->y][SingleSeg->x]==0)
			{
				free(SingleSeg);
				free(end);
				return 0;	
			}

			SingleSeg->endb_L_Num = 0;
			SingleSeg->enda_L_Num = 0;
			SingleSeg->enda_C_Num = 0;
			SingleSeg->endb_C_Num = 0;

			SingleSeg->type = 1;

			double searchRatio = 0.25;
			double g_Rio = 0; 
			double g_Rc = 0.;

			Bad_IO(SingleSeg, M,&g_Rio);
			Bad_EO(NEIGHBOORHOOD,searchRatio,SingleSeg,M,&g_Rc);


			int img_num = 0;
			double iterval = _PI/(RADIUS_SEGS);
			while (SingleSeg->theta > iterval )
			{
				img_num++;
				iterval+=_PI/RADIUS_SEGS;
			}
	
			if (img_num >= int(RADIUS_SEGS))
				img_num--;

			double dterm = dataterm_rec(img_seg, SingleSeg,patch,patch_len,SingleSeg->width,Matrix, prior_pen1,prior_pen2, img_height,img_width);
	
			double Echange_from_neighboor = Echange_from_neighboors_born(M,SingleSeg);

			Echange = -gamma_d*dterm-w_s-w_io*g_Rio-w_eo*g_Rc-Echange_from_neighboor;
			
			double exp_Echange = beta*exp(Echange);
			

			SingleSeg->engergy_for_transition = beta*exp(-gamma_d*dterm-w_s-w_io*g_Rio-w_eo*g_Rc);
			SingleSeg->dataterm = dterm;
			SingleSeg->img_num = img_num;


			double R = pow(exp_Echange,1.0/T) *(M->p_s_d)* (num)* (L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN)/((M->p_s_b)*(n_s+1));

	
			double r = random2();
			if (r < MIN(1,R))
			{
				M->n_s++;
				SingleSeg->type = 1;
		    	LinkedListInsert( M->link_s,M->n_s, SingleSeg);  
				UpdateItsNeighboors_born_SingleSeg(SingleSeg,M);

			   	#if INCLUDE_OPENCV_SINGLE_SEG
					printf("Single Segment Birth Energy: %.2f\n", Echange);		
					display_only_one_double(img, img_height, img_width, SingleSeg, 1);
				#endif
				
				M->Vo += dterm;
				M->VRio += g_Rio;
				M->VReo += g_Rc;
								

			}
			else
			{
				free(SingleSeg);
				Echange = 0;

			}
	}
	else
	{
		free(SingleSeg);
		Echange = 0;
	}
	}
	else
	{
		free(SingleSeg);
		Echange = 0;
	}

	free(end);
	return Echange;
}



/*******************************************************
Function: KillSingleSeg
Kills a Single Segment and returns the change in energy

Inputs: 													
Candy:   		Linked List with Candy Objectss
img  :   		Original Input Image in Double
img_height: 	Input Image Heigth
img_width:		Input Image Width
T: 				System Temperature
*********************************************************/
double KillSingleSeg(Candy *M, double **img, int img_height, int img_width, double T)
{
	double w_s = M->w_s;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_s = M->n_s;
	int n_f = M->n_f;
	double Echange = 0;
	

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
		Bad_IO(l, M,&g_Rio);
		double searchRatio = 0.25;
		double g_Rc = 0.;
		Bad_EO_death(NEIGHBOORHOOD,searchRatio,l,M,&g_Rc);

		double Echange_from_neighboor = Echange_from_neighboors_death(M,l);

		Echange =-(gamma_d*l->dataterm) -w_s - w_io*g_Rio - w_eo*g_Rc + Echange_from_neighboor;
		double exp_Echange = beta*exp(Echange);

		if (exp_Echange < 0.00000001)
			exp_Echange = 0.00000001;

		double R = (M->p_s_b)*n_s/(pow(exp_Echange,1.0/T) *(M->p_s_d) * (2*n_f + n_s)*(L_MAX-L_MIN) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN));
		double r = random2();

		if (r < MIN(1,R))
		{
		   pre->next = p->next;
		   M->n_s--;
		   UpdateItsNeighboors_death_SingleSeg(l,M);

		   	#if INCLUDE_OPENCV_SINGLE_SEG
				printf("Single Segment Death Energy: %.2f\n", -Echange);
				display_only_one_double(img, img_height, img_width, l, 1);
			#endif

			M->Vo += -l->dataterm;
			M->VRio += -g_Rio;
			M->VReo += -g_Rc;


		   free(l);
		   return -Echange;
		}
		return 0;
	}
	return 0;
}
