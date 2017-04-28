#include "QualityCandy.h"


#if INCLUDE_OPENCV
	#include "gui_functions.h"
#endif


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



/* Adds a Free Segment and returns the change in energy */
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
double AddFreeSeg (Candy *M, double **img, double **lm, int **img_seg, int img_height, int img_width, double T,double **patch, int patch_len, double ***Matrix, double prior_pen1, double prior_pen2)
{
	
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	double Echange = 0;


	/* Create Free Seg*/
	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));

	/* Assign Dimensions and angle*/
	int L_MIN_FREE_SEG = 1*L_MIN;
	freeSeg->width = (int)floor(random2()*(W_MAX-W_MIN)+W_MIN+0.5);
	freeSeg->len = (int)floor(random2()*(L_MAX-L_MIN_FREE_SEG)+L_MIN_FREE_SEG+0.5);
	freeSeg->theta = random2()*(THETA_MAX-THETA_MIN)+THETA_MIN;
	
	/* Random Location*/
	freeSeg->enda.x  = (int)floor(random2()*(img_width-2)+0.5);
	freeSeg->enda.y = (int)floor(random2()*(img_height-2)+0.5);

	while (freeSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}

	if (img_num >= int(RADIUS_SEGS))
		img_num--;

	freeSeg->endb = GenerateEndb(freeSeg->enda,freeSeg->len,freeSeg->theta,img_height,img_width);

	freeSeg->len = (int)sqrt(double(DistSquare(freeSeg->endb,freeSeg->enda)));

	if (freeSeg->len >= L_MIN_FREE_SEG)
	{
		freeSeg->x = (int)floor(0.5*(double(freeSeg->enda.x+freeSeg->endb.x))+0.5);
		freeSeg->y = (int)floor(0.5*(double(freeSeg->enda.y+freeSeg->endb.y))+0.5);
		if(lm[freeSeg->y][freeSeg->x]==0)
		{
			free(freeSeg);
			return 0;
		}

		freeSeg->endb_L_Num = 0;
		freeSeg->enda_L_Num = 0;
		freeSeg->enda_C_Num = 0;
		freeSeg->endb_C_Num = 0;
		freeSeg->img_num = img_num;

		double searchRatio = 0.25;
		double g_Rio = 0;
		double g_Rc =0.;
		
		Bad_IO(freeSeg, M,&g_Rio);
		Bad_EO_freeSeg(NEIGHBOORHOOD,searchRatio,freeSeg,M,&g_Rc);


		double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	
		Echange = -gamma_d*dterm - w_f - w_io*g_Rio -w_eo*g_Rc;	
		double exp_Echange = beta*exp(Echange);


		freeSeg->dataterm = dterm;
		freeSeg->engergy_for_transition = exp_Echange;

		double R = pow(exp_Echange,1.0/T) *(M->p_f_d)* (M->lambda) * (L_MAX-L_MIN_FREE_SEG) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN)/( (M->p_f_b)*(n_f+1) );


		double r = random2();
		if (r < MIN(1,R))
		{
			if (freeSeg->enda_C_Num == 0 &&  freeSeg->endb_C_Num == 0)
			{
				M->n_f++;
				freeSeg->type = 0;
		 	   	LinkedListInsert( M->link_f,M->n_f, freeSeg);  //this is a free segment	
				UpdateItsNeighboors_born_freeSeg(freeSeg,M);

				
				#if INCLUDE_OPENCV_FREE_SEG
					printf("Free Segment Birth Energy Change: %.2f\n", Echange);		
					display_only_one_double(img, img_height, img_width, freeSeg, 1);
				#endif

				M->Vo += dterm;
				M->VRio += g_Rio;
				M->VReo += g_Rc;
				return Echange;
			}
			else
			{
				free(freeSeg);
				return 0;
			}
		}
		else
		{
			free(freeSeg);
			return 0;
		}
	}
	else
	{
		free(freeSeg);
		return 0;
	}
}


/* Kills a Free Segment and returns the change in energy */
/*******************************************************
Inputs: 													
Candy:   		Linked List with Candy Objectss
img  :   		Original Input Image in Double
img_height: 	Input Image Heigth
img_width:		Input Image Width
T: 				System Temperature
*********************************************************/
double KillFreeSeg(Candy *M, double **img, int img_height, int img_width, double T)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;
	double Echange = 0;


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
		Bad_IO(l, M,&g_Rio);
		double searchRatio = 0.25;
		double g_Rc =0.;
		Bad_EO_freeSeg_death(NEIGHBOORHOOD,searchRatio,l,M,&g_Rc);

		Echange = -gamma_d*l->dataterm - w_f - w_io*g_Rio - w_eo*g_Rc;
		double exp_Echange = beta*exp(Echange);


		if (exp_Echange < 0.00000001)
			exp_Echange = 0.00000001;

		int L_MIN_FREE_SEG = 1*L_MIN;
		double R = (M->p_f_b)*n_f/( (pow(exp_Echange,1/T) *(M->p_f_d)* (M->lambda) * (L_MAX-L_MIN_FREE_SEG) * (W_MAX-W_MIN) * (THETA_MAX-THETA_MIN)) );
		
		 
		double r = random2();
		if (r < MIN(1,R))
		{
		   pre->next = p->next;
		   M->n_f--;
		   UpdateItsNeighboors_death_freeSeg(l,M);

		     
		   
		   	#if INCLUDE_OPENCV_FREE_SEG
				printf("Free Segment Death Energy Change: %.2f\n", -Echange);		
				display_only_one_double(img, img_height, img_width, l, 1);
			#endif

			free(l);

			M->Vo += -l->dataterm;
			M->VRio += -g_Rio;
			M->VReo += -g_Rc;

			return -Echange;

		}
		return 0.0;
	}

	return 0.0;
}
