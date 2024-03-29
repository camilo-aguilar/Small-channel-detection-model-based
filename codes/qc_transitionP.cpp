#include "QualityCandy.h"
#include "tiff.h"
#include "Parameters.h"


#if INCLUDE_OPENCV
	#include "gui_functions.h"	
#endif


lineObj* select_object_from_link(LinkedList link,int num, int choose_num)
{
	Node *p = link;
	for (int i = 0;i< choose_num;i++)
		p = p->next;
	return (p->index);
}

double FreeSeg_length_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;

	double g_Rio = 0; 
	double g_Rio_before = 0;
	double g_Rc =0.;
	double g_Rc_before =0.;

	double searchRatio = 0.25;

	if (n_f == 0)
		return 0.0;
	
	int choose_num = (int)floor(random2()*(n_f-1)+1+0.5);

	Node *pre = M->link_f;
	for (int i = 1;i<choose_num;i++)
	{
		pre = pre->next;
	}
	Node *p = pre->next;

	//CGA: Chose Object 
	lineObj *object = p->index;
	//CGA: Create New Object 
	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
	
	//CGA: d_len = L + [-delta_L delta_L]
	int d_len = (int)floor(random2()*2*DELTA_LEN - DELTA_LEN +0.5);
	double new_len = d_len + object->len;

	if(new_len < L_MIN || new_len > L_MAX)
	{
		free(freeSeg);
		return 0;
	}

	double theta = object->theta;
	int c_x = object->x;
	int c_y = object->y;

	//CGA: Create New End Points A and B
	site new_enda, new_endb;
	new_enda.x = int(c_x+new_len/2.0*cos(theta));
	new_enda.y = int(c_y-new_len/2.0*sin(theta));
	if (new_enda.x < 0 ||new_enda.x >= img_width ||new_enda.y < 0 ||new_enda.y >= img_height)
	{
			free(freeSeg);
			return 0;
	}

	new_endb.x = int(c_x-new_len/2.0*cos(theta));
	new_endb.y = int(c_y+new_len/2.0*sin(theta));
	if (new_endb.x < 0 ||new_endb.x >= img_width ||new_endb.y < 0 ||new_endb.y >= img_height)
	{
		free(freeSeg);
		return 0;
	}

	//CGA: Create New Free Segment and its Properties
	freeSeg->enda = new_enda;
	freeSeg->endb = new_endb;
	freeSeg->x = (int)floor(0.5*(double(freeSeg->enda.x+freeSeg->endb.x))+0.5);
	freeSeg->y = (int)floor(0.5*(double(freeSeg->enda.y+freeSeg->endb.y))+0.5);
	freeSeg->len = (int)new_len;
	freeSeg->theta = theta;
	freeSeg->width = object->width;
	freeSeg->img_num = object->img_num;
	freeSeg->endb_L_Num = 0;
	freeSeg->enda_L_Num = 0;
	freeSeg->enda_C_Num = 0;
	freeSeg->endb_C_Num = 0;

	
	Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	
	Bad_IO(object, M ,&g_Rio_before);
	Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);


	if (freeSeg->enda_C_Num != 0 || freeSeg->endb_C_Num != 0)
	{
		free(freeSeg);
		return 0;
	}

	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);

	double Echange = -gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc;
	double exp_Echange = beta*exp(Echange);
	
	double Echange_before = -gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before;
	double exp_Echange_before = beta*exp(Echange_before);

	double R = pow(exp_Echange/exp_Echange_before,1.0/T);
	double r = random2();
	if (r < MIN(1,R) && (Echange- Echange_before))
	{
		pre->next = p->next;  // kill the object, add the freeSeg
	    M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);

		M->n_f++;
		freeSeg->dataterm = dterm;
		freeSeg->engergy_for_transition = Echange;
		freeSeg->type = 0;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);

		

		M->Vo += (dterm - object->dataterm);
		M->VRio += (g_Rio- g_Rio_before);
		M->VReo += (g_Rc - g_Rc_before);

		#if INCLUDE_OPENCV_TRANSITION
			printf("Length Move Energy: %.2f\n", Echange- Echange_before);		
			display_only_one_double(img, img_height, img_width, freeSeg, 1);
		#endif


		return (Echange - Echange_before);
	}
	else
	{
		free(freeSeg);
		Echange = 0;
	}
	return 0;
}

double FreeSeg_width_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;

	if (n_f == 0)
		return 0;
	
	int choose_num = (int)floor(random2()*(n_f-1)+1+0.5);

	Node *pre = M->link_f;
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;
	Node *p = pre->next;

	lineObj *object = p->index;


	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
	int d_width = (int)floor(random2()*2*DELTA_WIDTH - DELTA_WIDTH +0.5);
	double new_width = d_width + object->width;

	if(new_width < W_MIN || new_width > W_MAX)
	{
			free(freeSeg);
			return 0;
	}

	double theta = object->theta;


	freeSeg->enda = object->enda;
	freeSeg->endb = object->endb;
	freeSeg->x = object->x;
	freeSeg->y = object->y;
	freeSeg->len = object->len;
	freeSeg->theta = theta;
	freeSeg->width = (int)new_width;
	freeSeg->img_num = object->img_num;
	freeSeg->endb_L_Num = 0;
	freeSeg->enda_L_Num = 0;
	freeSeg->enda_C_Num = 0;
	freeSeg->endb_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0; 
	double g_Rc =0.;
	Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio); // we don't need this operation
	Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);// we don't need this operation

	if (freeSeg->enda_C_Num != 0 || freeSeg->endb_C_Num != 0)
	{
		free(freeSeg);
		return 0;
	}

	
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	
	double Echange = -gamma_d*dterm - w_f - w_io*g_Rio - w_eo*g_Rc;
	double Echange_before = -gamma_d*object->dataterm - w_f - w_io*g_Rio - w_eo*g_Rc;

	double exp_Echange = beta*exp(Echange);
	double exp_Echange_before = beta*exp(Echange_before);

	double R = pow(exp_Echange/exp_Echange_before,1.0/T);
	double r = random2();
	if (r < MIN(1,R) && (Echange- Echange_before))
	{
		pre->next = p->next;  // kill the object, add the freeSeg
	    M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);

		M->n_f++;
		freeSeg->dataterm = dterm;
		freeSeg->engergy_for_transition = Echange;
		freeSeg->type = 0;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);

		M->Vo += (dterm - object->dataterm);
		M->VRio += (g_Rio- g_Rio);
		M->VReo += (g_Rc - g_Rc);


		#if INCLUDE_OPENCV_TRANSITION
			printf("Free Segment Width Move Energy: %.2f\n", Echange-Echange_before);		
			display_only_one_double(img, img_height, img_width, freeSeg, 1);
		#endif

		return (Echange - Echange_before);
	}
	else
	{
		free(freeSeg);
		return 0;
	}
	
}

double FreeSeg_theta_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;

	if (n_f == 0)
		return 0;
	

	int choose_num = (int)floor(random2()*(n_f-1)+1+0.5);	
	Node *pre = M->link_f;
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;
	Node *p = pre->next;

	lineObj *object = p->index;

	double d_theta = random2()*2*DELTA_THETA - DELTA_THETA;
	double new_theta = d_theta + object->theta;

	if(new_theta < THETA_MIN || new_theta > THETA_MAX)
		return 0;

	int c_x = object->x;
	int c_y = object->y;
	double len = double(object->len);

	site new_enda, new_endb;
	new_enda.x = int(c_x+len/2.0*cos(new_theta));
	new_enda.y = int(c_y-len/2.0*sin(new_theta));
	if (new_enda.x < 0 ||new_enda.x >= img_width ||new_enda.y < 0 ||new_enda.y >= img_height)
		return 0;
	new_endb.x = int(c_x-len/2.0*cos(new_theta));
	new_endb.y = int(c_y+len/2.0*sin(new_theta));
	if (new_endb.x < 0 ||new_endb.x >= img_width ||new_endb.y < 0 ||new_endb.y >= img_height)
		return 0;
	

	double searchRatio = 0.25;

	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (new_theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}

	if (img_num >= int(RADIUS_SEGS))
		img_num--;

	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
	freeSeg->enda = new_enda;
	freeSeg->endb = new_endb;
	freeSeg->x = (int)floor(0.5*(double(freeSeg->enda.x+freeSeg->endb.x))+0.5);
	freeSeg->y = (int)floor(0.5*(double(freeSeg->enda.y+freeSeg->endb.y))+0.5);
	freeSeg->len = (int)len;
	freeSeg->width = object->width;
	freeSeg->theta = new_theta;
	freeSeg->img_num = img_num;
	freeSeg->endb_L_Num = 0;
	freeSeg->enda_L_Num = 0;
	freeSeg->enda_C_Num = 0;
	freeSeg->endb_C_Num = 0;
 

	if((freeSeg->x == 62)&&(freeSeg->y==138))
		printf("..");
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;
	
	
	Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	Bad_IO(object, M,&g_Rio_before);
	Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

	if (freeSeg->enda_C_Num != 0 || freeSeg->endb_C_Num != 0)
	{
		free(freeSeg);
		return 0;
	}
	  

	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	

	double Echange =        -gamma_d*dterm            -w_f - w_io*g_Rio         -w_eo*g_Rc;
	double Echange_before = -gamma_d*object->dataterm -w_f - w_io*g_Rio_before  -w_eo*g_Rc_before;
	
	double exp_Echange = beta*exp(Echange);
	double exp_Echange_before = beta*exp(Echange_before);

	double R = pow(exp_Echange/exp_Echange_before,1.0/T);

	double r = random2();
	if (r < MIN(1,R))
	{

		pre->next = p->next;  // kill the object, add the freeSeg
	    M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);
		M->n_f++;
		freeSeg->dataterm = dterm;
		freeSeg->engergy_for_transition = Echange;
		freeSeg->type = 0;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);

		M->Vo += (dterm - object->dataterm);
		M->VRio += (g_Rio- g_Rio_before);
		M->VReo += (g_Rc - g_Rc_before);

		#if INCLUDE_OPENCV_TRANSITION
			printf("Free Segment Theta Move Energy: %.2f\n", Echange-Echange_before);		
			display_only_one_double(img, img_height, img_width, freeSeg, 1);
		#endif

		return ( (double)(Echange - Echange_before) );
	}
	else
	{
		free(freeSeg);
		return 0;
	}
	
}

double FreeSeg_freeEnd_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;

	if (n_f == 0)
		return 0.0;
	
	int choose_num = (int)floor(random2()*(n_f-1)+1+0.5);
	
	Node *pre = M->link_f;
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;
	Node *p = pre->next;

	lineObj *object = p->index;

	site end_candidate;
	site end_noncandidate;

	double rr = random2();

	if (rr<0.5)
	{
		end_candidate = object->enda;
		end_noncandidate = object->endb;
	}
	else
	{
		end_candidate = object->endb;
		end_noncandidate = object->enda;
	}

	int delta_x = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
	int delta_y = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);

	end_candidate.x += delta_x;
	end_candidate.y += delta_y;

	if(end_candidate.x <0 || end_candidate.x >= img_width ||end_candidate.y < 0 || end_candidate.y >= img_height)
		return 0;

	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
	freeSeg->enda = end_candidate;
	freeSeg->endb = end_noncandidate;
	freeSeg->x = (int)floor(0.5*(double(end_candidate.x+end_noncandidate.x))+0.5);
	freeSeg->y = (int)floor(0.5*(double(end_candidate.y+end_noncandidate.y))+0.5);
	freeSeg->theta= theta_from_two_ends(freeSeg->enda, freeSeg->endb);
	freeSeg->width = object->width;
	freeSeg->len = (int)sqrt(double(DistSquare(freeSeg->endb,freeSeg->enda)));
	if (freeSeg->len < L_MIN || freeSeg->len > L_MAX)
	{
			free(freeSeg);
			return 0;
	}
	
	freeSeg->endb_L_Num = 0;
	freeSeg->enda_L_Num = 0;
	freeSeg->endb_C_Num = 0;
	freeSeg->enda_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;
	Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	Bad_IO(object, M,&g_Rio_before);
	Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

	if ((freeSeg->endb_C_Num !=0)||(freeSeg->enda_C_Num != 0))
	{
		free(freeSeg);
		return 0;
	}
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (freeSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}


	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	

	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
	
	double Echange = -gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc;
	double Echange_before = -gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before;

	double exp_Echange = beta*exp(Echange);
	double exp_Echange_before = beta*exp(Echange_before);


	freeSeg->dataterm = dterm;
	freeSeg->img_num = img_num;
	freeSeg->engergy_for_transition = Echange;
	freeSeg->type = 0;

	double R = pow(exp_Echange/exp_Echange_before,1.0/T);
	
	double r = random2();
	if (r < MIN(1,R) && (Echange-Echange_before))
	{
		pre->next = p->next;  // kill the object, add the freeSeg
		M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);
		M->n_f++;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);

		M->Vo += (dterm - object->dataterm);
		M->VRio += (g_Rio- g_Rio_before);
		M->VReo += (g_Rc - g_Rc_before);

		#if INCLUDE_OPENCV_TRANSITION
			printf("Free Segment Theta Move Energy: %.2f\n", Echange-Echange_before);		
			display_only_one_double(img, img_height, img_width, freeSeg, 1);
		#endif


		return (Echange - Echange_before);
	}
	else
	{
			free(freeSeg);
			return 0.0;
	}

}

double FreeSeg_center_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;

	if (n_f == 0)
		return 0.0;
	
	int choose_num = (int)floor(random2()*(n_f-1)+1+0.5);
	Node *pre = M->link_f;
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;
	Node *p = pre->next;

	lineObj *object = p->index;

	int delta_x = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
	int delta_y = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);

	site new_enda, new_endb;
	new_enda.x = object->enda.x+delta_x;
	new_enda.y = object->enda.y+delta_y;
	if (new_enda.x < 0 ||new_enda.x >= img_width ||new_enda.y < 0 ||new_enda.y >= img_height)
		return 0.0;
	new_endb.x = object->endb.x+delta_x;
	new_endb.y = object->endb.y+delta_y;
	if (new_endb.x < 0 ||new_endb.x >= img_width ||new_endb.y < 0 ||new_endb.y >= img_height)
		return 0.0;
	

	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
	freeSeg->enda = new_enda;
	freeSeg->endb = new_endb;
	freeSeg->x = object->x+delta_x;
	freeSeg->y = object->y+delta_y;
	freeSeg->theta= object->theta;
	freeSeg->width = object->width;
	freeSeg->len = object->len;
	
	freeSeg->endb_L_Num = 0;
	freeSeg->enda_L_Num = 0;
	freeSeg->endb_C_Num = 0;
	freeSeg->enda_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;

	Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	Bad_IO(object, M,&g_Rio_before);
	Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

	if ((freeSeg->endb_C_Num !=0)||(freeSeg->enda_C_Num != 0))
	{
		free(freeSeg);
		return 0;
	}
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (freeSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);

	double Echange = -gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc;
	double Echange_before = -gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before;
	
	double exp_Echange = beta*exp(Echange);
	double exp_Echange_before = beta*exp(Echange_before);

	freeSeg->dataterm = dterm;
	freeSeg->img_num = img_num;
	freeSeg->engergy_for_transition = Echange;
	freeSeg->type = 0;

	double R = pow(exp_Echange/exp_Echange_before,1.0/T);
	
	double r = random2();
	if (r < MIN(1,R) && (Echange-Echange_before))
	{
		pre->next = p->next;  // kill the object, add the freeSeg
		M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);
		M->n_f++;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);

		
		M->Vo += (dterm - object->dataterm);
		M->VRio += (g_Rio- g_Rio_before);
		M->VReo += (g_Rc - g_Rc_before);
		
		#if INCLUDE_OPENCV_TRANSITION
			printf("Free Segment Free End Move Energy: %.2f\n", Echange-Echange_before);		
			display_only_one_double(img, img_height, img_width, freeSeg, 1);
		#endif

		return (Echange - Echange_before);
	}
	else
	{
			free(freeSeg);
			return 0;
	}
	return 0;

}


double SingleSeg_freeEnd_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_s = M->w_s;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_s = M->n_s;

	if (n_s == 0)
		return 0;
	
	int choose_num = (int)floor(random2()*(n_s-1)+1+0.5);
	//lineObj* object = select_object_from_link(M->link_s,n_s,choose_num);

	Node *pre = M->link_s;
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;
	Node *p = pre->next;

	lineObj *object = p->index;

	site end_candidate;
	site end_connected;

	if (object->enda_C_Num ==0 && object->endb_C_Num != 0)
		{
			end_candidate = object->enda;
			end_connected = object->endb;
	}
	else if (object->enda_C_Num !=0 && object->endb_C_Num == 0)
		{
			end_candidate = object->endb;
			end_connected = object->enda;
	}
	else
	{
		return 0;
	}

	int delta_x = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
	int delta_y = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);

	end_candidate.x += delta_x;
	end_candidate.y += delta_y;

	if(end_candidate.x <0 || end_candidate.x >= img_width ||end_candidate.y < 0 || end_candidate.y >= img_height)
		return 0;

	lineObj *SingleSeg = (lineObj *)malloc(sizeof(lineObj));
	SingleSeg->enda = end_candidate;
	SingleSeg->endb = end_connected;
	SingleSeg->x = (int)floor(0.5*(double(end_candidate.x+end_connected.x))+0.5);
	SingleSeg->y = (int)floor(0.5*(double(end_candidate.y+end_connected.y))+0.5);
	SingleSeg->theta= theta_from_two_ends(SingleSeg->enda, SingleSeg->endb);
	SingleSeg->width = object->width;
	SingleSeg->len = (int)sqrt(double(DistSquare(SingleSeg->endb,SingleSeg->enda)));
	if (SingleSeg->len < L_MIN || SingleSeg->len > L_MAX)
		{
			free(SingleSeg);
			return 0;
	}
	
	SingleSeg->endb_L_Num = 0;
	SingleSeg->enda_L_Num = 0;
	SingleSeg->endb_C_Num = 0;
	SingleSeg->enda_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;

	Bad_IO_transition(SingleSeg, M,choose_num-1,1,&g_Rio);
	Bad_EO_transition(NEIGHBOORHOOD,searchRatio,SingleSeg,M,choose_num-1,1,&g_Rc);
	Bad_IO(object, M,&g_Rio_before);
	Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

	if(fabs(g_Rc) == INF || fabs(g_Rc_before) == INF || fabs(g_Rio) == INF || fabs(g_Rio_before) == INF)
	{
		free(SingleSeg);
		return 0;
	}
	if (SingleSeg->endb_C_Num * SingleSeg->enda_C_Num != 0)
	{
		free(SingleSeg);
		return 0;
	}
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (SingleSeg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	

	double dterm = dataterm_rec(img_seg, SingleSeg,patch,patch_len,SingleSeg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);

	double Echange = -gamma_d*dterm-w_s-w_io*g_Rio-w_eo*g_Rc;
	double exp_Echange = beta*exp(Echange);
	
	double Echange_before = -gamma_d*object->dataterm-w_s-w_io*g_Rio_before-w_eo*g_Rc_before;
	double exp_Echange_before = beta*exp(Echange_before);

	SingleSeg->dataterm = dterm;
	SingleSeg->img_num = img_num;
	SingleSeg->engergy_for_transition = Echange;
	SingleSeg->type = 1;

	double R = pow(exp_Echange/exp_Echange_before,1.0/T);

	double r = random2();
	if (r < MIN(1,R)&& (Echange-Echange_before))
	{
		pre->next = p->next;  // kill the object, add the freeSeg
		M->n_s--;
		UpdateItsNeighboors_death_SingleSeg(object,M);
		free(object);
		M->n_s++;
		LinkedListInsert( M->link_s,M->n_s, SingleSeg);  	
		UpdateItsNeighboors_born_SingleSeg(SingleSeg,M);

		M->Vo += (dterm - object->dataterm);
		M->VRio += (g_Rio- g_Rio_before);
		M->VReo += (g_Rc - g_Rc_before);

		#if INCLUDE_OPENCV_TRANSITION
			printf("Single Segment Free End Move Energy: %.2f\n", Echange-Echange_before);		
			display_only_one_double(img, img_height, img_width, SingleSeg, 1);
		#endif

		return (Echange - Echange_before);
	}
	else
	{
			free(SingleSeg);
			return 0;
	}
	

}

double SingleDoubleSeg_Connection_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	//CGA: Get Parameters
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;


	int C_n = M->connection_n;
	int error = 1;
	double g_Rio = 0;
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;

	double return_energy = 0;

	double searchRatio = 0.25;
	double theta;
	int n_io = 0, n_eo = 0, n_io_before = 0, n_eo_before = 0;
	double dtmp, dterm = 0., dterm_before = 0.;
	double Echange = 0, Echange_before = 0., exp_Echange, exp_Echange_before;

	double total_gRc = 0, total_gRio = 0;
	double total_gRc_before = 0, total_gRio_before = 0;

	

	if(C_n == 0)
		return 0;

	// CGA: Choose any connection
	int choose_num = (int)floor(random2()*(C_n-1)+1+0.5);
	NClinks *pre = M->connection_l;

	//CGA: Choose the actual connection
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;

	//CGA: *p is the segment, ob1 and ob2 are its connections 
	NClinks *connected_segment_to_be_moved = pre->next;
	lineObj *first_connected_end = connected_segment_to_be_moved->object1;
	lineObj *second_connected_end = connected_segment_to_be_moved->object2;

	site end1,end2;

	//CGA find which end of each object is connected/close to the other object
	NC_pairs(first_connected_end,second_connected_end,&end1,&end2);  

	site new_endc;
	int x_offset;
	int y_offset;
	int n=0;
	
	//CGA: Find Obset to Move
	do{
		x_offset = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
		y_offset = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
		n++;
		if(n>1000) return 0;
	}while((x_offset==0)&&(y_offset==0));

	//CGA end1 is close to end2, so new end will be close as well
	new_endc.x = end1.x + x_offset;
	new_endc.y = end1.y + y_offset;

	if(new_endc.x <0 || new_endc.x >= img_width || new_endc.y <0 || new_endc.y >= img_height)
		return 0;

	//CGA: Create new object and new end
	lineObj **old_array_of_neighbors;
	int C_Num, *endc; // endc array: connection point 0:enda, 1:endb
	
	//If end1 was part of enda object1
	if ((first_connected_end->enda.x == end1.x)&&(first_connected_end->enda.y==end1.y))
	{
		//Add +1 to Connection Number
		C_Num = first_connected_end->enda_C_Num+1;
		
		//CGA: Create New Array of neighbors for that end
		old_array_of_neighbors = (lineObj **)malloc(C_Num*sizeof(lineObj *)); // first_connected_end and its connections(C_Num)
		
		//CGA: Array to keep track if its enda or endb
		endc = (int *)malloc(C_Num*sizeof(int));
		
		for(int i=0;i<C_Num-1;i++)
		{ 
			//CGA: Copy List of Neighbors attached to enda_C (connected to enda)
			old_array_of_neighbors[i] = first_connected_end->enda_C[i];
			
			//CGA: See if object1 is connected in enda or endb of each neighbor obj[i]
			NC_pairs(first_connected_end,old_array_of_neighbors[i],&end1,&end2);			
			if ((old_array_of_neighbors[i]->enda.x == end2.x)&&(old_array_of_neighbors[i]->enda.y==end2.y))
				endc[i] = 0; // enda
			else
				endc[i] = 1; // endb
		}
		endc[C_Num-1] = 0; // enda
	}

	//If end1 was part of endb of object1
	else
	{
		C_Num = first_connected_end->endb_C_Num+1;
		old_array_of_neighbors = (lineObj **)malloc(C_Num*sizeof(lineObj *)); // first_connected_end and its connections(C_Num)
		endc = (int *)malloc(C_Num*sizeof(int));
		for(int i=0;i<C_Num-1;i++)
		{ 
			old_array_of_neighbors[i] = first_connected_end->endb_C[i];
			NC_pairs(first_connected_end,old_array_of_neighbors[i],&end1,&end2);
			if ((old_array_of_neighbors[i]->enda.x == end2.x)&&(old_array_of_neighbors[i]->enda.y==end2.y))
				endc[i] = 0; // enda
			else
				endc[i] = 1; // endb
		}
		endc[C_Num-1] = 1; // endb
	}


	old_array_of_neighbors[C_Num-1] = first_connected_end;

	//CGA: Create Array of New Objects
	lineObj **new_array_of_neighbors = (lineObj **)malloc(C_Num*sizeof(lineObj*));
	
	for(int i=0;i<C_Num;i++) 
		new_array_of_neighbors[i] = (lineObj *)malloc(sizeof(lineObj));
	
	for(int i=0;i<C_Num;i++)
	{ 
		error = 1;
		if (endc[i]==0) // if they were connected trought enda, then the new segment connects to endb cause connection in enda is moving
		{
			new_array_of_neighbors[i]->enda = old_array_of_neighbors[i]->endb;
			new_array_of_neighbors[i]->enda_L_Num = old_array_of_neighbors[i]->endb_L_Num;
			new_array_of_neighbors[i]->enda_C_Num = old_array_of_neighbors[i]->endb_C_Num;

		}
		else			// endb
		{
			new_array_of_neighbors[i]->enda = old_array_of_neighbors[i]->enda;
			new_array_of_neighbors[i]->enda_L_Num = old_array_of_neighbors[i]->enda_L_Num;
			new_array_of_neighbors[i]->enda_C_Num = old_array_of_neighbors[i]->enda_C_Num;
		}

		new_array_of_neighbors[i]->endb  = new_endc;

		//CGA Make sure it is withing lenghts
		new_array_of_neighbors[i]->len = (int)sqrt(double(DistSquare(new_array_of_neighbors[i]->enda,new_array_of_neighbors[i]->endb)));
		if (new_array_of_neighbors[i]->len<L_MIN || new_array_of_neighbors[i]->len>L_MAX)
			break;

		//CGA find new coordinates for new object
		new_array_of_neighbors[i]->x = (int)floor(0.5*(double(new_array_of_neighbors[i]->enda.x+new_array_of_neighbors[i]->endb.x))+0.5);
		new_array_of_neighbors[i]->y = (int)floor(0.5*(double(new_array_of_neighbors[i]->enda.y+new_array_of_neighbors[i]->endb.y))+0.5);
		new_array_of_neighbors[i]->width = old_array_of_neighbors[i]->width;
		new_array_of_neighbors[i]->type = old_array_of_neighbors[i]->type;

		//CGA find new theta
		if(new_array_of_neighbors[i]->enda.y < new_array_of_neighbors[i]->endb.y)
		{
			if (new_array_of_neighbors[i]->enda.x-new_array_of_neighbors[i]->endb.x != 0)
			    theta = atan(double((new_array_of_neighbors[i]->endb.y-new_array_of_neighbors[i]->enda.y)/double(new_array_of_neighbors[i]->enda.x-new_array_of_neighbors[i]->endb.x)));
			else
				theta = _PI/2;
			if (theta < 0)
				theta += _PI;
			new_array_of_neighbors[i]->theta = theta;
		}
		else
		{
			if (new_array_of_neighbors[i]->enda.x-new_array_of_neighbors[i]->endb.x != 0)
			    theta = atan(double((new_array_of_neighbors[i]->enda.y-new_array_of_neighbors[i]->endb.y)/double(new_array_of_neighbors[i]->endb.x-new_array_of_neighbors[i]->enda.x)));
			else
				theta = _PI/2;
			if (theta < 0)
				theta += _PI;
			new_array_of_neighbors[i]->theta = theta;
		}

		new_array_of_neighbors[i]->endb_L_Num = 0;	
		new_array_of_neighbors[i]->endb_C_Num = 1;
		

		n_io += Bad_IO_connection_move(new_array_of_neighbors[i], M, old_array_of_neighbors,C_Num,&g_Rio);
		if(g_Rio == INF)
		{
			error =1;
			break;
		}
		total_gRio += g_Rio;

		n_eo += Bad_EO_connection_move(NEIGHBOORHOOD,searchRatio,new_array_of_neighbors[i], M, old_array_of_neighbors,C_Num, &g_Rc);
		total_gRc += g_Rc;
		
		n_io_before += Bad_IO(old_array_of_neighbors[i], M,&g_Rio_before);
		total_gRio_before += g_Rio_before;
		if(g_Rio_before == INF)
		{
			error =1;
			break;
		}

		n_eo_before += Bad_EO_death(NEIGHBOORHOOD,searchRatio,old_array_of_neighbors[i],M,&g_Rc_before);
		total_gRc_before += g_Rc_before;

		if (new_array_of_neighbors[i]->endb_C_Num != 0)
			break;
		
		dtmp = dataterm_rec(img_seg, new_array_of_neighbors[i],patch,patch_len,new_array_of_neighbors[i]->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
		new_array_of_neighbors[i]->dataterm = dtmp;
		
		dterm += dtmp;
		dterm_before += old_array_of_neighbors[i]->dataterm;
		error = 0;
	}
	if(fabs(total_gRc_before) >= INF || fabs(total_gRc) >= INF || fabs(total_gRio) >= INF || fabs(total_gRio_before) >= INF)
		error = 1;

	if (!error)
	{
		n_io += Bad_IO_objects(new_array_of_neighbors, C_Num, &g_Rio);
		n_eo += Bad_EO_objects(NEIGHBOORHOOD,searchRatio,new_array_of_neighbors, C_Num, &g_Rc);
		
		Echange = -gamma_d*dterm -w_io*g_Rio -w_eo*g_Rc;
		exp_Echange = beta*exp(Echange);

		Echange_before = -gamma_d*dterm_before -w_io*g_Rio_before -w_eo*g_Rc_before;
		exp_Echange_before = beta*exp(Echange_before);

		if(fabs(Echange_before) >= INF/10 || fabs(Echange)>= INF/10)
		{
			free(old_array_of_neighbors);
			free(endc);
			return 0;
		}


		double R = pow(exp_Echange/exp_Echange_before,1.0/T);
		
		double r = random2();
		if (r < MIN(1,R) && (Echange-Echange_before))
		{
			for(int i=0;i<C_Num;i++) 
				killlineSeg(old_array_of_neighbors[i],M);
			
			M->Vo += (dterm - dterm_before);
			M->VRio += (g_Rio- g_Rio_before);
			M->VReo += (g_Rc - g_Rc_before);

			for(int i=0;i<C_Num;i++)
			{
				

				if(i==0) new_array_of_neighbors[i]->type--; // first item has no connection yet in that connection point.
				new_array_of_neighbors[i]->endb_L_Num = 0;
				new_array_of_neighbors[i]->enda_L_Num = 0;
				new_array_of_neighbors[i]->endb_C_Num = 0;
				new_array_of_neighbors[i]->enda_C_Num = 0;
				Bad_EO(NEIGHBOORHOOD,searchRatio,new_array_of_neighbors[i],M,&g_Rc);
				AddlineSeg(new_array_of_neighbors[i],M);
			}

			//if(Echange - Echange_before < 100 && (Echange - Echange_before) > 0)
			//{
				return_energy = Echange - Echange_before;
			//}
			//else
			//{
			//	return_energy = 0;
			//}



			#if INCLUDE_OPENCV_TRANSITION
				printf("Dterm: %f\n", dterm - dterm_before);
				printf("Rc: %.2f\n", g_Rc - g_Rc_before);
				printf("Rio: %.2f\n", g_Rio - g_Rio_before);
				
				printf("Single/Double Connection Move Energy: %.2f\n", Echange-Echange_before);
				display_only_one_double(img, img_height, img_width, new_array_of_neighbors[0], 1);
			#endif


		}
		else
		{
			for(int i=0;i<C_Num;i++) 
				free(new_array_of_neighbors[i]);
			return_energy = 0;
		}
	}
	else
	{
		for(int i=0;i<C_Num;i++) 
			free(new_array_of_neighbors[i]);
		return_energy = 0;
	}
	
	free(old_array_of_neighbors);
	free(endc);
	return return_energy;
}
