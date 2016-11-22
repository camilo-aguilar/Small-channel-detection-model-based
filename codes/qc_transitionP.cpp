#include "QualityCandy.h"
#include "tiff.h"

#define DELTA_LEN 2 //KDW 3
#define DELTA_THETA _PI/20
#define DELTA_WIDTH 2 //KDW

#define END_MOVE_RANGE 2

lineObj* select_object_from_link(LinkedList link,int num, int choose_num)
{
	Node *p = link;
	for (int i = 0;i< choose_num;i++)
		p = p->next;
	return (p->index);
}

int FreeSeg_length_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
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

	//lineObj* object = select_object_from_link(M->link_f,n_f,choose_num);

	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
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

	site old_enda = object->enda;
	site old_endb = object->endb;
	int old_len = object->len;
//	double Echange_before = object->engergy_for_transition;

	double searchRatio = 0.25;

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

	double g_Rio = 0; 
	double g_Rio_before = 0;
	double g_Rc =0.;
	double g_Rc_before =0.;
	int w_num = freeSeg->width-W_MIN;
	int n_io = Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	int n_eo = Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	int n_io_before = Bad_IO(object, M ,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);


	if (freeSeg->enda_C_Num != 0 || freeSeg->endb_C_Num != 0)
	{
		free(freeSeg);
		return 0;
	}

	int img_num = object->img_num;
	//double dterm = dataterm(mid_img[img_num],new_enda,new_endb);
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	//double dterm = dataterm(img,freeSeg->enda,freeSeg->endb,freeSeg->theta,freeSeg->len,img_num,patch,patch_len);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*n_io-w_eo*n_eo);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*n_io_before-w_eo*n_eo_before);
#endif
	double R = pow(Echange/Echange_before,T);
	double r = random2();
	if (r < MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
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
	}
	else
	{
		free(freeSeg);
	}
	return 0;
}

int FreeSeg_width_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
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

	//lineObj* object = select_object_from_link(M->link_f,n_f,choose_num);

	lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));
	int d_width = (int)floor(random2()*2*DELTA_WIDTH - DELTA_WIDTH +0.5);
	double new_width = d_width + object->width;

	if(new_width < W_MIN || new_width > W_MAX)
	{
			free(freeSeg);
			return 0;
	}

	double theta = object->theta;
	int c_x = object->x;
	int c_y = object->y;

	site new_enda = object->enda;
	site new_endb = object->endb;

	site old_enda = object->enda;
	site old_endb = object->endb;
	int old_len = object->len;


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
//	double Echange_before = object->engergy_for_transition;
	double g_Rio = 0; 
	double g_Rc =0.;
	int w_num = freeSeg->width-W_MIN;
	int n_io = Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio); // we don't need this operation
	int n_eo = Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);// we don't need this operation

	if (freeSeg->enda_C_Num != 0 || freeSeg->endb_C_Num != 0)
	{
		free(freeSeg);
		return 0;
	}

	int img_num = object->img_num;
	//double dterm = dataterm(mid_img[img_num],new_enda,new_endb);
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	//double dterm = dataterm(img,freeSeg->enda,freeSeg->endb,freeSeg->theta,freeSeg->len,img_num,patch,patch_len);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*g_Rio-w_eo*g_Rc);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*n_io-w_eo*n_eo);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*n_io-w_eo*n_eo);
#endif
	double R = pow(Echange/Echange_before,T);
	double r = random2();
	if (r < MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
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
	}
	else
	{
		free(freeSeg);
	}
	return 0;
}

int FreeSeg_theta_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;

	if (n_f == 0)
		return 0;
	
//KDW	int choose_num = floor(random2()*(n_f-1)+1+0.5);
	int choose_num = (int)floor(random2()*(n_f-1)+1+0.5); //KDW
	
	Node *pre = M->link_f;
	for (int i = 1;i<choose_num;i++)
		pre = pre->next;
	Node *p = pre->next;

	lineObj *object = p->index;

	//lineObj* object = select_object_from_link(M->link_f,n_f,choose_num);

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
	
	site old_enda = object->enda;
	site old_endb = object->endb;
	double old_theta = object->theta;

	double searchRatio = 0.25;

	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (new_theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

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
 
//	double Echange_before = object->engergy_for_transition;
	if((freeSeg->x == 62)&&(freeSeg->y==138))
		printf("..");
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;
	int w_num = freeSeg->width-W_MIN;
	int n_io = Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	int n_eo = Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	int n_io_before = Bad_IO(object, M,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

	if (freeSeg->enda_C_Num != 0 || freeSeg->endb_C_Num != 0)
	{
		free(freeSeg);
		return 0;
	}
	  
	//double dterm = dataterm(mid_img[img_num],new_enda,new_endb);
//	double dterm = dataterm_rec(freeSeg,img_mpp_l[w_num][img_num],img_seg_l[w_num][img_num],patch,patch_len,freeSeg->width);
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width, Matrix,prior_pen1, prior_pen2, img_height,img_width);
	
	//double dterm = dataterm(img,freeSeg->enda,freeSeg->endb,freeSeg->theta,freeSeg->len,img_num,patch,patch_len);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*n_io-w_eo*n_eo);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*n_io_before-w_eo*n_eo_before);
#endif
	double R = pow(Echange/Echange_before,T);

	double r = random2();
	if (r < MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
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
	}
	else
	{
		free(freeSeg);
	}
	return 0;
}

int FreeSeg_freeEnd_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
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
	//lineObj* object = select_object_from_link(M->link_f,n_f,choose_num);

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
//	double Echange_before = object->engergy_for_transition;

//	if((freeSeg->x == 155)&&(freeSeg->y==126))
//		printf("..");
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;
	int n_io = Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	int n_eo = Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	int n_io_before = Bad_IO(object, M,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

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
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	int w_num = freeSeg->width-W_MIN;
//	double dterm = dataterm(mid_img[img_num],freeSeg->enda,freeSeg->endb);
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*n_io-w_eo*n_eo);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*n_io_before-w_eo*n_eo_before);
#endif
	freeSeg->dataterm = dterm;
	freeSeg->img_num = img_num;
	freeSeg->engergy_for_transition = Echange;
	freeSeg->type = 0;

	double R = pow(Echange/Echange_before,T);
	
	double r = random2();
	if (r < MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
		pre->next = p->next;  // kill the object, add the freeSeg
		M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);
		M->n_f++;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);
	}
	else
	{
			free(freeSeg);
	}
	return 0;

}

int FreeSeg_center_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
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
	//lineObj* object = select_object_from_link(M->link_f,n_f,choose_num);

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
		return 0;
	new_endb.x = object->endb.x+delta_x;
	new_endb.y = object->endb.y+delta_y;
	if (new_endb.x < 0 ||new_endb.x >= img_width ||new_endb.y < 0 ||new_endb.y >= img_height)
		return 0;
	

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
//	double Echange_before = object->engergy_for_transition;
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;

	int n_io = Bad_IO_transition(freeSeg, M,choose_num-1,0,&g_Rio);
	int n_eo = Bad_EO_transition(NEIGHBOORHOOD,searchRatio,freeSeg,M,choose_num-1,0,&g_Rc);
	int n_io_before = Bad_IO(object, M,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

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
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	int w_num = freeSeg->width-W_MIN;
//	double dterm = dataterm(mid_img[img_num],freeSeg->enda,freeSeg->endb);
	double dterm = dataterm_rec(img_seg, freeSeg,patch,patch_len,freeSeg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*g_Rio-w_eo*g_Rc);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_f-w_io*n_io-w_eo*n_eo);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_f-w_io*n_io_before-w_eo*n_eo_before);
#endif
	freeSeg->dataterm = dterm;
	freeSeg->img_num = img_num;
	freeSeg->engergy_for_transition = Echange;
	freeSeg->type = 0;

	double R = pow(Echange/Echange_before,T);
	
	double r = random2();
	if (r < MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
		pre->next = p->next;  // kill the object, add the freeSeg
		M->n_f--;
		UpdateItsNeighboors_death_freeSeg(object,M);
		free(object);
		M->n_f++;
		LinkedListInsert( M->link_f,M->n_f, freeSeg);  	
		UpdateItsNeighboors_born_freeSeg(freeSeg,M);
	}
	else
	{
			free(freeSeg);
	}
	return 0;

}

int SingleSeg_freeEnd_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
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
	//	printf("single seg number wrong\n");
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

	int n_io = Bad_IO_transition(SingleSeg, M,choose_num-1,1,&g_Rio);
	int n_eo = Bad_EO_transition(NEIGHBOORHOOD,searchRatio,SingleSeg,M,choose_num-1,1,&g_Rc);
	int n_io_before = Bad_IO(object, M,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,object,M,&g_Rc_before);

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
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	int w_num = SingleSeg->width-W_MIN;
//	double dterm = dataterm(mid_img[img_num],SingleSeg->enda,SingleSeg->endb);
	double dterm = dataterm_rec(img_seg, SingleSeg,patch,patch_len,SingleSeg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w_s-w_io*g_Rio-w_eo*g_Rc);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_s-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
	double Echange = beta*exp(-gamma_d*dterm-w_s-w_io*n_io-w_eo*n_eo);
	double Echange_before = beta*exp(-gamma_d*object->dataterm-w_s-w_io*n_io_before-w_eo*n_eo_before);
#endif
	SingleSeg->dataterm = dterm;
	SingleSeg->img_num = img_num;
	SingleSeg->engergy_for_transition = Echange;
	SingleSeg->type = 1;

	double R = pow(Echange/Echange_before,T);
	
	double r = random2();
	if (r < MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
		pre->next = p->next;  // kill the object, add the freeSeg
		M->n_s--;
		UpdateItsNeighboors_death_SingleSeg(object,M);
		free(object);
		M->n_s++;
		LinkedListInsert( M->link_s,M->n_s, SingleSeg);  	
		UpdateItsNeighboors_born_SingleSeg(SingleSeg,M);
	}
	else
		{
			free(SingleSeg);
	}
	return 0;

}

#if 1
int SingleDoubleSeg_Connection_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_s = M->w_s;
	double w_d = M->w_d;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;
	int n_s = M->n_s;
	int n_d = M->n_d;
	int C_n = M->connection_n;
	int error = 1;
	double g_Rio = 0;
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;

	if(C_n == 0)
		return 0;

	int choose_num = (int)floor(random2()*(C_n-1)+1+0.5);
	NClinks *pre = M->connection_l;

	for (int i = 1;i<choose_num;i++)
		pre = pre->next;

	NClinks *p = pre->next;
	lineObj *obj1 = p->object1;
	lineObj *obj2 = p->object2;

	site end1,end2;

	NC_pairs(obj1,obj2,&end1,&end2);  //find the connection set end1 & end2

	site new_endc;
	int x_offset;
	int y_offset;
	int n=0;
	do{
		x_offset = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
		y_offset = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
		n++;
		if(n>1000) return 0;
	}while((x_offset==0)&&(y_offset==0));
	new_endc.x = end1.x + x_offset;
	new_endc.y = end1.y + y_offset;
	if(new_endc.x <0 || new_endc.x >= img_width || new_endc.y <0 || new_endc.y >= img_height)
		return 0;

	lineObj **obj;
	int C_Num, *endc; // endc : connection point 0:enda, 1:endb
	if ((obj1->enda.x == end1.x)&&(obj1->enda.y==end1.y)){
//		if(obj1->enda_C_Num>1) return 0;
		C_Num = obj1->enda_C_Num+1;
		obj = (lineObj **)malloc(C_Num*sizeof(lineObj *)); // obj1 and its connections(C_Num)
		endc = (int *)malloc(C_Num*sizeof(int));
		for(int i=0;i<C_Num-1;i++){ 
			obj[i] = obj1->enda_C[i];
			NC_pairs(obj1,obj[i],&end1,&end2);
			if ((obj[i]->enda.x == end2.x)&&(obj[i]->enda.y==end2.y))
				endc[i] = 0; // enda
			else
				endc[i] = 1; // endb
		}
		endc[C_Num-1] = 0; // enda
	}
	else{
//		if(obj1->endb_C_Num>1) return 0;
		C_Num = obj1->endb_C_Num+1;
		obj = (lineObj **)malloc(C_Num*sizeof(lineObj *)); // obj1 and its connections(C_Num)
		endc = (int *)malloc(C_Num*sizeof(int));
		for(int i=0;i<C_Num-1;i++){ 
			obj[i] = obj1->endb_C[i];
			NC_pairs(obj1,obj[i],&end1,&end2);
			if ((obj[i]->enda.x == end2.x)&&(obj[i]->enda.y==end2.y))
				endc[i] = 0; // enda
			else
				endc[i] = 1; // endb
		}
		endc[C_Num-1] = 1; // endb
	}
	obj[C_Num-1] = obj1;

	double searchRatio = 0.25;
	double theta;
	int n_io = 0, n_eo = 0, n_io_before = 0, n_eo_before = 0;
	double dtmp, dterm = 0., dterm_before = 0.;
	double Echange = 0, Echange_before = 0.;
	lineObj **new_obj = (lineObj **)malloc(sizeof(lineObj*));

	for(int i=0;i<C_Num;i++) 
		new_obj[i] = (lineObj *)malloc(sizeof(lineObj));
	
	for(int i=0;i<C_Num;i++){ 
		error = 1;
		if (endc[i]==0) // enda
			new_obj[i]->enda = obj[i]->endb;
		else			// endb
			new_obj[i]->enda = obj[i]->enda;
		new_obj[i]->endb  = new_endc;
		new_obj[i]->len = (int)sqrt(double(DistSquare(new_obj[i]->enda,new_obj[i]->endb)));
		if (new_obj[i]->len<L_MIN || new_obj[i]->len>L_MAX)
			break;
		new_obj[i]->x = (int)floor(0.5*(double(new_obj[i]->enda.x+new_obj[i]->endb.x))+0.5);
		new_obj[i]->y = (int)floor(0.5*(double(new_obj[i]->enda.y+new_obj[i]->endb.y))+0.5);
		new_obj[i]->width = obj[i]->width;
		new_obj[i]->type = obj[i]->type;

		if(new_obj[i]->enda.y < new_obj[i]->endb.y)
		{
			if (new_obj[i]->enda.x-new_obj[i]->endb.x != 0)
			    theta = atan(double((new_obj[i]->endb.y-new_obj[i]->enda.y)/double(new_obj[i]->enda.x-new_obj[i]->endb.x)));
			else
				theta = _PI/2;
			if (theta < 0)
				theta += _PI;
			new_obj[i]->theta = theta;
		}
		else
		{
			if (new_obj[i]->enda.x-new_obj[i]->endb.x != 0)
			    theta = atan(double((new_obj[i]->enda.y-new_obj[i]->endb.y)/double(new_obj[i]->endb.x-new_obj[i]->enda.x)));
			else
				theta = _PI/2;
			if (theta < 0)
				theta += _PI;
			new_obj[i]->theta = theta;
		}

		new_obj[i]->endb_L_Num = 0;
		new_obj[i]->enda_L_Num = 0;
		new_obj[i]->endb_C_Num = 0;
		new_obj[i]->enda_C_Num = 0;

		n_io += Bad_IO_connection_move(new_obj[i], M, obj,C_Num,&g_Rio);
		n_eo += Bad_EO_connection_move(NEIGHBOORHOOD,searchRatio,new_obj[i], M, obj,C_Num, &g_Rc);
		n_io_before += Bad_IO(obj[i], M,&g_Rio_before);
		n_eo_before += Bad_EO_death(NEIGHBOORHOOD,searchRatio,obj[i],M,&g_Rc_before);
		if (new_obj[i]->endb_C_Num != 0)
			break;
		dtmp = dataterm_rec(img_seg, new_obj[i],patch,patch_len,new_obj[i]->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
		new_obj[i]->dataterm = dtmp;
		dterm += dtmp;
		dterm_before += obj[i]->dataterm;
		error = 0;
	}
//	if (0){
	if (!error){
		n_io += Bad_IO_objects(new_obj, C_Num, &g_Rio);
		n_eo += Bad_EO_objects(NEIGHBOORHOOD,searchRatio,new_obj, C_Num, &g_Rc);
#ifdef QUALITY_CANDY
		Echange = beta*exp(-gamma_d*dterm-w_io*g_Rio-w_eo*g_Rc);
		Echange_before = beta*exp(-gamma_d*dterm_before-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
		Echange = beta*exp(-gamma_d*dterm-w_io*n_io-w_eo*n_eo);
		Echange_before = beta*exp(-gamma_d*dterm_before-w_io*n_io_before-w_eo*n_eo_before);
#endif
//		for(int i=1;i<C_Num;i++){ 
//			new_obj[i]->engergy_for_transition = beta*exp(-gamma_d*new_obj[i]->dterm-w_io*n_io-w_eo*n_eo);

		double R = pow(Echange/Echange_before,T);
		
		double r = random2();
		if (r < MIN(1,R))
		{
	if(n_io>10000)
		printf("..");
			for(int i=0;i<C_Num;i++) 
				killlineSeg(obj[i],M);
			for(int i=0;i<C_Num;i++){
				if(i==0) new_obj[i]->type--; // first item has no connection yet in that connection point.
				new_obj[i]->endb_L_Num = 0;
				new_obj[i]->enda_L_Num = 0;
				new_obj[i]->endb_C_Num = 0;
				new_obj[i]->enda_C_Num = 0;
				Bad_EO(NEIGHBOORHOOD,searchRatio,new_obj[i],M,&g_Rc);
				AddlineSeg(new_obj[i],M);
			}
		}
		else{
			for(int i=0;i<C_Num;i++) 
				free(new_obj[i]);
			// I don't know why this doesn't work  // free(new_obj);
		}
	}
	else{
		for(int i=0;i<C_Num;i++) 
			free(new_obj[i]);
		//free(new_obj);
	}
	//free(new_obj);
	free(obj);
	free(endc);

	return !error;
}
#else
int SingleDoubleSeg_Connection_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2)
{
	double w_f = M->w_f;
	double w_s = M->w_s;
	double w_d = M->w_d;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int n_f = M->n_f;
	int n_s = M->n_s;
	int n_d = M->n_d;
	int C_n = M->connection_n;
	int error = 1;

	if(C_n == 0)
		return 0;

	int choose_num = (int)floor(random2()*(C_n-1)+1+0.5);
	NClinks *pre = M->connection_l;

	for (int i = 1;i<choose_num;i++)
		pre = pre->next;

	NClinks *p = pre->next;
	lineObj *obj1 = p->object1;
	lineObj *obj2 = p->object2;

	site end1,end2;

	NC_pairs(obj1,obj2,&end1,&end2);  //find the connection set end1 & end2

	site new_endc;
	int x_offset;
	int y_offset;
	int n=0;
	do{
		x_offset = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
		y_offset = (int)floor(random2()*2*END_MOVE_RANGE-END_MOVE_RANGE+0.5);
		n++;
		if(n>1000) return 0;
	}while((x_offset==0)&&(y_offset==0));
	new_endc.x = end1.x + x_offset;
	new_endc.y = end1.y + y_offset;
	if(new_endc.x <0 || new_endc.x >= img_width || new_endc.y <0 || new_endc.y >= img_height)
		return 0;

	lineObj **obj;
	int C_Num, *endc; // endc : connection point 0:enda, 1:endb
	if ((obj1->enda.x == end1.x)&&(obj1->enda.y==end1.y)){
		C_Num = obj1->enda_C_Num+1;
		obj = (lineObj **)malloc(C_Num*sizeof(lineObj *)); // obj1 and its connections(C_Num)
		endc = (int *)malloc(C_Num*sizeof(int));
		for(int i=0;i<C_Num-1;i++){ 
			obj[i] = obj1->enda_C[i];
			NC_pairs(obj1,obj[i],&end1,&end2);
			if ((obj[i]->enda.x == end2.x)&&(obj[i]->enda.y==end2.y))
				endc[i] = 0; // enda
			else
				endc[i] = 1; // endb
		}
		endc[C_Num-1] = 0; // enda
	}
	else{
		C_Num = obj1->endb_C_Num+1;
		obj = (lineObj **)malloc(C_Num*sizeof(lineObj *)); // obj1 and its connections(C_Num)
		endc = (int *)malloc(C_Num*sizeof(int));
		for(int i=0;i<C_Num-1;i++){ 
			obj[i] = obj1->endb_C[i];
			NC_pairs(obj1,obj[i],&end1,&end2);
			if ((obj[i]->enda.x == end2.x)&&(obj[i]->enda.y==end2.y))
				endc[i] = 0; // enda
			else
				endc[i] = 1; // endb
		}
		endc[C_Num-1] = 1; // endb
	}
	obj[C_Num-1] = obj1;

	double searchRatio = 0.25;
	double theta;
	int n_io = 0, n_eo = 0, n_io_before = 0, n_eo_before = 0;
	double dtmp, dterm = 0., dterm_before = 0.;
	double Echange = 0, Echange_before = 0.;
	lineObj *new_obj = (lineObj *)malloc(C_Num*sizeof(lineObj));

	for(int i=0;i<C_Num;i++){ 
		error = 1;
		if (endc[i]==0) // enda
			new_obj[i].enda = obj[i]->endb;
		else			// endb
			new_obj[i].enda = obj[i]->enda;
		new_obj[i].endb  = new_endc;
		new_obj[i].len = (int)sqrt(double(DistSquare(new_obj[i].enda,new_obj[i].endb)));
		if (new_obj[i].len<L_MIN || new_obj[i].len>L_MAX)
			break;
		new_obj[i].x = (int)floor(0.5*(double(new_obj[i].enda.x+new_obj[i].endb.x))+0.5);
		new_obj[i].y = (int)floor(0.5*(double(new_obj[i].enda.y+new_obj[i].endb.y))+0.5);
		new_obj[i].width = obj[i]->width;
		new_obj[i].type = obj[i]->type;

		if(new_obj[i].enda.y < new_obj[i].endb.y)
		{
			if (new_obj[i].enda.x-new_obj[i].endb.x != 0)
			    theta = atan(double((new_obj[i].endb.y-new_obj[i].enda.y)/double(new_obj[i].enda.x-new_obj[i].endb.x)));
			else
				theta = _PI/2;
			if (theta < 0)
				theta += _PI;
			new_obj[i].theta = theta;
		}
		else
		{
			if (new_obj[i].enda.x-new_obj[i].endb.x != 0)
			    theta = atan(double((new_obj[i].enda.y-new_obj[i].endb.y)/double(new_obj[i].endb.x-new_obj[i].enda.x)));
			else
				theta = _PI/2;
			if (theta < 0)
				theta += _PI;
			new_obj[i].theta = theta;
		}

		new_obj[i].endb_L_Num = 0;
		new_obj[i].enda_L_Num = 0;
		new_obj[i].endb_C_Num = 0;
		new_obj[i].enda_C_Num = 0;

		n_io += Bad_IO_connection_move(&new_obj[i], M, obj,C_Num);
		n_eo += Bad_EO_connection_move(NEIGHBOORHOOD,searchRatio,&new_obj[i], M, obj,C_Num);
		n_io_before += Bad_IO(obj[i], M);
		n_eo_before += Bad_EO_death(NEIGHBOORHOOD,searchRatio,obj[i],M);
		if (new_obj[i].endb_C_Num != 0)
			break;
		dtmp = dataterm_rec(img_seg, &new_obj[i],patch,patch_len,new_obj[i].width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
		new_obj[i].dataterm = dtmp;
		dterm += dtmp;
		dterm_before += obj[i]->dataterm;
		error = 0;
	}
//	if (0){
	if (!error){
		n_io += Bad_IO_objects(new_obj, C_Num);
		n_eo += Bad_EO_objects(NEIGHBOORHOOD,searchRatio,new_obj, C_Num);
		Echange = beta*exp(-gamma_d*dterm-w_io*n_io-w_eo*n_eo);
		Echange_before = beta*exp(-gamma_d*dterm_before-w_io*n_io_before-w_eo*n_eo_before);
//		for(int i=1;i<C_Num;i++){ 
//			new_obj[i].engergy_for_transition = beta*exp(-gamma_d*new_obj[i].dterm-w_io*n_io-w_eo*n_eo);

		double R = pow(Echange/Echange_before,T);
		
		double r = random2();
		if (r < MIN(1,R))
		{
			for(int i=0;i<C_Num;i++) 
				killlineSeg(obj[i],M);
			for(int i=0;i<C_Num;i++){
				if(i==0) new_obj[i].type--; // first item has no connection yet in that connection point.
				new_obj[i].endb_L_Num = 0;
				new_obj[i].enda_L_Num = 0;
				new_obj[i].endb_C_Num = 0;
				new_obj[i].enda_C_Num = 0;
				Bad_EO(NEIGHBOORHOOD,searchRatio,&new_obj[i],M);
				AddlineSeg(&new_obj[i],M);
			}
		}
		else
			free(new_obj); 
	}
	else
		free(new_obj); 
	free(obj);
	free(endc);

	return !error;
}


#endif