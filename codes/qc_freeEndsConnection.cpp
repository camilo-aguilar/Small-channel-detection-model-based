#include "QualityCandy.h"
#include "tiff.h"

NClinks *NClinksInit()
{
    NClinks *L;
    L = (NClinks *)malloc(sizeof(NClinks));  
    if(L == NULL)                     
        printf("fail to allocate NClink\n");
    L->next = NULL;   

	return L;
}


void AddNClinks(NClinks *l, int index, lineObj *obj1, lineObj *obj2)
{
	NClinks *pre = l;
	for (int i = 1; i <index;i++)
		pre=pre->next;
	NClinks *p;

	p = (NClinks *)malloc(sizeof(NClinks));

	p->object1 = obj1;
	p->object2 = obj2;
	p->next = pre->next;
	pre->next = p;
}

void KillNClinks(NClinks *l, lineObj *obj1, lineObj *obj2)
{
	NClinks *p,*pre;
	pre = l;
	p = l->next;
	
	int if_found = ((p->object1 == obj1)&&(p->object2 == obj2))||((p->object2 == obj1)&&(p->object1 == obj2));
	while(!if_found)
	{
		pre = p;
		p = p->next;
		if_found = ((p->object1 == obj1)&&(p->object2 == obj2))||((p->object2 == obj1)&&(p->object1 == obj2));
	}
	pre->next = p->next;

	free(p);
}

void NC_pairs(lineObj *obj1, lineObj *obj2,site *end1, site *end2)
{
	if(DIST(obj1->enda.x,obj1->enda.y,obj2->enda.x,obj2->enda.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
	{
		*end1 = obj1->enda;
		*end2 = obj2->enda;
	}
	if(DIST(obj1->endb.x,obj1->endb.y,obj2->enda.x,obj2->enda.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
	{
		*end1 = obj1->endb;
		*end2 = obj2->enda;
	}
	if(DIST(obj1->enda.x,obj1->enda.y,obj2->endb.x,obj2->endb.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
	{
		*end1 = obj1->enda;
		*end2 = obj2->endb;
	}
	if(DIST(obj1->endb.x,obj1->endb.y,obj2->endb.x,obj2->endb.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
	{
		*end1 = obj1->endb;
		*end2 = obj2->endb;
	}
}

Node* findObj(LinkedList l,lineObj *obj)
{
	Node *pre = l;
	Node *p = pre->next;
	while(p->index != obj)
	{
		pre = pre->next;
		p=pre->next;
	}
	return(pre);
}
void killlineSeg(lineObj *obj, Candy *M)
{
	int type = obj->type;
	if(type == 0)
	{
		Node *pre = findObj(M->link_f,obj);
		Node *p = pre->next;
		pre->next = p->next;
		M->n_f--;
		UpdateItsNeighboors_death_freeSeg(obj,M);
		free(obj);
	}
	else if(type == 1)
	{
		Node *pre = findObj(M->link_s,obj);
		Node *p = pre->next;
		pre->next = p->next;
		M->n_s--;
		UpdateItsNeighboors_death_SingleSeg(obj,M);
		free(obj);
	}
	else if(type == 2)
	{
		Node *pre = findObj(M->link_d,obj);
		Node *p = pre->next;
		pre->next = p->next;
		M->n_d--;
		UpdateItsNeighboors_death_DoubleSeg(obj,M);
		free(obj);
	}
	else
	{
//KDW		printf("forget to assign type to this obj\n");
	}
}

void AddlineSeg(lineObj *obj, Candy *M)
{
	int type = obj->type;
	if(type == 0)
	{
		M->n_f++;
		LinkedListInsert( M->link_f,M->n_f, obj);  //this is a free segment	
		UpdateItsNeighboors_born_freeSeg(obj,M);
	}
	else if(type == 1)
	{
		M->n_s++;
		LinkedListInsert( M->link_s,M->n_s, obj);  //this is a single segment
		UpdateItsNeighboors_born_SingleSeg(obj,M);
	}
	else if(type == 2)
	{
		M->n_d++;
		LinkedListInsert( M->link_d,M->n_d, obj);  //this is a double segment
		UpdateItsNeighboors_born_DoubleSeg(obj,M);
	}
	else
	{
//KDW		printf("forget to assign type to this obj\n");
	}
}

int Conncet_freeEnds(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, double **patch, int patch_len,int img_height, int img_width, double ***Matrix, double prior_pen1, double prior_pen2, double T)
{
	double w_f = M->w_f;
	double w_s = M->w_s;
	double w_d = M->w_d;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int C_n = M->connection_n;
	int N_n = M->neighbor_n;

	if(N_n == 0)
		return 0;

	int choose_num = (int)floor(random2()*(N_n-1)+1+0.5);
	NClinks *pre = M->neighbor_l;

	for (int i = 1;i<choose_num;i++)
		pre = pre->next;

	NClinks *p = pre->next;

	lineObj *obj1,*obj2;
	if(random2()<0.5)
	{
		obj1 = p->object1;
		obj2 = p->object2;
	}
	else
	{
		obj1 = p->object2;
		obj2 = p->object1;
	}

	site end1,end2;

	NC_pairs(obj1,obj2,&end1,&end2);  //find the neighoor set end1 & end2

//	site tmp;
//	tmp.x = 161; tmp.y = 128;
//	if(DistSquare(end1,tmp)<5)
//		printf("..");
	site new_enda, new_endb;  //move obj1's end poit to its neighbor

	if ((obj1->enda.x == end1.x)&&(obj1->enda.y==end1.y))
		new_enda = obj1->endb;
	else
		new_enda = obj1->enda;
	new_endb = end2;
//KDW	int new_length = sqrt(double(DistSquare(new_enda,new_endb)));
	int new_length = (int)sqrt(double(DistSquare(new_enda,new_endb))); //KDW

	if (new_length<L_MIN || new_length>L_MAX)
		return 0;

	lineObj *new_seg=(lineObj *)malloc(sizeof(lineObj));
	new_seg->enda = new_enda;
	new_seg->endb = new_endb;
	new_seg->len = new_length;
	new_seg->width = obj1->width;
	new_seg->x = (int)floor(0.5*(double(new_seg->enda.x+new_seg->endb.x))+0.5);
	new_seg->y = (int)floor(0.5*(double(new_seg->enda.y+new_seg->endb.y))+0.5);

	double theta;

	if(new_seg->enda.y < new_seg->endb.y)
	{
		if (new_seg->enda.x-new_seg->endb.x != 0)
		 theta = atan(double((new_seg->endb.y-new_seg->enda.y)/double(new_seg->enda.x-new_seg->endb.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
		new_seg->theta = theta;
	}
	else
	{
		if (new_seg->enda.x-new_seg->endb.x != 0)
		 theta = atan(double((new_seg->enda.y-new_seg->endb.y)/double(new_seg->endb.x-new_seg->enda.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
		new_seg->theta = theta;
	}

	new_seg->endb_L_Num = 0;
	new_seg->enda_L_Num = 0;
	new_seg->enda_C_Num = 0;
	new_seg->endb_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0;
	double g_Rio_before = 0;
	double g_Rc =0.;
	double g_Rc_before =0.;
	int n_io = Bad_IO_freeEndsC(new_seg, obj1, M,&g_Rio);
	int n_eo = Bad_EO_freeEndsC(NEIGHBOORHOOD,searchRatio,new_seg,M,obj1,&g_Rc);
	int n_io_before = Bad_IO(obj1, M,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,obj1,M,&g_Rc_before);

	//should add data term here, assume homogeneous Poisson 
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (new_seg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	int w_num = new_seg->width-W_MIN;
	//double dterm = dataterm_double(mid_img[img_num],new_seg->enda,new_seg->endb);
	double dterm = dataterm_rec(img_seg, new_seg,patch,patch_len,new_seg->width, Matrix, prior_pen1,prior_pen2, img_height,img_width);
	double w = 0;
	double Echange_before;
	if(new_seg->enda_C_Num != 0 && new_seg->endb_C_Num != 0)
	{
		new_seg->type = 2;   // this is a double seg
		w = M->w_d;
#ifdef QUALITY_CANDY
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_s-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_s-w_io*n_io_before-w_eo*n_eo_before);
#endif
	}
	else if(new_seg->enda_C_Num != 0 || new_seg->endb_C_Num != 0)
	{
		new_seg->type = 1;    //this is a single seg
		w = M->w_s;
#ifdef QUALITY_CANDY
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_f-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_f-w_io*n_io_before-w_eo*n_eo_before);
#endif
	}
	//KDW else
	//KDW	printf("free connection wrong\n");
	double Echange_from_neighboor = 0;
	if(obj2->type == 0)
		Echange_from_neighboor = M->w_s - M->w_f;
	if(obj2->type == 1)
		Echange_from_neighboor = M->w_d - M->w_s;

#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w-w_io*g_Rio-w_eo*g_Rc-Echange_from_neighboor);
	new_seg->engergy_for_transition = beta*exp(-gamma_d*dterm-w-w_io*g_Rio-w_eo*g_Rc);
#else
	double Echange = beta*exp(-gamma_d*dterm-w-w_io*n_io-w_eo*n_eo-Echange_from_neighboor);
	new_seg->engergy_for_transition = beta*exp(-gamma_d*dterm-w-w_io*n_io-w_eo*n_eo);
#endif
	new_seg->dataterm = dterm;
	new_seg->img_num = img_num;

//	double Echange_before = obj1->engergy_for_transition;

	double R = pow(Echange/Echange_before,T) *(M->p_c_CtoF)* N_n /((M->p_c_FtoC)*(C_n+1)*(M->lambda*NEIGHBOORHOOD*NEIGHBOORHOOD/(img_height*img_width))); //should be  lambda*NEIGHBOORHOOD*NEIGHBOORHOOD/(img_height*img_width)
	//printf("r= %f\n",R);
	double r= random2();
	if(r<MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
		//connect the two ends
		killlineSeg(obj1,M);
		AddlineSeg(new_seg,M);
	}
	else
	{
		free(new_seg);
	}

	return 0;
}


int Seperate_connectedEnds(Candy *M, double **img, int **img_seg, double ****img_mpp_l, double ****img_seg_l, double **patch, int patch_len,int img_height, int img_width, double ***Matrix, double prior_pen1, double prior_pen2 ,double T)
{
	double w_s = M->w_s;
	double w_d = M->w_d;
	double w_io = M->w_io;
	double w_eo = M->w_eo;
	double beta = M->beta;
	double gamma_d = M->gamma_d;
	int C_n = M->connection_n;
	int N_n = M->neighbor_n;

	if(C_n == 0)
		return 0;

	int choose_num = (int)floor(random2()*(C_n-1)+1+0.5);
	NClinks *pre = M->connection_l;

	for (int i = 1;i<choose_num;i++)
		pre = pre->next;

	NClinks *p = pre->next;

	lineObj *obj1,*obj2;
	if(random2()<0.5)
	{
		obj1 = p->object1;
		obj2 = p->object2;
	}
	else
	{
		obj1 = p->object2;
		obj2 = p->object1;
	}

	site end1,end2;

	NC_pairs(obj1,obj2,&end1,&end2);  //find the connection set end1 & end2

	site new_enda, new_endb;  //move obj1's end poit to its neighbor

	if ((obj1->enda.x == end1.x)&&(obj1->enda.y==end1.y))
		new_enda = obj1->endb;
	else
		new_enda = obj1->enda;

	int x_offset = (int)floor(random2()*2*NEIGHBOORHOOD-NEIGHBOORHOOD+0.5);
	int y_offset = (int)floor(random2()*2*NEIGHBOORHOOD-NEIGHBOORHOOD+0.5);
	int dist = x_offset*x_offset + y_offset*y_offset;

	while(dist== 0 || dist >= NEIGHBOORHOOD*NEIGHBOORHOOD )
	{
		x_offset = (int)floor(random2()*2*NEIGHBOORHOOD-NEIGHBOORHOOD+0.5);
		y_offset = (int)floor(random2()*2*NEIGHBOORHOOD-NEIGHBOORHOOD+0.5);
		dist = x_offset*x_offset + y_offset*y_offset;
	}

	new_endb.x = end2.x+x_offset;
	new_endb.y = end2.y + y_offset;

	if(new_endb.x <0 || new_endb.x >= img_width || new_endb.y <0 || new_endb.y >= img_height)
		return 0;

	int new_length = (int)sqrt(double(DistSquare(new_enda,new_endb)));

	if (new_length<L_MIN || new_length>L_MAX)
		return 0;

	lineObj *new_seg=(lineObj *)malloc(sizeof(lineObj));
	new_seg->enda = new_enda;
	new_seg->endb = new_endb;
	new_seg->len = new_length;
	new_seg->x = (int)floor(0.5*(double(new_seg->enda.x+new_seg->endb.x))+0.5);
	new_seg->y = (int)floor(0.5*(double(new_seg->enda.y+new_seg->endb.y))+0.5);

	double theta;

	if(new_seg->enda.y < new_seg->endb.y)
	{
		if (new_seg->enda.x-new_seg->endb.x != 0)
		 theta = atan(double((new_seg->endb.y-new_seg->enda.y)/double(new_seg->enda.x-new_seg->endb.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
		new_seg->theta = theta;
	}
	else
	{
		if (new_seg->enda.x-new_seg->endb.x != 0)
		 theta = atan(double((new_seg->enda.y-new_seg->endb.y)/double(new_seg->endb.x-new_seg->enda.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
		new_seg->theta = theta;
	}
	new_seg->width = obj1->width;
	new_seg->endb_L_Num = 0;
	new_seg->enda_L_Num = 0;
	new_seg->enda_C_Num = 0;
	new_seg->endb_C_Num = 0;

	double searchRatio = 0.25;
	double g_Rio = 0; 
	double g_Rio_before = 0; 
	double g_Rc =0.;
	double g_Rc_before =0.;
	int n_io = Bad_IO_freeEndsC(new_seg, obj1, M,&g_Rio);
	int n_eo = Bad_EO_freeEndsC(NEIGHBOORHOOD,searchRatio,new_seg,M,obj1,&g_Rc);
	int n_io_before = Bad_IO(obj1, M,&g_Rio_before);
	int n_eo_before = Bad_EO_death(NEIGHBOORHOOD,searchRatio,obj1,M,&g_Rc_before);

	//should add data term here, assume homogeneous Poisson 
	int img_num = 0;
	double iterval = _PI/(RADIUS_SEGS);
	while (new_seg->theta > iterval )
	{
		img_num++;
		iterval+=_PI/RADIUS_SEGS;
	}
	//if (random2()<0.5)
	//	img_num += int(RADIUS_SEGS/2.0);

	if (img_num >= int(RADIUS_SEGS))
		img_num--;
	int w_num = new_seg->width-W_MIN;
	//double dterm = dataterm_double(mid_img[img_num],new_seg->enda,new_seg->endb);
	double dterm = dataterm_rec(img_seg,new_seg,patch,patch_len,new_seg->width,Matrix, prior_pen1, prior_pen2, img_height,img_width);
	double w = 0;
	double Echange_before;
	if(new_seg->enda_C_Num != 0 && new_seg->endb_C_Num != 0)
	{
		//free(new_seg);
		//return 0; // this can not be happen
		new_seg->type = 2;   // this is a double seg
		w = M->w_d;
#ifdef QUALITY_CANDY
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_d-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_d-w_io*n_io_before-w_eo*n_eo_before);
#endif
	}
	else if(new_seg->enda_C_Num != 0 || new_seg->endb_C_Num != 0)
	{
		new_seg->type = 1;    //this is a single seg
		w = M->w_s;
#ifdef QUALITY_CANDY
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_d-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_d-w_io*n_io_before-w_eo*n_eo_before);
#endif
	}
	else
	{
		new_seg->type = 0;
		w = M->w_f;
#ifdef QUALITY_CANDY
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_s-w_io*g_Rio_before-w_eo*g_Rc_before);
#else
		Echange_before = beta*exp(-gamma_d*obj1->dataterm-w_s-w_io*n_io_before-w_eo*n_eo_before);
#endif
	}

	double Echange_from_neighboor = 0;
	if(obj2->type == 2)
		Echange_from_neighboor = M->w_s - M->w_d;
	if(obj2->type == 1)
		Echange_from_neighboor = M->w_f - M->w_s;

#ifdef QUALITY_CANDY
	double Echange = beta*exp(-gamma_d*dterm-w-w_io*g_Rio-w_eo*g_Rc-Echange_from_neighboor);
	new_seg->engergy_for_transition = beta*exp(-gamma_d*dterm-w-w_io*g_Rio-w_eo*g_Rc);
#else
	double Echange = beta*exp(-gamma_d*dterm-w-w_io*n_io-w_eo*n_eo-Echange_from_neighboor);
	new_seg->engergy_for_transition = beta*exp(-gamma_d*dterm-w-w_io*n_io-w_eo*n_eo);
#endif
	new_seg->dataterm = dterm;
	new_seg->img_num = img_num;

//	double Echange_before = obj1->engergy_for_transition;

	double R = pow(Echange/Echange_before,T) *((M->p_c_FtoC)*C_n*(M->lambda*NEIGHBOORHOOD*NEIGHBOORHOOD/(img_height*img_width)))/(M->p_c_CtoF)* (N_n+1); //should be  lambda*NEIGHBOORHOOD*NEIGHBOORHOOD/(img_height*img_width)
	
	//if (R< 0.0001)
	//	R = 0.0001;

	//R = 1/R;
	if (R>0.1)
	{
	//	printf("r = %f, echange = %f, echange_before = %f\n",R,Echange,Echange_before);
	}
	double r= random2();
	if(r<MIN(1,R))
	{
	if(n_io>10000)
		printf("..");
		//connect the two ends
		killlineSeg(obj1,M);
		AddlineSeg(new_seg,M);
	}
	else
	{
		free(new_seg);
	}

	return 0;
}