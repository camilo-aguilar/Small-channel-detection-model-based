#include "QualityCandy.h"

#define IO_RANGE_MARGIN 1.4


void check_io(site c_enda, site c_endb, site enda, site endb, int c_x, int c_y, int x, int y, 
			  double tau_ij, double c_theta, double theta, int *n_io, double *g_Rio)
{
	#define PI_OVER2_MINUS_DELTA_MIN_SQURE  1.3879131192656386475625 
	
	if (ABS(tau_ij - _PI/2) > DELTA_IO_MIN) // > Pi/8
		(*g_Rio)=INF;
	else
		(*g_Rio)=1-sigma(tau_ij*tau_ij,PI_OVER2_MINUS_DELTA_MIN_SQURE);

}



int Bad_IO(lineObj *Seg, Candy *M, double *g_Rio)
{
	int n_io = 0;
	(*g_Rio) = 0;
	double grio;

	Bad_IO_from_Link(Seg,M->link_f,M->n_f,&grio);
	(*g_Rio) += grio;
	if((*g_Rio)<INF)
	{
		Bad_IO_from_Link(Seg,M->link_s,M->n_s,&grio);
		(*g_Rio) += grio;
	}
	if((*g_Rio)<INF)
	{
		Bad_IO_from_Link(Seg,M->link_d,M->n_d,&grio);
		(*g_Rio) += grio;
	}
	
	return n_io;
}

/* 
1 - For each discovered line object
2     if the distance between each object's center is less than half len max(L1,L2) 
3     Check that their angles are at least more than Pi/8 difference.   
4     Assigns a potential difference  
 */
int Bad_IO_from_Link(lineObj *Seg, LinkedList link, int n, double *g_Rio)
{
	int n_io = 0;
	Node *p;
	p = link;

	int c_x = Seg->x;
	int c_y = Seg->y;
	int c_len = Seg->len;
	double c_theta = Seg->theta;
	site c_enda = Seg->enda;
	site c_endb = Seg->endb;
	double grio;
	(*g_Rio) = 0;

	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		int x= l->x;
		int y = l->y;
		int len = l->len;
		double theta = l->theta;
		site enda = l->enda;
		site endb = l->endb;

		int cret = 0;
		//if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= MAX(10,MAX(c_len,len)/2)|| sqrt(double((c_enda.x-x)*(c_enda.x-x)+(c_enda.y-y)*(c_enda.y-y))) <= (MAX(c_len,len)/2-2)||sqrt(double((c_endb.x-x)*(c_endb.x-x)+(c_endb.y-y)*(c_endb.y-y))) <= (MAX(c_len,len)/2-2))
		if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= (double)MAX(c_len,len)/2.*IO_RANGE_MARGIN)
		cret = 1;
		
		 if (l!= Seg &&  cret == 1)
		{
			if(l->enda.x == Seg->enda.x && l->enda.y == Seg->enda.y && l->endb.x == Seg->endb.x && l->endb.y == Seg->endb.y)
				n_io = 4; //KDW identical n_io = 4 INF 
			if(l->enda.x == Seg->endb.x && l->enda.y == Seg->endb.y && l->endb.x == Seg->enda.x && l->endb.y == Seg->enda.y)
				n_io = 4; //KDW identical n_io = 4 			
			double tau_ij = MIN(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
			check_io(c_enda, c_endb, enda, endb, c_x, c_y, x, y, tau_ij, c_theta, theta, &n_io, &grio);
			if((*g_Rio)<INF) (*g_Rio) += grio;
		}
	}
		return n_io;
}

int Bad_IO_from_Link_transition(lineObj *Seg, LinkedList link, int n, int choose_num , double *g_Rio)
{
	int n_io = 0;
	Node *p;
	p = link;

	int c_x = Seg->x;
	int c_y = Seg->y;
	int c_len = Seg->len;
	double c_theta = Seg->theta;
	site c_enda = Seg->enda;
	site c_endb = Seg->endb;
	double grio;
	(*g_Rio) = 0;

	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		if (i != choose_num)
		{
			int x= l->x;
			int y = l->y;
			int len = l->len;
			double theta = l->theta;
			site enda = l->enda;
			site endb = l->endb;

			int cret = 0;
			//if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= MAX(10,MAX(c_len,len)/2)|| sqrt(double((c_enda.x-x)*(c_enda.x-x)+(c_enda.y-y)*(c_enda.y-y))) <= (MAX(c_len,len)/2-2)||sqrt(double((c_endb.x-x)*(c_endb.x-x)+(c_endb.y-y)*(c_endb.y-y))) <= (MAX(c_len,len)/2-2))
			if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= (double)MAX(c_len,len)/2.*IO_RANGE_MARGIN)
				cret = 1;
			
			if (l!= Seg && cret == 1)
			{
				if(l->enda.x == Seg->enda.x && l->enda.y == Seg->enda.y && l->endb.x == Seg->endb.x && l->endb.y == Seg->endb.y)
					n_io =4;
				if(l->enda.x == Seg->endb.x && l->enda.y == Seg->endb.y && l->endb.x == Seg->enda.x && l->endb.y == Seg->enda.y)
					n_io =4;	

				double tau_ij = MAX(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
				check_io(c_enda, c_endb, enda, endb, c_x, c_y, x, y, tau_ij, c_theta, theta, &n_io,&grio);
				if((*g_Rio)<INF) (*g_Rio) += grio;
			}
		}
	}
	return n_io;
}

int Bad_IO_transition(lineObj *Seg, Candy *M, int choose_num, int end_type, double *g_Rio)
{
	int n_io = 0;
	(*g_Rio) = 0;
	double grio;
#ifdef QUALITY_CANDY
	if (end_type == 0)
	{
		Bad_IO_from_Link_transition(Seg,M->link_f,M->n_f,choose_num,&grio);
		(*g_Rio) += grio;
		if((*g_Rio)<INF){
			Bad_IO_from_Link(Seg,M->link_s,M->n_s,&grio);
			(*g_Rio) += grio;
		}
		if((*g_Rio)<INF){
			Bad_IO_from_Link(Seg,M->link_d,M->n_d,&grio);
			(*g_Rio) += grio;
		}
	}
	else if (end_type == 1)
	{
		Bad_IO_from_Link_transition(Seg,M->link_f,M->n_f,choose_num,&grio);
		(*g_Rio) += grio;
		if((*g_Rio)<INF){
			Bad_IO_from_Link(Seg,M->link_f,M->n_f,&grio);
			(*g_Rio) += grio;
		}
		if((*g_Rio)<INF){
			Bad_IO_from_Link_transition(Seg,M->link_s,M->n_s,choose_num,&grio);
			(*g_Rio) += grio;
		}
		if((*g_Rio)<INF){
			Bad_IO_from_Link(Seg,M->link_d,M->n_d,&grio);
			(*g_Rio) += grio;
		}
	}
	else
	{
		Bad_IO_from_Link(Seg,M->link_f,M->n_f,&grio);
		(*g_Rio) += grio;
		if((*g_Rio)<INF){
			Bad_IO_from_Link(Seg,M->link_s,M->n_s,&grio);
			(*g_Rio) += grio;
		}
		if((*g_Rio)<INF){
			Bad_IO_from_Link_transition(Seg,M->link_d,M->n_d,choose_num,&grio);
			(*g_Rio) += grio;
		}
	}
#else
	if (end_type == 0)
	{
		n_io+= Bad_IO_from_Link_transition(Seg,M->link_f,M->n_f,choose_num,g_Rio);
		if(n_io<INF) n_io+= Bad_IO_from_Link(Seg,M->link_s,M->n_s,g_Rio);
		if(n_io<INF) n_io+= Bad_IO_from_Link(Seg,M->link_d,M->n_d,g_Rio);
	}
	else if (end_type == 1)
	{
		n_io+= Bad_IO_from_Link(Seg,M->link_f,M->n_f,g_Rio);
		if(n_io<INF) n_io+= Bad_IO_from_Link_transition(Seg,M->link_s,M->n_s,choose_num,g_Rio);
		if(n_io<INF) n_io+= Bad_IO_from_Link(Seg,M->link_d,M->n_d,g_Rio);
	}
	else
	{
		n_io+= Bad_IO_from_Link(Seg,M->link_f,M->n_f,g_Rio);
		if(n_io<INF) n_io+= Bad_IO_from_Link(Seg,M->link_s,M->n_s,g_Rio);
		if(n_io<INF) n_io+= Bad_IO_from_Link_transition(Seg,M->link_d,M->n_d,choose_num,g_Rio);
	}
#endif	
	
	return n_io;
}


int Bad_IO_from_Link_connection_move(lineObj *Seg, LinkedList link, int n, lineObj **obj, int obj_num, double *g_Rio)
{
	int n_io = 0;
	Node *p;
	p = link;

	int c_x = Seg->x;
	int c_y = Seg->y;
	int c_len = Seg->len;
	double c_theta = Seg->theta;
	site c_enda = Seg->enda;
	site c_endb = Seg->endb;
	int flag;

	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;

		flag = 0;
		for (int j = 0;j < obj_num; j++)
		{
			if (l == obj[j])
			{
				flag = 1;
				break;
			}
		}
		if (!flag)	
		{
			int x= l->x;
			int y = l->y;
			int len = l->len;
			double theta = l->theta;
			site enda = l->enda;
			site endb = l->endb;

			int cret = 0;
			if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= (double)MAX(c_len,len)/2.*IO_RANGE_MARGIN)
				cret = 1;
			
			 if (l!= Seg && cret == 1)
			{
				if(l->enda.x == Seg->enda.x && l->enda.y == Seg->enda.y && l->endb.x == Seg->endb.x && l->endb.y == Seg->endb.y)
					n_io =4;
				if(l->enda.x == Seg->endb.x && l->enda.y == Seg->endb.y && l->endb.x == Seg->enda.x && l->endb.y == Seg->enda.y)
					n_io =4;	

				double tau_ij = MAX(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
				check_io(c_enda, c_endb, enda, endb, c_x, c_y, x, y, tau_ij, c_theta, theta, &n_io, g_Rio);
			}
		} // if(!flag)
	} // for(int i)
	return n_io;
}

int Bad_IO_connection_move(lineObj *Seg, Candy *M, lineObj **obj, int obj_num, double *g_Rio)
{
	int n_io = 0;
	(*g_Rio) = 0;
#ifdef QUALITY_CANDY
	*g_Rio+= Bad_IO_from_Link_connection_move(Seg,M->link_f,M->n_f,obj,obj_num,g_Rio);
	if(*g_Rio<INF) n_io+= Bad_IO_from_Link_connection_move(Seg,M->link_s,M->n_s,obj,obj_num,g_Rio);
	if(*g_Rio<INF) n_io+= Bad_IO_from_Link_connection_move(Seg,M->link_d,M->n_d,obj,obj_num,g_Rio);
#else
	n_io+= Bad_IO_from_Link_connection_move(Seg,M->link_f,M->n_f,obj,obj_num,g_Rio);
	if(n_io<INF) n_io+= Bad_IO_from_Link_connection_move(Seg,M->link_s,M->n_s,obj,obj_num,g_Rio);
	if(n_io<INF) n_io+= Bad_IO_from_Link_connection_move(Seg,M->link_d,M->n_d,obj,obj_num,g_Rio);
#endif	
	
	return n_io;
}

int Bad_IO_objects(lineObj **obj, int obj_num, double *g_Rio)
{
	int n_io = 0;

	for (int i = 0;i < obj_num; i++)
	{
		int c_x = obj[i]->x;
		int c_y = obj[i]->y;
		int c_len = obj[i]->len;
		double c_theta = obj[i]->theta;
		site c_enda = obj[i]->enda;
		site c_endb = obj[i]->endb;
		for (int j = 0;j < obj_num; j++){
			if (i == j) j++;
			if(j == obj_num) break;
			int x= obj[j]->x;
			int y = obj[j]->y;
			int len = obj[j]->len;
			double theta = obj[j]->theta;
			site enda = obj[j]->enda;
			site endb = obj[j]->endb;

			int cret = 0;
		if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= (double)MAX(c_len,len)/2.*IO_RANGE_MARGIN)
			cret = 1;
			
			if (cret == 1)
			{
				if(obj[j]->enda.x == obj[i]->enda.x && obj[j]->enda.y == obj[i]->enda.y 
					&& obj[j]->endb.x == obj[i]->endb.x && obj[j]->endb.y == obj[i]->endb.y)
					n_io =4;
				if(obj[j]->enda.x == obj[i]->endb.x && obj[j]->enda.y == obj[i]->endb.y 
					&& obj[j]->endb.x == obj[i]->enda.x && obj[j]->endb.y == obj[i]->enda.y)
					n_io =4;	

				double tau_ij = MAX(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
				check_io(c_enda, c_endb, enda, endb, c_x, c_y, x, y, tau_ij, c_theta, theta, &n_io, g_Rio);
			}
		} // for(int j)
	} // for(int i)
	return n_io;
}


int Bad_IO_from_Link_freeEndsC(lineObj *Seg, LinkedList link, int n, lineObj *obj , double *g_Rio)
{
	int n_io = 0;
	Node *p;
	p = link;

	int c_x = Seg->x;
	int c_y = Seg->y;
	int c_len = Seg->len;
	double c_theta = Seg->theta;
	site c_enda = Seg->enda;
	site c_endb = Seg->endb;

		for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
	//	if(l == obj)
	//	{
	//		printf("find obj to connect its end\n");
	//	}
		if (l != obj)
		{
		int x= l->x;
		int y = l->y;
		int len = l->len;
		double theta = l->theta;
		site enda = l->enda;
		site endb = l->endb;

		int cret = 0;
		if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) <= (double)MAX(c_len,len)/2.*IO_RANGE_MARGIN)
		cret = 1;
		
		 if (l!= Seg && cret == 1)
		{
			if(l->enda.x == Seg->enda.x && l->enda.y == Seg->enda.y && l->endb.x == Seg->endb.x && l->endb.y == Seg->endb.y)
				n_io =4;
			if(l->enda.x == Seg->endb.x && l->enda.y == Seg->endb.y && l->endb.x == Seg->enda.x && l->endb.y == Seg->enda.y)
				n_io =4;	

			double tau_ij = MAX(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
			check_io(c_enda, c_endb, enda, endb, c_x, c_y, x, y, tau_ij, c_theta, theta, &n_io, g_Rio);
		}
		}
		}
		return n_io;
}


int Bad_IO_freeEndsC(lineObj *Seg, lineObj *obj1, Candy *M, double *g_Rio)
{
	int n_io = 0;
	(*g_Rio) = 0;
#ifdef QUALITY_CANDY
	*g_Rio+= Bad_IO_from_Link_freeEndsC(Seg,M->link_f,M->n_f,obj1,g_Rio);
	if(*g_Rio<INF) n_io+= Bad_IO_from_Link_freeEndsC(Seg,M->link_s,M->n_s,obj1,g_Rio);
	if(*g_Rio<INF) n_io+= Bad_IO_from_Link_freeEndsC(Seg,M->link_d,M->n_d,obj1,g_Rio);
#else
	n_io+= Bad_IO_from_Link_freeEndsC(Seg,M->link_f,M->n_f,obj1,g_Rio);
	if(n_io<INF) n_io+= Bad_IO_from_Link_freeEndsC(Seg,M->link_s,M->n_s,obj1,g_Rio);
	if(n_io<INF) n_io+= Bad_IO_from_Link_freeEndsC(Seg,M->link_d,M->n_d,obj1,g_Rio);
#endif	
	
	return n_io;
}
