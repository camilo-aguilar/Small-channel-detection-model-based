#include "QualityCandy.h"

void check_eo(int c_x, int c_y, int x, int y, int c_len, int len,
			  double c_theta, double theta, double dist, int *n_eo, double *g_Rc)
{
#ifdef QUALITY_CANDY
	double g_tau, g_dist;
	if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) >= MAX(c_len,len)/2)
	{
		double tau_ij = MIN(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
		if (tau_ij < TAU_MAX)
			g_tau = -sigma(tau_ij*tau_ij,TAU_MAX*TAU_MAX);
		else
			g_tau = 2;
		g_dist = -sigma(dist*dist, NEIGHBOORHOOD*NEIGHBOORHOOD);
		//if (ABS(tau_ij-TAU_MAX) > _PI/20)
				*g_Rc += (g_tau+g_dist)/2.;

		// 임시로테스트
		if (tau_ij > TAU_MAX)
				(*n_eo)++;
	}
#else // quility candy
	if (sqrt(double((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))) >= MAX(c_len,len)/2)
	{
		double tau_ij = MIN(ABS(c_theta-theta),_PI-ABS(c_theta-theta));
		if (tau_ij > TAU_MAX)
		//if (ABS(tau_ij-TAU_MAX) > _PI/20)
				(*n_eo)++;
	}
#endif
}


int Bad_EO(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0.;
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_f,M->n_f, g_Rc);
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_s,M->n_s, g_Rc);
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_d,M->n_d, g_Rc);

	return n_eo;


}
#if 1
int Bad_EO_from_Link(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, double *g_Rc)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
	
	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;
#if 1
		dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo, &if_connect, &which_end_is_neighboor);
		if (l!= line && if_eo != 0)
		//if (if_eo != 0)
		{
			int len = l->len;
			int x = l->x;
			int y = l->y;
			if (which_end_is_neighboor ==1 && if_connect != 1) // only nearby, no connection
			{
				line->enda_L[line->enda_L_Num] = l;
				line->enda_L_Num++;
			}
			else if (which_end_is_neighboor == 2 && if_connect != 2)
			{
				line->endb_L[line->endb_L_Num] = l;
				line->endb_L_Num++;
			}
			else if (which_end_is_neighboor == 3)
				n_eo = 5;
			if(if_connect == 1)
			{
				line->enda_C[line->enda_C_Num] = l;
				line->enda_C_Num++;
			}
			if (if_connect == 2)
			{
				line->endb_C[line->endb_C_Num] = l;
				line->endb_C_Num++;

			}
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist,&n_eo,g_Rc);

		}
}
#endif
		return n_eo;
}
#endif

#if 1
int Bad_EO_from_Link_freeSeg(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, double *g_Rc)//, int *dn_f, int *dn_s, int *dn_d)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
//	(*dn_f) = 0;
//	(*dn_s) = 0;
//	(*dn_d) = 0;
	
	
	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;


		dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo,&if_connect, &which_end_is_neighboor);
		if (l!= line && if_eo != 0)
		//if (if_eo != 0)
		{
			int len = l->len;
			int x = l->x;
			int y = l->y;
			if (which_end_is_neighboor == 1 && if_connect != 1)  //only nearby, no connection
			{
				line->enda_L[line->enda_L_Num] = l;
				line->enda_L_Num++;
//				(*dn_f)++;

			}
			else if (which_end_is_neighboor == 2 && if_connect != 2)
			{
				line->endb_L[line->endb_L_Num] = l;
				line->endb_L_Num++;
				//line->enda_L[line->enda_L_Num] = l;
				//line->enda_L_Num++;
//				(*dn_f)++;
			}
//			else if (which_end_is_neighboor == 3 && if_connect != 3)  //only nearby, no connection
//			{
//				line->endb_L[line->endb_L_Num] = l;
//				line->endb_L_Num++;
////				(*dn_f)++;
//			}
//			else if (which_end_is_neighboor == 4 && if_connect != 4)
//			{
//				line->endb_L[line->endb_L_Num] = l;
//				line->endb_L_Num++;
////				(*dn_f)++;
//			}
			else if (which_end_is_neighboor == 3)
			{
				n_eo=5;
//				(*dn_f)++;
			}
			if(if_connect == 1)
			{
				line->enda_C[line->enda_C_Num] = l;
				line->enda_C_Num++;
//				if(l->type == 0) //free end
			}
			if (if_connect == 2)
			{
				line->endb_C[line->endb_C_Num] = l;
				line->endb_C_Num++;
//				line->enda_C[line->enda_C_Num] = l;
//				line->enda_C_Num++;

			}
			//if(if_connect == 3)
			//{
			//	line->endb_C[line->endb_C_Num] = l;
			//	line->endb_C_Num++;
			//}
			//if (if_connect == 4)
			//{
			//	line->endb_C[line->endb_C_Num] = l;
			//	line->endb_C_Num++;

			//}
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist,&n_eo, g_Rc);
		}
		}
		return n_eo;
}
#endif
int Bad_EO_freeSeg(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0;

	n_eo += Bad_EO_from_Link_freeSeg(line,neighborR,searchRatio,M->link_f,M->n_f, g_Rc);
	n_eo += Bad_EO_from_Link_freeSeg(line,neighborR,searchRatio,M->link_s,M->n_s, g_Rc);
	n_eo += Bad_EO_from_Link_freeSeg(line,neighborR,searchRatio,M->link_d,M->n_d, g_Rc);

	return n_eo;


}

int Bad_EO_death(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0;

	n_eo += Bad_EO_from_Link_death(line,neighborR,searchRatio,M->link_f,M->n_f, g_Rc);
	n_eo += Bad_EO_from_Link_death(line,neighborR,searchRatio,M->link_s,M->n_s, g_Rc);
	n_eo += Bad_EO_from_Link_death(line,neighborR,searchRatio,M->link_d,M->n_d, g_Rc);

	return n_eo;


}

#if 1
int Bad_EO_from_Link_death(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, double *g_Rc)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
	
		for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;

		dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo, &if_connect,&which_end_is_neighboor);
		if (l!= line && if_eo != 0)
		//if (if_eo != 0)
		{
			int len = l->len;
			int x = l->x;
			int y = l->y;
			
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist,&n_eo,g_Rc);
		}
		}
		return n_eo;
}

#endif

#if 1
int Bad_EO_from_Link_freeSeg_death(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, double *g_Rc)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
	
	
	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;

		dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo,&if_connect, &which_end_is_neighboor);
		if (l!= line && if_eo != 0)
		//if (if_eo != 0)
		{
			int len = l->len;
			int x = l->x;
			int y = l->y;
			
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist,&n_eo,g_Rc);
		}
		}

		return n_eo;
}
#endif
int Bad_EO_freeSeg_death(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0;

	n_eo += Bad_EO_from_Link_freeSeg_death(line,neighborR,searchRatio,M->link_f,M->n_f,g_Rc);
	n_eo += Bad_EO_from_Link_freeSeg_death(line,neighborR,searchRatio,M->link_s,M->n_s,g_Rc);
	n_eo += Bad_EO_from_Link_freeSeg_death(line,neighborR,searchRatio,M->link_d,M->n_d,g_Rc);

	return n_eo;


}

int Bad_EO_from_Link_transition(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, int choose_num, double *g_Rc)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
	
		for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;

		if (i!= choose_num)
		{
#if 1
			dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo, &if_connect,&which_end_is_neighboor);
		if (l!= line && if_eo != 0)
		//if (if_eo != 0)
		{
			int len = l->len;
			int x = l->x;
			int y = l->y;
			if (which_end_is_neighboor ==1 && if_connect != 1) // only nearby, no connection
			{
				line->enda_L[line->enda_L_Num] = l;
				line->enda_L_Num++;
			}
			else if (which_end_is_neighboor == 2 && if_connect != 2)
			{
				line->endb_L[line->endb_L_Num] = l;
				line->endb_L_Num++;
			}
			else if (which_end_is_neighboor == 3)
				n_eo = 5;
			if(if_connect == 1)
			{
				line->enda_C[line->enda_C_Num] = l;
				line->enda_C_Num++;
			}
			if (if_connect == 2)
			{
				line->endb_C[line->endb_C_Num] = l;
				line->endb_C_Num++;

			}
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist, &n_eo, g_Rc);
		}
#endif
		}
		}
		return n_eo;
}

int Bad_EO_transition(int neighborR,double searchRatio, lineObj *line, Candy *M, int choose_num, int end_type, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0;

	if (end_type == 0)
	{
	n_eo += Bad_EO_from_Link_transition(line,neighborR,searchRatio,M->link_f,M->n_f,choose_num,g_Rc);
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_s,M->n_s,g_Rc);
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_d,M->n_d,g_Rc);
	}
	else if (end_type == 1)
	{
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_f,M->n_f,g_Rc);
	n_eo += Bad_EO_from_Link_transition(line,neighborR,searchRatio,M->link_s,M->n_s,choose_num,g_Rc);
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_d,M->n_d,g_Rc);
	}
	else
	{
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_f,M->n_f,g_Rc);
	n_eo += Bad_EO_from_Link(line,neighborR,searchRatio,M->link_s,M->n_s,g_Rc);
	n_eo += Bad_EO_from_Link_transition(line,neighborR,searchRatio,M->link_d,M->n_d,choose_num,g_Rc);
	}

	return n_eo;


}

int Bad_EO_from_Link_connection_move(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, lineObj **obj, int obj_num, double *g_Rc)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	int flag;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
	
	for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;

		flag = 0;
		for (int j = 0;j < obj_num; j++){
			if (l == obj[j]){
				flag = 1;
				break;
			}
		}
		if (!flag){
#if 1
			dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo, &if_connect,&which_end_is_neighboor);
			if (l!= line && if_eo != 0)
			//if (if_eo != 0)
			{
				int len = l->len;
				int x = l->x;
				int y = l->y;
				if (which_end_is_neighboor ==1 && if_connect != 1) // only nearby, no connection
				{
					line->enda_L[line->enda_L_Num] = l;
					line->enda_L_Num++;
				}
				else if (which_end_is_neighboor == 2 && if_connect != 2)
				{
					line->endb_L[line->endb_L_Num] = l;
					line->endb_L_Num++;
				}
				else if (which_end_is_neighboor == 3)
					n_eo = 5;
				if(if_connect == 1)
				{
					line->enda_C[line->enda_C_Num] = l;
					line->enda_C_Num++;
				}
				if (if_connect == 2)
				{
					line->endb_C[line->endb_C_Num] = l;
					line->endb_C_Num++;

				}
				check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist,&n_eo,g_Rc);
			}
#endif
		}
	}
	return n_eo;
}

int Bad_EO_connection_move(int neighborR,double searchRatio, lineObj *line, Candy *M, lineObj **obj, int obj_num, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0;

	n_eo += Bad_EO_from_Link_connection_move(line,neighborR,searchRatio,M->link_f,M->n_f,obj,obj_num,g_Rc);
	n_eo += Bad_EO_from_Link_connection_move(line,neighborR,searchRatio,M->link_s,M->n_s,obj,obj_num,g_Rc);
	n_eo += Bad_EO_from_Link_connection_move(line,neighborR,searchRatio,M->link_d,M->n_d,obj,obj_num,g_Rc);

	return n_eo;
}


int Bad_EO_objects(int neighboorR,double searchRatio, lineObj **obj, int obj_num, double *g_Rc)
{
	int n_eo = 0;
	double dist;
	(*g_Rc) = 0;
	
	for (int i = 0;i < obj_num; i++)
	{
		int c_x = obj[i]->x;
		int c_y = obj[i]->y;
		site c_enda = obj[i]->enda;
		site c_endb = obj[i]->endb;
		double c_theta = obj[i]->theta;
		int c_len=obj[i]->len;
		for (int j = 0;j < obj_num; j++){
			if (i == j) j++;
			if(j == obj_num) break;
			site enda = obj[j]->enda;
			site endb = obj[j]->endb;

			int if_eo = 0;
			int which_end_is_neighboor = 0;
			int if_connect = 0;
#if 1
			dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo, &if_connect,&which_end_is_neighboor);
			if (if_eo != 0)
			//if (if_eo != 0)
			{
				int len = obj[j]->len;
				int x = obj[j]->x;
				int y = obj[j]->y;
				if (which_end_is_neighboor ==1 && if_connect != 1) // only nearby, no connection
				{
					obj[i]->enda_L[obj[i]->enda_L_Num] = obj[j];
					obj[i]->enda_L_Num++;
				}
				else if (which_end_is_neighboor == 2 && if_connect != 2)
				{
					obj[i]->endb_L[obj[i]->endb_L_Num] = obj[j];
					obj[i]->endb_L_Num++;
				}
				else if (which_end_is_neighboor == 3)
					n_eo = 5;
				if(if_connect == 1)
				{
					obj[i]->enda_C[obj[i]->enda_C_Num] = obj[j];
					obj[i]->enda_C_Num++;
				}
				if (if_connect == 2)
				{
					obj[i]->endb_C[obj[i]->endb_C_Num] = obj[j];
					obj[i]->endb_C_Num++;

				}
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, obj[j]->theta, dist,&n_eo, g_Rc);
			}
#endif
		}
	}
	return n_eo;
}


int Bad_EO_from_Link_freeEndsC(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, lineObj *obj1, double *g_Rc)
{
	int c_x = line->x;
	int c_y = line->y;
	site c_enda = line->enda;
	site c_endb = line->endb;
	double c_theta = line->theta;
	int c_len=line->len;
	double dist;

	int n_eo = 0;
	Node *p;
	p = link;
	
		for (int i = 0;i < n; i++)
	{
		p = p->next;
		lineObj *l = p->index;
		site enda = l->enda;
		site endb = l->endb;

		int if_eo = 0;
		int which_end_is_neighboor = 0;
		int if_connect = 0;

	//	if(l == obj1)
	//	{
	//		printf("find obj to connect its end\n ");
	//	}
		if (l != obj1)
		{
#if 1
			dist = CheckConnection(c_enda,c_endb,enda,endb,searchRatio,neighboorR,&if_eo, &if_connect,&which_end_is_neighboor);
		if (l!= line && if_eo != 0)
		//if (if_eo != 0)
		{
			int len = l->len;
			int x = l->x;
			int y = l->y;
			if (which_end_is_neighboor ==1 && if_connect != 1) // only nearby, no connection
			{
				line->enda_L[line->enda_L_Num] = l;
				line->enda_L_Num++;
			}
			else if (which_end_is_neighboor == 2 && if_connect != 2)
			{
				line->endb_L[line->endb_L_Num] = l;
				line->endb_L_Num++;
			}
			else if (which_end_is_neighboor == 3)
				n_eo = 5;
			if(if_connect == 1)
			{
				line->enda_C[line->enda_C_Num] = l;
				line->enda_C_Num++;
			}
			if (if_connect == 2)
			{
				line->endb_C[line->endb_C_Num] = l;
				line->endb_C_Num++;

			}
			check_eo(c_x, c_y, x, y, c_len, len, c_theta, l->theta, dist,&n_eo,g_Rc);
		}
#endif
		}
		}
		return n_eo;
}


int Bad_EO_freeEndsC(int neighborR,double searchRatio, lineObj *line, Candy *M, lineObj *obj1, double *g_Rc)
{
	int n_eo = 0;
	(*g_Rc) = 0;

	n_eo += Bad_EO_from_Link_freeEndsC(line,neighborR,searchRatio,M->link_f,M->n_f,obj1,g_Rc);
	n_eo += Bad_EO_from_Link_freeEndsC(line,neighborR,searchRatio,M->link_s,M->n_s,obj1,g_Rc);
	n_eo += Bad_EO_from_Link_freeEndsC(line,neighborR,searchRatio,M->link_d,M->n_d,obj1,g_Rc);

	return n_eo;
}

