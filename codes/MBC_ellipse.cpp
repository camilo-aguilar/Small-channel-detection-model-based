/* Created by Dae Woo Kim 01/05/15 */

#include<iostream>
#include<time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#include "tiff.h"
#include "allocate.h"
#include "util.h"
#include "common.h"

#include "randlib.h"
#include "typeutil.h"
#ifdef _WIN64
	#include <cstdio>
	#include <opencv2\opencv.hpp>
#else
	#include "cv.h"
	#include "cxcore.h"
	#include "highgui.h"
#endif
#include "graph.h"
#include "MBD.h"
#include "MBD_ellipse.h"
#include "MBD_ellipsoid.h"
#include "ellipsoid.h"

using namespace std;

extern double *cosine;
extern double *sine;
extern double *cosine2;
extern double *sine2;


LinkedList3D LinkedList3DInit(void)
{
    Node3D *L;
    L = (Node3D *)malloc(sizeof(Node3D));  
    if(L == NULL)                     
        printf("fail to allocate node\n");
    L->next = NULL;   

	return L;
}

void LinkedList3DInsert(LinkedList3D L,int i, ellipsoidObj* x)
{
    Node3D *pre;                     
    pre = L;
    int tempi = 0;
    for (tempi = 1; tempi < i; tempi++)
        pre = pre->next;                
    Node3D *p;                                
    p = (Node3D *)malloc(sizeof(Node3D));
    p->index = x; 
    p->next = pre->next;
    pre->next = p;
     
   // return L;                           
} 

void LinkedList3DSortInsert(LinkedList3D L,int k, ellipsoidObj* x)
{
    int i;
    Node3D *pre;
	ellipsoidObj* obj;
    pre = L;
	if(k>1){
		for(i=1;i<k;i++){
			obj = pre->next->index;
			if(x->dataterm >= obj->dataterm){
				break;
			}
			pre = pre->next;
		}
	}
    Node3D *p;                                
    p = (Node3D *)malloc(sizeof(Node3D));
    p->index = x; 
    p->next = pre->next;
    pre->next = p;
     
   // return L;                           
} 

void LinkedList3DDelete(LinkedList3D L,ellipsoidObj* x)
{
    Node3D *p,*pre;
	pre = L;
    p = L->next;
    while(p->index != x)              
    {   
        pre = p; 
        p = p->next;
    }
    pre->next = p->next;        
	free(p);
} 

void LinkedList3DDelete(LinkedList3D L, int k)
{
	int i;
    Node3D *p,*pre;
	pre = L;
    p = L->next;
	for(i=1;i<k;i++){             
        pre = p; 
		p = p->next;
	}
    pre->next = p->next;        
	free(p); // node p is deleted, object p->index will not be deleted
} 

Node3D *LinkedList3DGetObj(LinkedList3D L, int k)
{
	int i;
    Node3D *p;
    p = L->next;
	for(i=1;i<k;i++)             
		p = p->next;
	return p; 
} 

Node3D *LinkedList3DGetNode(LinkedList3D L, ellipsoidObj *x)
{
    Node3D *p;
    p = L->next;
    while(p->index != x)              
		p = p->next;
	return p; 
}

#if 1
ElipsoidConfig *ElipsoidConfigInit(MPP_Parameters par_elsoid)
{
	ElipsoidConfig *config;
	int k;
	double n = par_elsoid.n;
	double cosine_u[POINT_NUM2+1][POINT_NUM],sine_u[POINT_NUM2+1][POINT_NUM];
	double cosine_v[POINT_NUM2+1],sine_v[POINT_NUM2+1];
	double cosv2_n[POINT_NUM2+1],sinv2_n[POINT_NUM2+1];
	double cosu2_n[POINT_NUM2+1][POINT_NUM],sinu2_n[POINT_NUM2+1][POINT_NUM];
	double cosv2mn_n[POINT_NUM2+1],sinv2mn_n[POINT_NUM2+1];
	double cosu2mn_n[POINT_NUM2+1][POINT_NUM],sinu2mn_n[POINT_NUM2+1][POINT_NUM];
	double dxdv[POINT_NUM2+1],dydv[POINT_NUM2+1];
	double dxdu[POINT_NUM2+1][POINT_NUM],dydu[POINT_NUM2+1][POINT_NUM];

	config = (ElipsoidConfig *)malloc(sizeof(ElipsoidConfig));
	config->elink = LinkedList3DInit();
	config->mp_num = 0;

	/* Allocate and initialize memory	*/
	config->cosine = (double *)mget_spc(POINT_NUM,sizeof(double));
	config->sine = (double *)mget_spc(POINT_NUM,sizeof(double));
	config->cosine2 = (double *)mget_spc(POINT_NUM3,sizeof(double));
	config->sine2 = (double *)mget_spc(POINT_NUM3,sizeof(double));
	config->coscos = (double **)get_img(POINT_NUM,POINT_NUM,sizeof(double));
	config->cossin = (double **)get_img(POINT_NUM,POINT_NUM,sizeof(double));
	config->coscos2 = (double **)get_img(POINT_NUM3,POINT_NUM3,sizeof(double));
	config->cossin2 = (double **)get_img(POINT_NUM3,POINT_NUM3,sizeof(double));
//	config->nx = (double *)mget_spc(POINT_NUM,sizeof(double));
//	config->ny = (double *)mget_spc(POINT_NUM,sizeof(double));
	config->num_ufrm = (int *)mget_spc(POINT_NUM2+1,sizeof(int));
	config->sine_ufrm = (double *)mget_spc(POINT_NUM2+1,sizeof(double));
	config->coscos_ufrm = (double **)get_img(POINT_NUM,POINT_NUM2+1,sizeof(double));
	config->cossin_ufrm = (double **)get_img(POINT_NUM,POINT_NUM2+1,sizeof(double));
	config->nx3d = (double **)get_img(POINT_NUM,POINT_NUM2+1,sizeof(double));
	config->ny3d = (double **)get_img(POINT_NUM,POINT_NUM2+1,sizeof(double));
	config->nz3d = (double **)get_img(POINT_NUM,POINT_NUM2+1,sizeof(double));

	k = POINT_NUM/POINT_NUM3;
	for (int i = 0;i<POINT_NUM;i++)
	{
		config->cosine[i]=SIGN(cos(i*M_PI/POINT_NUM2))*pow((cos(i*M_PI/POINT_NUM2)*cos(i*M_PI/POINT_NUM2)),1/n);
		config->sine[i]=SIGN(sin(i*M_PI/POINT_NUM2))*pow((sin(i*M_PI/POINT_NUM2)*sin(i*M_PI/POINT_NUM2)),1/n);
		if(i%k==0)
		{
			config->cosine2[i/k]=config->cosine[i];
			config->sine2[i/k]=config->sine[i];
		}
	}
	for (int i = 0;i<POINT_NUM;i++)
		for (int j = 0;j<POINT_NUM;j++)
		{
			config->coscos[i][j]=config->cosine[i]*config->cosine[j];
			config->cossin[i][j]=config->cosine[i]*config->sine[j];
			if((i%k==0)&&(j%k==0))
			{
				config->coscos2[i/k][j/k]=config->coscos[i][j];
				config->cossin2[i/k][j/k]=config->cossin[i][j];
			}
		}

	// for surface normal
	double u1, u2, u3=0.;
	double v1, v2, v3;
	config->total_num_ufrm = 0;
	for (int i = 0;i<=POINT_NUM2;i++){ // 0~30 : -pi/2 ~ pi/2
		double v = (double)(i-POINT_NUM4)*M_PI/POINT_NUM2;//-(i-POINT_NUM4)*(2pi/POINT_NUM)
		config->num_ufrm[i] = POINT_NUM*cos(v); 
		cosine_v[i] = cos(v);
		sine_v[i] = sin(v);
		cosv2_n[i] = SIGN(cosine_v[i])*pow(fabs(cosine_v[i]),2./n);
		sinv2_n[i] = SIGN(sine_v[i])*pow(fabs(sine_v[i]),2./n);
		config->sine_ufrm[i] = sinv2_n[i];
		if(cosine_v[i]==0)
			cosv2mn_n[i] = INF;
		else
			cosv2mn_n[i] = pow(fabs(cosine_v[i]),(2.-n)/n);
		if(sine_v[i]==0)
			sinv2mn_n[i] = INF;
		else
			sinv2mn_n[i] = pow(fabs(sine_v[i]),(2.-n)/n);
		dxdv[i] = -cosv2mn_n[i]*sine_v[i];
		dydv[i] = sinv2mn_n[i]*cosine_v[i];
		v3 = dydv[i];
		double du = 2*M_PI/config->num_ufrm[i];
		for (int j = 0;j<config->num_ufrm[i];j++){
			config->total_num_ufrm++;
			double u = (double)j*du;
			cosine_u[i][j] = cos(u);
			sine_u[i][j] = sin(u);
			cosu2_n[i][j] = SIGN(cosine_u[i][j])*pow(fabs(cosine_u[i][j]),2/n);
			sinu2_n[i][j] = SIGN(sine_u[i][j])*pow(fabs(sine_u[i][j]),2/n);
			if(cosine_u[i][j]==0)
				cosu2mn_n[i][j] = INF;
			else
				cosu2mn_n[i][j] = pow(fabs(cosine_u[i][j]),(2.-n)/n);
			if(sine_u[i][j]==0)
				sinu2mn_n[i][j] = INF;
			else
				sinu2mn_n[i][j] = pow(fabs(sine_u[i][j]),(2.-n)/n);
			config->coscos_ufrm[i][j] = cosv2_n[i]*cosu2_n[i][j];
			config->cossin_ufrm[i][j] = cosv2_n[i]*sinu2_n[i][j];
			dxdu[i][j] = -cosu2mn_n[i][j]*sine_u[i][j];
			dydu[i][j] = sinu2mn_n[i][j]*cosine_u[i][j];
			u1 = cosv2_n[i]*dxdu[i][j];
			u2 = cosv2_n[i]*dydu[i][j];
			v1 = dxdv[i]*cosu2_n[i][j];
			v2 = dxdv[i]*sinu2_n[i][j];
			config->nx3d[i][j] = u2*v3;
			config->ny3d[i][j] = -u1*v3;
			config->nz3d[i][j] = u1*v2-u2*v1;
		}
	}

	// for convergence info
	config->record_step = 100000;
	config->exe_time[0] = 0;
	config->num_obj[0] = 0;
	config->energy[0] = 0;
		
	return config;
}

#endif


void free_obj_in_ElipsoidConfig(ElipsoidConfig *config)
{
	Node3D *p;

	while(config->mp_num>0){
		p = config->elink->next;
		if(p->index->ptmpmask!=NULL){
			free_volume((void***)(p->index->ptmpmask));
			p->index->ptmpmask = NULL;
		}
		free(p->index);
		config->elink->next = p->next;
		free(p);
		config->mp_num--;
	}

}


void freeElipsoidConfig(ElipsoidConfig *config)
{
	free_obj_in_ElipsoidConfig(config);
	free((void*)config->cosine);
	free((void*)config->sine);
	free((void*)config->cosine2);
	free((void*)config->sine2);
	free_img((void**)config->nx3d);
	free_img((void**)config->ny3d);
	free_img((void**)config->nz3d);
	free_img((void**)config->coscos);
	free_img((void**)config->cossin);
	free_img((void**)config->coscos2);
	free_img((void**)config->cossin2);
	free((void*)config->num_ufrm);
	free((void*)config->sine_ufrm);
	free_img((void**)config->coscos_ufrm);
	free_img((void**)config->cossin_ufrm);
	free(config->elink);
	free(config);
}

void Qdelete_elipsoid3(ElipsoidConfig *config, int k)
{
	Node3D *p;
	p = LinkedList3DGetObj(config->elink, k);
	if (p->index->ptmpmask != NULL)
		free_volume((void***)(p->index->ptmpmask));
	p->index->ptmpmask = NULL;
	free(p->index);
	LinkedList3DDelete(config->elink, k);
	config->mp_num--;
	// LinkedListDelete(config->elink, p); // this also works
}

void Qdelete_elipsoid3(LinkedList3D L, int k, int *num)
{
	Node3D *p;
	p = LinkedList3DGetObj(L, k);
	if (p->index->ptmpmask != NULL)
		free_volume((void***)(p->index->ptmpmask));
	p->index->ptmpmask = NULL;
	free(p->index);
	LinkedList3DDelete(L, k);
	(*num)--;
	// LinkedListDelete(config->elink, p); // this also works
}

#define ELLIPSE_OVERLAP_TH 0
int elips_overlap_check(int y1,int x1, int y2, int x2, TmpmaskSt **tmpmaskarry)
{
	int sx1, sy1, sx2, sy2;
	int minx1, miny1, minx2, miny2;
	int maxx1, maxy1, maxx2, maxy2;
	int minii, minjj;
	int maxii, maxjj;
	int overlap;
	double min_area, ratio = 0.;

	sx1 = tmpmaskarry[y1][x1].sizex/2;
	sy1 = tmpmaskarry[y1][x1].sizey/2;
	sx2 = tmpmaskarry[y2][x2].sizex/2;
	sy2 = tmpmaskarry[y2][x2].sizey/2;

	if((sx1+sx2>abs(x2-x1))&&(sy1+sy2>abs(y2-y1))&&((x1!=x2)||(y1!=y2))){ 
		overlap = 0;
		minx1 = x1-sx1;
		minx2 = x2-sx2;
		miny1 = y1-sy1;
		miny2 = y2-sy2;
		maxx1 = x1+sx1;
		maxx2 = x2+sx2;
		maxy1 = y1+sy1;
		maxy2 = y2+sy2;
		minjj = max(minx1,minx2)-minx1;
		maxjj = min(maxx1,maxx2)-minx1;
		minii = max(miny1,miny2)-miny1;
		maxii = min(maxy1,maxy2)-miny1;
		for (int ii = minii;ii<maxii;ii++)
			for (int jj = minjj;jj<maxjj;jj++)
			{
				if((tmpmaskarry[y1][x1].ptmpmask[ii][jj]!=0)&&(tmpmaskarry[y2][x2].ptmpmask[ii+miny1-miny2][jj+minx1-minx2]!=0))
					overlap ++;
			}
		min_area = (double)min(tmpmaskarry[y1][x1].area,tmpmaskarry[y2][x2].area);
		ratio = (double)overlap/min_area;
		if (ratio>ELLIPSE_OVERLAP_TH)
			return 1;
	}
	return 0;
}

int elips_overlap_check3(ellipseObj *obj1, ellipseObj *obj2)
{
	int x1 = obj1->x;
	int y1 = obj1->y;
	int x2 = obj2->x;
	int y2 = obj2->y;
	int sx1 = obj1->sizex/2;
	int sy1 = obj1->sizey/2;
	int sx2 = obj2->sizex/2;
	int sy2 = obj2->sizey/2;
	int overlap;
	if(obj1->z==obj2->z){
		if((sx1+sx2>abs(x2-x1))&&(sy1+sy2>abs(y2-y1))){ 
			overlap = 0;
			int minx1 = x1-sx1;
			int minx2 = x2-sx2;
			int miny1 = y1-sy1;
			int miny2 = y2-sy2;
			int maxx1 = x1+sx1;
			int maxx2 = x2+sx2;
			int maxy1 = y1+sy1;
			int maxy2 = y2+sy2;
			int minjj = max(minx1,minx2)-minx1;
			int maxjj = min(maxx1,maxx2)-minx1;
			int minii = max(miny1,miny2)-miny1;
			int maxii = min(maxy1,maxy2)-miny1;
			for (int ii = minii;ii<maxii;ii++)
				for (int jj = minjj;jj<maxjj;jj++)
				{
					if((obj1->ptmpmask[ii][jj]!=0)&&(obj2->ptmpmask[ii+miny1-miny2][jj+minx1-minx2]!=0))
						overlap ++;
				}
			double	min_area = (double)min(obj1->area,obj2->area);
			double ratio = (double)overlap/min_area;
			if (ratio>ELLIPSE_OVERLAP_TH)
				return 1; // overlapped
		}
	}
	return 0; // not overlapped
}

int elips_overlap_check_all(ellipseObj *obj, LinkedList link, int num)
{
	int overlap;
	Node *p;
	p = link;

	for(int i=1;i<=num;i++){
		p = p->next;
		if((obj->z==p->index->z)&&(obj!=p->index)){
			overlap = elips_overlap_check3(obj, p->index);
			if(overlap)
				return 1;
		}
	}
	return 0;
}

#define ELLIPSOID_OVERLAP_TH 0
#define ELLIPSOID_CY_OVERLAP_TH 0.15

int elipsoid_overlap_check3(ellipsoidObj *obj1, ellipsoidObj *obj2)
{
	int x1 = obj1->x;
	int y1 = obj1->y;
	int z1 = obj1->z;
	int x2 = obj2->x;
	int y2 = obj2->y;
	int z2 = obj2->z;
	int sx1 = obj1->sizex/2;
	int sy1 = obj1->sizey/2;
	int sz1 = obj1->sizez/2;
	int sx2 = obj2->sizex/2;
	int sy2 = obj2->sizey/2;
	int sz2 = obj2->sizez/2;
	int overlap;
	if((sx1+sx2>abs(x2-x1))&&(sy1+sy2>abs(y2-y1))&&(sz1+sz2>abs(z2-z1))){ 
		overlap = 0;
		int minx1 = x1-sx1;
		int minx2 = x2-sx2;
		int miny1 = y1-sy1;
		int miny2 = y2-sy2;
		int minz1 = z1-sz1;
		int minz2 = z2-sz2;
		int maxx1 = x1+sx1;
		int maxx2 = x2+sx2;
		int maxy1 = y1+sy1;
		int maxy2 = y2+sy2;
		int maxz1 = z1+sz1;
		int maxz2 = z2+sz2;
		int minjj = max(minx1,minx2)-minx1;
		int maxjj = min(maxx1,maxx2)-minx1;
		int minii = max(miny1,miny2)-miny1;
		int maxii = min(maxy1,maxy2)-miny1;
		int minkk = max(minz1,minz2)-minz1;
		int maxkk = min(maxz1,maxz2)-minz1;
		for (int kk = minkk;kk<maxkk;kk++)
			for (int ii = minii;ii<maxii;ii++)
				for (int jj = minjj;jj<maxjj;jj++)
				{
					if((obj1->ptmpmask[kk][ii][jj]!=0)&&(obj2->ptmpmask[kk+minz1-minz2][ii+miny1-miny2][jj+minx1-minx2]!=0))
						overlap ++;
				}
		double	min_volume = (double)min(obj1->volume,obj2->volume);
		double ratio = (double)overlap/min_volume;
		if(obj1->type == OBJ3D_ELLIPSOID){
			if (ratio>ELLIPSOID_OVERLAP_TH)
				return 1; // overlapped
		}
		else{
			if (ratio>ELLIPSOID_CY_OVERLAP_TH)
				return 1; // overlapped
		}
	}
	return 0; // not overlapped
}

int elipsoid_overlap_check_all(ellipsoidObj *obj, LinkedList3D link, int num)
{
	int overlap;
	Node3D *p;
	p = link;

	for(int i=1;i<=num;i++){
		p = p->next;
		overlap = elipsoid_overlap_check3(obj, p->index);
		if(overlap)
			return 1;
	}
	return 0;
}

int elipsoid_overlap_check_disk(ellipsoidObj *obj1, ellipsoidObj *obj2)
{
	double x,y,z,dist;
	x = (obj1->x - obj2->x);
	y = (obj1->y - obj2->y);
	z = (obj1->z - obj2->z);
	dist = x*x+y*y+z*z;
	if(dist<10)
		return 1;
	else
		return 0;
}

int elipsoid_overlap_check_all_disk(ellipsoidObj *obj, LinkedList3D link, int num)
{
	int overlap;
	Node3D *p;
	p = link;

	for(int i=1;i<=num;i++){
		p = p->next;
		overlap = elipsoid_overlap_check_disk(obj, p->index);
		if(overlap)
			return 1;
	}
	return 0;
}

/****************************************************/
/* sorted marked point list                         */
/* (x,y): coordinate, l:single energy               */
/* Qmark[0][0] : number of marked point             */
/****************************************************/
void MBC_Qadd(int **Qmark, int y, int x, int l, int live)
{
	int flag = 0;
	Qmark[0][0]++;

	for (int i = 0;i<Qmark[0][0];i++)
	{
		if (l>=Qmark[i][2])
		{
			flag = 1;
			for (int j = Qmark[0][0];j>i;j--)
			{
				Qmark[j][0]= Qmark[j-1][0];
				Qmark[j][1]= Qmark[j-1][1];
				Qmark[j][2]= Qmark[j-1][2];
				Qmark[j][3]= Qmark[j-1][3];

			}
			Qmark[i][0] = y;
			Qmark[i][1] = x;
			Qmark[i][2] = l;
			Qmark[i][3] = live;
			break;
		}
	}
	if (flag == 0)
	{
		int i = Qmark[0][0];
		Qmark[i][0] = y;
		Qmark[i][1] = x;
		Qmark[i][2] = l;
		Qmark[i][3] = live;
	}
}

void MBC_Qdelete(int **Qmark, TmpmaskSt **tmpmaskarry, int k)
{
	int y = Qmark[k][0];
	int x = Qmark[k][1];

	free_img((void**)tmpmaskarry[y][x].ptmpmask);
	tmpmaskarry[y][x].ptmpmask = NULL;

	for (int j = k;j<Qmark[0][0];j++)
		{
			Qmark[j][0] = Qmark[j+1][0];
			Qmark[j][1] = Qmark[j+1][1];
			Qmark[j][2] = Qmark[j+1][2];
			Qmark[j][3] = Qmark[j+1][3];
	}

	Qmark[0][0]--;
}

void graph_cut(int **Qmark1, int **Qmark2, TmpmaskSt **tmpmaskarry)
{
///////////////////////////////////////////////////
//
//		        SOURCE
//		       /       \
//		     1/         \2
//		     /      3    \
//		   node0 -----> node1
//		     |   <-----   |
//		     |      4     |
//		     \            /
//		     5\          /6
//		       \        /
//		          SINK
//
///////////////////////////////////////////////////
//	typedef Graph<int,int,int> GraphType;
//	GraphType *g = new GraphType(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1); 
//	g -> add_node(); 
//	g -> add_node(); 

//	g -> add_tweights( 0,   /* capacities */  1, 5 );
//	g -> add_tweights( 1,   /* capacities */  2, 6 );
//	g -> add_edge( 0, 1,    /* capacities */  3, 4 );

//	int flow = g -> maxflow();
//
//	printf("Flow = %d\n", flow);
//	printf("Minimum cut:\n");
//	if (g->what_segment(0) == GraphType::SOURCE)
//		printf("node0 is in the SOURCE set\n");
//	else
//		printf("node0 is in the SINK set\n");
//	if (g->what_segment(1) == GraphType::SOURCE)
//		printf("node1 is in the SOURCE set\n");
//	else
//		printf("node1 is in the SINK set\n");
//	delete g;

	typedef Graph<int,int,int> GraphType;
	int num_nodes = Qmark1[0][0]+Qmark2[0][0];
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ Qmark2[0][0]*4); 

	for(int i = 0;i<num_nodes;i++)
		g -> add_node(); 
	for(int i = 0;i<Qmark1[0][0];i++){
		g -> add_tweights( i, Qmark1[i+1][2],  1000-Qmark1[i+1][2]);
	}
	for(int i = 0;i<Qmark2[0][0];i++){
		g -> add_tweights( i+Qmark1[0][0], 1000-Qmark2[i+1][2], Qmark2[i+1][2] );
	}
	for(int i = 0;i<Qmark1[0][0];i++){
		for(int j = 0;j<Qmark2[0][0];j++){
			int y1 = Qmark1[i+1][0];
			int x1 = Qmark1[i+1][1];
			int y2 = Qmark2[j+1][0];
			int x2 = Qmark2[j+1][1];
			int overlap = elips_overlap_check(y1,x1,y2,x2,tmpmaskarry);
			if(overlap)
				g -> add_edge( i, j+Qmark1[0][0], 0, INF );
		}
	}
//	printf("maxflow alorithm..\n");
	int flow = g -> maxflow();
	for(int i = 0;i<Qmark1[0][0];i++){
		if (g->what_segment(i) == GraphType::SOURCE)
			Qmark1[i+1][3] = 0;
		else
			Qmark1[i+1][3] = 1;
	}
	for(int i = 0;i<Qmark2[0][0];i++){
		if (g->what_segment(i+Qmark1[0][0]) == GraphType::SOURCE)
			Qmark2[i+1][3] = 1;
		else
			Qmark2[i+1][3] = 0;
	}
}

double MBC_ellipseMain(unsigned char **img, int **mask, double **Pmask, double **outP, Config_Ellipses *config, int width, int height,
					 double *cosine, double *sine, double *cosine2, double *sine2, MPP_Parameters mpp)
{
	int max_a  = mpp.max_a;
	int max_b = mpp.max_b;
	int min_a  = mpp.min_a;
	int min_b = mpp.min_b;
	int new_a, new_b;
	double n = mpp.n;

	double pi = CV_PI;

	double max_theta1 = mpp.max_theta;
	double min_theta1 = mpp.min_theta;
	double new_theta;

	double pen  = mpp.pen;		// overlap penalty
	double thredh = mpp.thredh; // Eq. (5) T

	double b_zero = mpp.b_zero;
	double delta = mpp.delta;	// 0.9 sigma for death step 
	double beta = mpp.beta_mpp;	// 10  alpha for death step
	double F = mpp.F;			// sigma and alpha decreasing and increasing factor
	double l_th = mpp.l_th;
	int iter_num = mpp.iter_num;
	int overlap, birth;

	struct TIFF_img inter_img;
	int tag = 0;
	int select = 0;
	int select2 = 0;
	clock_t total_time = 0;
	get_TIFF ( &inter_img, height, width, 'g' );

	// 0(y),1(x),2(likelihood),3(alive/dead)
	int **Qmark1 = (int **)get_img(4,width*height,sizeof(int));
	int **Qmark2 = (int **)get_img(4,width*height,sizeof(int));
	int tmpmask_len = 2*(max_a+9);
	int **tmpmask=(int **)get_img(tmpmask_len,tmpmask_len,sizeof(int));
//	int **mask=(int **)get_img(width,height,sizeof(int));

	double ***img_mark = (double ***)malloc(height*sizeof(double **));   
    for(int i=0; i<height; i++)   
    {     
        img_mark[i] = (double**)malloc(width*sizeof(double*));   
        for(int j=0; j<width; j++)   
		{   // 0(mark exist),1(a),2(b),3(theta),4(likelihood)
            img_mark[i][j] = (double*)malloc(5*sizeof(double));   
        }   
    }
	TmpmaskSt **tmpmaskarry = (TmpmaskSt **)get_img(width,height,sizeof(TmpmaskSt));   

	for (int i =0;i<width*height;i++)
		for(int j = 0;j<4;j++){
			Qmark1[i][j] = 0;
			Qmark2[i][j] = 0;
		}

	for (int i = 0;i<height;i++)
		for (int j = 0;j<width;j++)
			for (int k = 0;k<5;k++){
				img_mark[i][j][k] = 0;
			}

    for(int i=0; i<height; i++) {
        for(int j=0; j<width; j++) {   
            tmpmaskarry[i][j].ptmpmask = NULL; 
        }   
    }       

	//set the end tag for Qmark
	Qmark1[0][0] = 0;  // store the length of the queue
	Qmark1[0][1] = 0;
	Qmark1[0][2] = 10000; // max likelihood

	for (int iter = 1;iter < iter_num; iter++)
	{
		printf("iter = %d\n",iter);
		Qmark2[0][0] = 0;  // store the length of the queue
		Qmark2[0][1] = 0;
		Qmark2[0][2] = 10000; // max likelihood
		// birth step
		double p = b_zero*delta;

		clock_t start_time=clock();

		for (int i = 4;i<height-5;i++)
			for(int j = 4;j<width-5;j++)
			{
				if (img_mark[i][j][0]!=1)
				{
					double pp = rand()/(double)(RAND_MAX);

					if (pp < p*Pmask[i][j])
					{
						

						new_a = (int)floor((rand()/(double)(RAND_MAX)*(max_a-min_a)+min_a)+0.5);
						new_b = (int)floor((rand()/(double)(RAND_MAX)*(max_b-min_b)+min_b)+0.5);
						new_theta = rand()/(double)(RAND_MAX)*(max_theta1-min_theta1)+min_theta1;
						
						for (int ii = 0;ii<tmpmask_len;ii++)
							for (int jj = 0;jj<tmpmask_len;jj++)
							{
								tmpmask[ii][jj] = 0;
							}					

						double d;
						double l = ellipse_likelyhood(tmpmask,tmpmask_len,img,width,height,i,j,new_a,new_b,new_theta,n,&d,thredh,cosine,sine);

						if (l < l_th)//-0.55)
						{
							int likelihood = (int)(500.*(1+l));

							if(iter==1){
								overlap = 0;
								birth = 1;
								drawSuperEllipse2(tmpmaskarry, i, j, new_a, new_b,new_theta,n,1,tmpmask_len,cosine2,sine2);
								if(Qmark1[0][0]!=0){
									for(int m=1;m<=Qmark1[0][0];m++){
										int y = Qmark1[m][0];
										int x = Qmark1[m][1];
										overlap = elips_overlap_check(i,j,y,x,tmpmaskarry);
										if(overlap){
											birth = 0;
											break;
										}
									}
								}
								if(birth){
									img_mark[i][j][0] = 1;   //give birth to this point
									img_mark[i][j][1]=new_a;
									img_mark[i][j][2]=new_b;
									img_mark[i][j][3]=new_theta;
									img_mark[i][j][4]=l;
									MBC_Qadd(Qmark1,i,j,likelihood,birth);
								}
								else{
									free_img((void**)tmpmaskarry[i][j].ptmpmask);
									tmpmaskarry[i][j].ptmpmask = NULL;
								}
							}
							else{
								overlap = 0;
								birth = 1;
								drawSuperEllipse2(tmpmaskarry, i, j, new_a, new_b,new_theta,n,1,tmpmask_len,cosine2,sine2);
								if(Qmark2[0][0]!=0){
									for(int m=1;m<=Qmark2[0][0];m++){
										int y = Qmark2[m][0];
										int x = Qmark2[m][1];
										overlap = elips_overlap_check(i,j,y,x,tmpmaskarry);
										if(overlap){
											birth = 0;
											break;
										}
									}
								}
								if(birth){
									img_mark[i][j][0] = 1;   //give birth to this point
									img_mark[i][j][1]=new_a;
									img_mark[i][j][2]=new_b;
									img_mark[i][j][3]=new_theta;
									img_mark[i][j][4]=l;
									MBC_Qadd(Qmark2,i,j,likelihood,birth);
								}
								else{
									free_img((void**)tmpmaskarry[i][j].ptmpmask);
									tmpmaskarry[i][j].ptmpmask = NULL;
								}
							}
						
						}

					}
				}
			}
		clock_t mid_time=clock();

		if(iter!=1){
			graph_cut(Qmark1,Qmark2,tmpmaskarry);

			for(int i = 1;i<=Qmark1[0][0];i++){
				if (Qmark1[i][3] == 0){ // delete
					int y = Qmark1[i][0];
					int x = Qmark1[i][1];
					img_mark[y][x][0] = 0;	
					MBC_Qdelete(Qmark1, tmpmaskarry, i);
					i--;
				}
			}
			for(int i = 1;i<=Qmark2[0][0];i++){
				if (Qmark2[i][3] == 0){ // delete
					int y = Qmark2[i][0];
					int x = Qmark2[i][1];
					img_mark[y][x][0] = 0;	
					MBC_Qdelete(Qmark2, tmpmaskarry, i);
					i--;
				}
				else{
					int y = Qmark2[i][0];
					int x = Qmark2[i][1];
					int likelihood = Qmark2[i][2];
					img_mark[y][x][0] = 1;	
					MBC_Qadd(Qmark1,y,x,likelihood,1);
				}
			}
		}
		 clock_t end_time=clock();
		 total_time += (end_time-start_time);
	}
	cout<< "Running time is: "<<total_time<<"ms"<<endl;
	total_time /= CLOCKS_PER_SEC;
	
	int k = 0;
	for (int i = 0;i<height;i++)
		for (int j = 0;j < width;j++)
		{
			if (img_mark[i][j][0] == 1)  // there is an object here
			{
				config->ellipse[k].center.x = j;
				config->ellipse[k].center.y = i;
				config->ellipse[k].a = img_mark[i][j][1];
				config->ellipse[k].b = img_mark[i][j][2];
				config->ellipse[k].theta = img_mark[i][j][3];
				config->ellipse[k].single_E = img_mark[i][j][4];
				config->ellipse[k].obj_id = k;
				config->ellipse[k].n = n; // ELLIPSE_POWER
				k++;
			}
		}
	config->num_obj = k;

	    //free the memory   
    for(int i=0; i<height; i++) 
    {
        for(int j=0; j<width; j++) 
        {   
            free(img_mark[i][j]);   
        }   
    }       
    for(int i=0; i<height; i++)   
    {       
        free(img_mark[i]);   
    }   
    free(img_mark);

//	free_img((void**)mask);
//	for(int i=1; i<= Qmark1[0][0]; i++) 
//	{
//		int y = (int)Qmark1[i][0]; 
//		int x = (int)Qmark1[i][1]; 
//		free_img((void**)tmpmaskarry[y][x].ptmpmask);   
//	}
    for(int i=0; i<height; i++) {
        for(int j=0; j<width; j++) {   
            if(tmpmaskarry[i][j].ptmpmask!=NULL) 
				free_img((void**)tmpmaskarry[i][j].ptmpmask);   
        }   
    }       
	free_img((void**)tmpmaskarry );
	free_img((void**)Qmark1 );
	free_img((void**)Qmark2 );
	free_img((void**)tmpmask );
	free_TIFF(&inter_img);
	return (double)total_time;
}

void err_fn(char *s)
{
	printf(s);
}

void graph_cut3(ElipsConfig *config, LinkedList newlink, int newlink_num)
{
///////////////////////////////////////////////////
//
//		        SOURCE
//		       /       \
//		     1/         \2
//		     /      3    \
//		   node0 -----> node1
//		     |   <-----   |
//		     |      4     |
//		     \            /
//		     5\          /6
//		       \        /
//		          SINK
//
///////////////////////////////////////////////////
//	typedef Graph<int,int,int> GraphType;
//	GraphType *g = new GraphType(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1); 
//	g -> add_node(); 
//	g -> add_node(); 

//	g -> add_tweights( 0,   /* capacities */  1, 5 );
//	g -> add_tweights( 1,   /* capacities */  2, 6 );
//	g -> add_edge( 0, 1,    /* capacities */  3, 4 );

//	int flow = g -> maxflow();
//
//	printf("Flow = %d\n", flow);
//	printf("Minimum cut:\n");
//	if (g->what_segment(0) == GraphType::SOURCE)
//		printf("node0 is in the SOURCE set\n");
//	else
//		printf("node0 is in the SINK set\n");
//	if (g->what_segment(1) == GraphType::SOURCE)
//		printf("node1 is in the SOURCE set\n");
//	else
//		printf("node1 is in the SINK set\n");
//	delete g;

	Node *p1, *p2;
	typedef Graph<int,int,int> GraphType;
	int num_nodes = config->mp_num+newlink_num;
	int num_edges = 2*newlink_num;
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges, &err_fn); 

	for(int i = 0;i<num_nodes;i++)
		g -> add_node(); 
	p1 = config->elink;
	for(int i = 0;i<config->mp_num;i++){
		p1 = p1->next;
		g -> add_tweights( i, p1->index->gcut_weight,  1000-p1->index->gcut_weight);
	}
	p2 = newlink;
	for(int i = 0;i<newlink_num;i++){
		p2 = p2->next;
		g -> add_tweights( i+config->mp_num, 1000-p2->index->gcut_weight, p2->index->gcut_weight );
	}
	p1 = config->elink;
	for(int i = 0;i<config->mp_num;i++){
		p1 = p1->next;
		p2 = newlink;
		for(int j = 0;j<newlink_num;j++){
			p2 = p2->next;
			int overlap = elips_overlap_check3(p1->index, p2->index);
			if(overlap)
				g -> add_edge( i, j+config->mp_num, 0, INF );
		}
	}
//	printf("maxflow alorithm..\n");
	int flow = g -> maxflow();
	p1 = config->elink;
	for(int i = 0;i<config->mp_num;i++){
		p1 = p1->next;
		if (g->what_segment(i) == GraphType::SOURCE)
			p1->index->my_type = 0;
		else
			p1->index->my_type = 1;
	}
	p2 = newlink;
	for(int i = 0;i<newlink_num;i++){
		p2 = p2->next;
		if (g->what_segment(i+config->mp_num) == GraphType::SOURCE)
			p2->index->my_type = 1;
		else
			p2->index->my_type = 0;
	}
	delete g;
}

void graph_cut33(ElipsoidConfig *config, LinkedList3D newlink, int newlink_num)
{
///////////////////////////////////////////////////
//
//		        SOURCE
//		       /       \
//		     1/         \2
//		     /      3    \
//		   node0 -----> node1
//		     |   <-----   |
//		     |      4     |
//		     \            /
//		     5\          /6
//		       \        /
//		          SINK
//
///////////////////////////////////////////////////
//	typedef Graph<int,int,int> GraphType;
//	GraphType *g = new GraphType(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1); 
//	g -> add_node(); 
//	g -> add_node(); 

//	g -> add_tweights( 0,   /* capacities */  1, 5 );
//	g -> add_tweights( 1,   /* capacities */  2, 6 );
//	g -> add_edge( 0, 1,    /* capacities */  3, 4 );

//	int flow = g -> maxflow();
//
//	printf("Flow = %d\n", flow);
//	printf("Minimum cut:\n");
//	if (g->what_segment(0) == GraphType::SOURCE)
//		printf("node0 is in the SOURCE set\n");
//	else
//		printf("node0 is in the SINK set\n");
//	if (g->what_segment(1) == GraphType::SOURCE)
//		printf("node1 is in the SOURCE set\n");
//	else
//		printf("node1 is in the SINK set\n");
//	delete g;

	Node3D *p1, *p2;
	typedef Graph<int,int,int> GraphType;
	int num_nodes = config->mp_num+newlink_num;
	int num_edges = 2*newlink_num;
	int overlap;
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges, &err_fn); 

	for(int i = 0;i<num_nodes;i++)
		g -> add_node(); 
	p1 = config->elink;
	for(int i = 0;i<config->mp_num;i++){
		p1 = p1->next;
		g -> add_tweights( i, p1->index->gcut_weight,  1000-p1->index->gcut_weight);
	}
	p2 = newlink;
	for(int i = 0;i<newlink_num;i++){
		p2 = p2->next;
		g -> add_tweights( i+config->mp_num, 1000-p2->index->gcut_weight, p2->index->gcut_weight );
	}
	p1 = config->elink;
	for(int i = 0;i<config->mp_num;i++){
		p1 = p1->next;
		p2 = newlink;
		for(int j = 0;j<newlink_num;j++){
			p2 = p2->next;
			if(p1->index->type != OBJ3D_DISK)
				overlap = elipsoid_overlap_check3(p1->index, p2->index);
			else
				overlap = elipsoid_overlap_check_disk(p1->index, p2->index);
			if(overlap)
				g -> add_edge( i, j+config->mp_num, 0, INF );
		}
	}
//	printf("maxflow alorithm..\n");
	int flow = g -> maxflow();
	p1 = config->elink;
	for(int i = 0;i<config->mp_num;i++){
		p1 = p1->next;
		if (g->what_segment(i) == GraphType::SOURCE)
			p1->index->my_type = 0;
		else
			p1->index->my_type = 1;
	}
	p2 = newlink;
	for(int i = 0;i<newlink_num;i++){
		p2 = p2->next;
		if (g->what_segment(i+config->mp_num) == GraphType::SOURCE)
			p2->index->my_type = 1;
		else
			p2->index->my_type = 0;
	}
	delete g;
}

double ellipse_likelyhood_sn(unsigned char **img,int width, int height,int y,int x,int a,int b,
						  double theta,double n,double *d,int grad_dir, double thredh,ElipsConfig *config)
{
	double px, py;
	double nx, ny, nrx, nry;
	double gx, gy, dtmp;
	double sum, sum1, sum2, sum3, sum4;
	int	   num;
	double cos_t = cos(theta);
	double sin_t = sin(theta);
	double dataterm;
	int ix, iy;
	double mean;
	double vari;

	num = 0;
	sum = 0.;
	sum1 = 0.;
	sum2 = 0.;
	sum3 = 0.;
	sum4 = 0.;
	for (int k = 0;k<SN_NUM;k++)
	{
		// gradient
		px = (double)a*config->cos2_n[k];
		py = (double)b*config->sin2_n[k];
		ix = (int)floor(px*cos_t - py*sin_t + x+0.5);
		iy = (int)floor(px*sin_t + py*cos_t + y+0.5);
		if((ix>0)&&(ix<=width-2)&&
			(iy>0)&&(iy<=height-2)){
			gx = img[iy-1][ix+1] - img[iy-1][ix-1];
			gx += img[iy][ix+1] - img[iy][ix-1];
			gx += img[iy+1][ix+1] - img[iy+1][ix-1];
			gy = img[iy+1][ix-1] - img[iy-1][ix-1];
			gy += img[iy+1][ix] - img[iy-1][ix];
			gy += img[iy+1][ix+1] - img[iy-1][ix+1];
			dtmp = sqrt(gx*gx+gy*gy+GRAD_TH);
			gx /= dtmp;
			gy /= dtmp;
			// normal
			nx = (double)b*config->nx[k];
			ny = (double)a*config->ny[k];
			dtmp = hypot(nx,ny);
			nx /= dtmp;
			ny /= dtmp;
			nrx = nx*cos_t - ny*sin_t;
			nry = nx*sin_t + ny*cos_t;
	//		normal[k].x = nrx + curve[k].x;
	//		normal[k].y = nry + curve[k].y;
			// gradient.normal
#if 0
			if((k>=SN_NUM_S1)&&(k<SN_NUM_S2)){
				sum2 += nrx*gx+nry*gy;
				num2++;
			}
			else if((k>=SN_NUM_S2)&&(k<SN_NUM_S3)){
				sum3 += nrx*gx+nry*gy;
				num3++;
			}
			else if((k>=SN_NUM_S3)&&(k<SN_NUM_S4)){
				sum4 += nrx*gx+nry*gy;
				num4++;
			}
			else{
				sum1 += nrx*gx+nry*gy;
				num1++;
			}
#endif
			dtmp = nrx*gx+nry*gy; 
			sum += dtmp;
			sum2 += dtmp*dtmp;
			num++;
		}
		else
			sum += -0.5;
	}
	if(num>SN_NUM*0.4){
#if 1
		mean = sum/(double)SN_NUM;
		vari = sum2/(double)SN_NUM-mean*mean;
		dataterm = (double)grad_dir*mean/(vari+1);
		(*d) = dataterm;
#else
		dataterm = (double)grad_dir*sum/(double)SN_NUM;
#endif
		dtmp = (dataterm-thredh)/(1.+thredh);
		dataterm = min(dtmp,1);
	}
	else{
		(*d) = 1;
		dataterm = 1;
	}
	return dataterm;
}

double ellipse_likelyhood_sn_test(unsigned char **img,int width, int height,int y,int x,int a,int b,
						  double theta,double n,double *d,double thredh,ElipsConfig *config)
{
	double px, py;
	double nx, ny, nrx, nry;
	double gx, gy, dtmp;
	double sum, sum1, sum2, sum3, sum4;
	int	   num;
	double cos_t = cos(theta);
	double sin_t = sin(theta);
	double dataterm;
	int ix, iy;
	double mean;
	double vari;
//	DPoint curve[60];

	num = 0;
	sum = 0.;
	sum1 = 0.;
	sum2 = 0.;
	sum3 = 0.;
	sum4 = 0.;
	for (int k = 0;k<SN_NUM;k++)
	{
		// gradient
		px = (double)a*config->cos2_n[k];
		py = (double)b*config->sin2_n[k];
		ix = (int)floor(px*cos_t - py*sin_t + x+0.5);
		iy = (int)floor(px*sin_t + py*cos_t + y+0.5);
		if((ix>0)&&(ix<=width-2)&&
			(iy>0)&&(iy<=height-2)){
			gx = img[iy-1][ix+1] - img[iy-1][ix-1];
			gx += img[iy][ix+1] - img[iy][ix-1];
			gx += img[iy+1][ix+1] - img[iy+1][ix-1];
			gy = img[iy+1][ix-1] - img[iy-1][ix-1];
			gy += img[iy+1][ix] - img[iy-1][ix];
			gy += img[iy+1][ix+1] - img[iy-1][ix+1];
			dtmp = sqrt(gx*gx+gy*gy+GRAD_TH);
			gx /= dtmp;
			gy /= dtmp;
			// normal
			nx = (double)b*config->nx[k];
			ny = (double)a*config->ny[k];
			dtmp = hypot(nx,ny);
			nx /= dtmp;
			ny /= dtmp;
			nrx = nx*cos_t - ny*sin_t;
			nry = nx*sin_t + ny*cos_t;
	//		normal[k].x = nrx + curve[k].x;
	//		normal[k].y = nry + curve[k].y;
			// gradient.normal
#if 0
			if((k>=SN_NUM_S1)&&(k<SN_NUM_S2)){
				sum2 += nrx*gx+nry*gy;
				num2++;
			}
			else if((k>=SN_NUM_S2)&&(k<SN_NUM_S3)){
				sum3 += nrx*gx+nry*gy;
				num3++;
			}
			else if((k>=SN_NUM_S3)&&(k<SN_NUM_S4)){
				sum4 += nrx*gx+nry*gy;
				num4++;
			}
			else{
				sum1 += nrx*gx+nry*gy;
				num1++;
			}
#endif
			dtmp = nrx*gx+nry*gy; 
			sum += dtmp;
			sum2 += dtmp*dtmp;
			num++;
		}
		else
			sum += -0.5;
	}
	if(num>SN_NUM*0.4){
#if 1
		mean = sum/(double)SN_NUM;
		vari = sum2/(double)SN_NUM-mean*mean;
		dataterm = mean/(vari+1);
#else
		dataterm = sum/(double)SN_NUM;
#endif
		dtmp = (dataterm-thredh)/(1.+thredh);
		dataterm = min(dtmp,1);
	}
	else{
		(*d) = 1;
		dataterm = 1;
	}
//	dtmp2 = (mean-thredh)/(1.+thredh);
//	(*d) = min(dtmp2,1);
	dataterm = mean;
	(*d) = vari;
	return dataterm;
}

double ellipsoid_likelyhood_sn(unsigned char ***img, int width, int height,int depth, int y,int x, int z, 
							   int a, int b, int c,	double alpha, double beta, double gamma, double n,
							   double *d, int grad_dir, double thredh, ElipsoidConfig *config)
{
	double pi = 3.1415926;
	double dtmpx, dtmpy, dtmpz;
	double nx, ny, nz;
	double gx, gy, gz, dtmp,sum,sum2;
	double dataterm;
	int xx,yy,zz,num;
	CvMat *V = cvCreateMat(3,1,CV_64FC1);
	CvMat *NV = cvCreateMat(3,1,CV_64FC1);
	CvMat *RV = cvCreateMat(3,1,CV_64FC1);
	CvMat *RVT = cvCreateMat(3,1,CV_64FC1);
	CvMat *T = cvCreateMat(3,1,CV_64FC1);
	CvMat *Rz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ry = cvCreateMat(3,3,CV_64FC1);
	CvMat *Rx = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryz = cvCreateMat(3,3,CV_64FC1);
	CvMat *R = cvCreateMat(3,3,CV_64FC1);

	cvmSet(T,0,0,(double)x);
	cvmSet(T,1,0,(double)y);
	cvmSet(T,2,0,(double)z);

	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,-sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,-sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,-sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rz,Ryz);
	cvMatMul(Rx,Ryz,R); // R = RxRyRz

	num = 0;
	sum = 0.;
	sum2 = 0.;
	for(int i = 0;i<POINT_NUM2+1;i++)
	{
		dtmpz = (double)c*config->sine_ufrm[i];
		cvmSet(V,2,0,dtmpz);
		for(int j = 0;j<config->num_ufrm[i];j++)
		{
			dtmpx = (double)a*config->coscos_ufrm[i][j];
			dtmpy = (double)b*config->cossin_ufrm[i][j];
			cvmSet(V,0,0,dtmpx);
			cvmSet(V,1,0,dtmpy);

			// rotation and translation
			cvMatMul(R,V,RV);
			cvScaleAdd(T, cvScalar(1), RV, RVT );

			xx = (int)(cvmGet(RVT,0,0)+0.5);
			yy = (int)(cvmGet(RVT,1,0)+0.5);
			zz = (int)(cvmGet(RVT,2,0)+0.5);
			if (xx>0 && xx<width-1&&yy>0 && yy<height-1 && zz>0 && zz<depth-1)
			{
				gx = 0; gy = 0; gz = 0;
				for(int aa=-1;aa<=1;aa++)
					for(int bb=-1;bb<=1;bb++){
						//dtmp = (3.-abs(aa)-abs(bb))*0.6;
						//gx += dtmp*(img[zz+aa][yy+bb][xx+1] - img[zz+aa][yy+bb][xx-1]);
						//gy += dtmp*(img[zz+aa][yy+1][xx+bb] - img[zz+aa][yy-1][xx+bb]);
						//gz += dtmp*(img[zz+1][yy+aa][xx+bb] - img[zz-1][yy+aa][xx+bb]);
						gx += img[zz+aa][yy+bb][xx+1] - img[zz+aa][yy+bb][xx-1];
						gy += img[zz+aa][yy+1][xx+bb] - img[zz+aa][yy-1][xx+bb];
						gz += img[zz+1][yy+aa][xx+bb] - img[zz-1][yy+aa][xx+bb];
					}

				dtmp = sqrt(gx*gx+gy*gy+gz*gz+GRAD_TH_3D);
				gx /= dtmp;
				gy /= dtmp;
				gz /= dtmp;
				nx = (double)b*(double)c*config->nx3d[i][j];
				ny = (double)a*(double)c*config->ny3d[i][j];
				nz = (double)a*(double)b*config->nz3d[i][j];
				dtmp = sqrt(nx*nx+ny*ny+nz*nz);
				nx = nx/dtmp;
				ny = ny/dtmp;
				nz = nz/dtmp;
				cvmSet(NV,0,0,nx);
				cvmSet(NV,1,0,ny);
				cvmSet(NV,2,0,nz);

				// rotation
				cvMatMul(R,NV,RV);

				nx = cvmGet(RV,0,0);
				ny = cvmGet(RV,1,0);
				nz = cvmGet(RV,2,0);

				// gradient.normal
				dtmp = nx*gx+ny*gy+nz*gz;
				sum += dtmp;
				sum2 += dtmp*dtmp;
				num++;
			}
			else
				sum += -0.5;
		}
	}

	if(num>config->total_num_ufrm*0.4){
#if 1
		double mean = sum/(double)config->total_num_ufrm;
		double vari = sum2/(double)config->total_num_ufrm-mean*mean;
		dataterm = (double)grad_dir*mean/(vari+1);
#else
		dataterm = (double)grad_dir*sum/(double)config->total_num_ufrm;
#endif
		(*d) = dataterm;
		dtmp = (dataterm-thredh)/(1.+thredh);
		dataterm = min(dtmp,1);
	}
	else{
		(*d) = 1;
		dataterm = 1;
	}

	cvReleaseMat(&V);
	cvReleaseMat(&NV);
	cvReleaseMat(&RV);
	cvReleaseMat(&RVT);
	cvReleaseMat(&T);
	cvReleaseMat(&Rz);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rx);
	cvReleaseMat(&Ryz);
	cvReleaseMat(&R);
	return dataterm;
}

double ellipsoid_likelyhood_cylinder(unsigned char ***img, int width, int height,int depth, int y,int x, int z, 
							   int a, int b, int c,	double alpha, double beta, double gamma, double n,
							   double *d, int grad_dir, double thredh, ElipsoidConfig *config)
{
	double pi = 3.1415926;
	double dtmpx, dtmpy, dtmpz;
	double nx, ny, nz;
	double gx, gy, gz, dtmp,sum,sum2;
	double dataterm;
	int xx,yy,zz,num;
	CvMat *V = cvCreateMat(3,1,CV_64FC1);
	CvMat *NV = cvCreateMat(3,1,CV_64FC1);
	CvMat *RV = cvCreateMat(3,1,CV_64FC1);
	CvMat *RVT = cvCreateMat(3,1,CV_64FC1);
	CvMat *T = cvCreateMat(3,1,CV_64FC1);
	CvMat *Rz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ry = cvCreateMat(3,3,CV_64FC1);
	CvMat *Rx = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryz = cvCreateMat(3,3,CV_64FC1);
	CvMat *R = cvCreateMat(3,3,CV_64FC1);

	cvmSet(T,0,0,(double)x);
	cvmSet(T,1,0,(double)y);
	cvmSet(T,2,0,(double)z);

	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,-sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,-sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,-sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rz,Ryz);
	cvMatMul(Rx,Ryz,R); // R = RxRyRz

	num = 0;
	sum = 0.;
	sum2 = 0.;
	for(int i = -c;i<=c;i++)
	{
		cvmSet(V,2,0,(double)i);
		for(int j = 0;j<config->num_ufrm[POINT_NUM4];j++)
		{
			dtmpx = (double)a*config->coscos_ufrm[POINT_NUM4][j]; // u = 0, v = 0~2pi
			dtmpy = (double)b*config->cossin_ufrm[POINT_NUM4][j]; // u = 0, v = 0~2pi
			cvmSet(V,0,0,dtmpx);
			cvmSet(V,1,0,dtmpy);

			// rotation and translation
			cvMatMul(R,V,RV);
			cvScaleAdd(T, cvScalar(1), RV, RVT );

			xx = (int)(cvmGet(RVT,0,0)+0.5);
			yy = (int)(cvmGet(RVT,1,0)+0.5);
			zz = (int)(cvmGet(RVT,2,0)+0.5);
			if (xx>0 && xx<width-1&&yy>0 && yy<height-1 && zz>0 && zz<depth-1)
			{
				gx = 0; gy = 0; gz = 0;
				for(int aa=-1;aa<=1;aa++)
					for(int bb=-1;bb<=1;bb++){
						//dtmp = (3.-abs(aa)-abs(bb))*0.6;
						//gx += dtmp*(img[zz+aa][yy+bb][xx+1] - img[zz+aa][yy+bb][xx-1]);
						//gy += dtmp*(img[zz+aa][yy+1][xx+bb] - img[zz+aa][yy-1][xx+bb]);
						//gz += dtmp*(img[zz+1][yy+aa][xx+bb] - img[zz-1][yy+aa][xx+bb]);
						gx += img[zz+aa][yy+bb][xx+1] - img[zz+aa][yy+bb][xx-1];
						gy += img[zz+aa][yy+1][xx+bb] - img[zz+aa][yy-1][xx+bb];
						gz += img[zz+1][yy+aa][xx+bb] - img[zz-1][yy+aa][xx+bb];
					}

				dtmp = sqrt(gx*gx+gy*gy+gz*gz+GRAD_TH_3D);
				gx /= dtmp;
				gy /= dtmp;
				gz /= dtmp;
				nx = (double)b*config->nx3d[POINT_NUM4][j];// u = 0, v = 0~2pi
				ny = (double)a*config->ny3d[POINT_NUM4][j];// u = 0, v = 0~2pi
				nz = 0;// u = 0, v = 0~2pi
//				nx = (double)b*config->nx[j];
//				ny = (double)a*config->ny[j];
//				nz = 0;
				dtmp = sqrt(nx*nx+ny*ny+nz*nz);
				nx = nx/dtmp;
				ny = ny/dtmp;
				nz = nz/dtmp;
				cvmSet(NV,0,0,nx);
				cvmSet(NV,1,0,ny);
				cvmSet(NV,2,0,nz);

				// rotation
				cvMatMul(R,NV,RV);

				nx = cvmGet(RV,0,0);
				ny = cvmGet(RV,1,0);
				nz = cvmGet(RV,2,0);

				// gradient.normal
				dtmp = nx*gx+ny*gy+nz*gz;
				sum += dtmp;
				sum2 += dtmp*dtmp;
				num++;
			}
			else
				sum += -0.5;
		}
	}

	int total_num = config->num_ufrm[POINT_NUM4]*(2*c+1);
	if(num>total_num*0.4){
#if 1
		double mean = sum/(double)total_num;
		double vari = sum2/(double)total_num-mean*mean;
		dataterm = (double)grad_dir*mean/(vari+1);
#else
		dataterm = (double)grad_dir*sum/(double)total_num;
#endif
		(*d) = dataterm;
		dtmp = (dataterm-thredh)/(1.+thredh);
		dataterm = min(dtmp,1);
	}
	else{
		(*d) = 1;
		dataterm = 1;
	}

	cvReleaseMat(&V);
	cvReleaseMat(&NV);
	cvReleaseMat(&RV);
	cvReleaseMat(&RVT);
	cvReleaseMat(&T);
	cvReleaseMat(&Rz);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rx);
	cvReleaseMat(&Ryz);
	cvReleaseMat(&R);
	return dataterm;
}

double ellipsoid_likelyhood_disk(unsigned char ***img, int width, int height,int depth, int y,int x, int z, 
							   int a, int b, int c,	double alpha, double beta, double gamma, double n,
							   double *d, int grad_dir, double thredh, ElipsoidConfig *config)
{
	double pi = 3.1415926;
	double dtmpx, dtmpy, dtmpz;
	double nx, ny, nz;
	double gx, gy, gz, dtmp,sum,sum2;
	double dataterm;
	int xx,yy,zz,num;
	CvMat *V = cvCreateMat(3,1,CV_64FC1);
	CvMat *NV = cvCreateMat(3,1,CV_64FC1);
	CvMat *RV = cvCreateMat(3,1,CV_64FC1);
	CvMat *RVT = cvCreateMat(3,1,CV_64FC1);
	CvMat *T = cvCreateMat(3,1,CV_64FC1);
	CvMat *Rz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ry = cvCreateMat(3,3,CV_64FC1);
	CvMat *Rx = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryz = cvCreateMat(3,3,CV_64FC1);
	CvMat *R = cvCreateMat(3,3,CV_64FC1);

	cvmSet(T,0,0,(double)x);
	cvmSet(T,1,0,(double)y);
	cvmSet(T,2,0,(double)z);

	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,-sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,-sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,-sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rz,Ryz);
	cvMatMul(Rx,Ryz,R); // R = RxRyRz

	num = 0;
	sum = 0.;
	sum2 = 0.;
	cvmSet(V,2,0,0);
	for(int j = 0;j<config->num_ufrm[POINT_NUM4];j++)
	{
		dtmpx = (double)a*config->coscos_ufrm[POINT_NUM4][j]; // u = 0, v = 0~2pi
		dtmpy = (double)b*config->cossin_ufrm[POINT_NUM4][j]; // u = 0, v = 0~2pi
		cvmSet(V,0,0,dtmpx);
		cvmSet(V,1,0,dtmpy);

		// rotation and translation
		cvMatMul(R,V,RV);
		cvScaleAdd(T, cvScalar(1), RV, RVT );

		xx = (int)(cvmGet(RVT,0,0)+0.5);
		yy = (int)(cvmGet(RVT,1,0)+0.5);
		zz = (int)(cvmGet(RVT,2,0)+0.5);
		if (xx>0 && xx<width-1&&yy>0 && yy<height-1 && zz>0 && zz<depth-1)
		{
			gx = 0; gy = 0; gz = 0;
			for(int aa=-1;aa<=1;aa++)
				for(int bb=-1;bb<=1;bb++){
					//dtmp = (3.-abs(aa)-abs(bb))*0.6;
					//gx += dtmp*(img[zz+aa][yy+bb][xx+1] - img[zz+aa][yy+bb][xx-1]);
					//gy += dtmp*(img[zz+aa][yy+1][xx+bb] - img[zz+aa][yy-1][xx+bb]);
					//gz += dtmp*(img[zz+1][yy+aa][xx+bb] - img[zz-1][yy+aa][xx+bb]);
					gx += img[zz+aa][yy+bb][xx+1] - img[zz+aa][yy+bb][xx-1];
					gy += img[zz+aa][yy+1][xx+bb] - img[zz+aa][yy-1][xx+bb];
					gz += img[zz+1][yy+aa][xx+bb] - img[zz-1][yy+aa][xx+bb];
				}

			dtmp = sqrt(gx*gx+gy*gy+gz*gz+GRAD_TH_3D);
			gx /= dtmp;
			gy /= dtmp;
			gz /= dtmp;
			nx = (double)b*(double)c*config->nx3d[POINT_NUM4][j];// u = 0, v = 0~2pi
			ny = (double)a*(double)c*config->ny3d[POINT_NUM4][j];// u = 0, v = 0~2pi
			nz = (double)a*(double)b*config->nz3d[POINT_NUM4][j];// u = 0, v = 0~2pi
			dtmp = sqrt(nx*nx+ny*ny+nz*nz);
			nx = nx/dtmp;
			ny = ny/dtmp;
			nz = nz/dtmp;
			cvmSet(NV,0,0,nx);
			cvmSet(NV,1,0,ny);
			cvmSet(NV,2,0,nz);

			// rotation
			cvMatMul(R,NV,RV);

			nx = cvmGet(RV,0,0);
			ny = cvmGet(RV,1,0);
			nz = cvmGet(RV,2,0);

			// gradient.normal
			dtmp = nx*gx+ny*gy+nz*gz;
			sum += dtmp;
			sum2 += dtmp*dtmp;
			num++;
		}
		else
			sum += -0.5;
	}

	int total_num = config->num_ufrm[POINT_NUM4];
	if(num>total_num*0.4){
#if 1
		double mean = sum/(double)total_num;
		double vari = sum2/(double)total_num-mean*mean;
		dataterm = (double)grad_dir*mean/(vari+1);
#else
		dataterm = (double)grad_dir*sum/(double)total_num;
#endif
		(*d) = dataterm;
		dtmp = (dataterm-thredh)/(1.+thredh);
		dataterm = min(dtmp,1);
	}
	else{
		(*d) = 1;
		dataterm = 1;
	}

	cvReleaseMat(&V);
	cvReleaseMat(&NV);
	cvReleaseMat(&RV);
	cvReleaseMat(&RVT);
	cvReleaseMat(&T);
	cvReleaseMat(&Rz);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rx);
	cvReleaseMat(&Ryz);
	cvReleaseMat(&R);
	return dataterm;
}

#define SYM_TH 10.
#define FIT_MARGIN 5
#define FIT_VARI 1.//0.1
#define NEAR_PNT_RANGE 6
#define Z_FACTOR 2
#define Z_FACTOR2 4 // Z_FACTOR*Z_FACTOR
#define NEAR_PNT_ZRANGE 3 // NEAR_PNT_RANGE/Z_FACTOR
#define IMG_SIZE 40
double ellipsoid_fitlikelyhood(int ***pnt, int width, int height,int depth,int y,int x, int z,
						int a, int b, int c, double alpha, double beta, double gamma, double n, int id, 
						double *d, double thredh, ElipsoidConfig *config, int draw)
{
	int nx,ny,nz;
	int mink, maxk;
	int mini, maxi;
	int minj, maxj;
	int minr;
	int tmp;
	int max_dist = 3*NEAR_PNT_RANGE*NEAR_PNT_RANGE;
	double var = (double)max_dist*0.05;
	double dtmp, dsum;
	double likely;
	CvMat *V = cvCreateMat(3,1,CV_64FC1);
	CvMat *RV = cvCreateMat(3,1,CV_64FC1);
	CvMat *TV = cvCreateMat(3,1,CV_64FC1);
	CvMat *T = cvCreateMat(3,1,CV_64FC1);
	CvMat *Rz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ry = cvCreateMat(3,3,CV_64FC1);
	CvMat *Rx = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryx = cvCreateMat(3,3,CV_64FC1);
	CvMat *R = cvCreateMat(3,3,CV_64FC1);
	CvMat *IR = cvCreateMat(3,3,CV_64FC1);
	unsigned char ***img;

	if(draw){
		img = (unsigned char ***)get_volume(IMG_SIZE, IMG_SIZE, IMG_SIZE, sizeof(unsigned char));
		for (int kk = 0; kk < IMG_SIZE; kk++ )
			for (int ii = 0; ii < IMG_SIZE; ii++ )
				for (int jj = 0; jj < IMG_SIZE; jj++ ) {
					int zz = z+kk-IMG_SIZE/2;
					int yy = y+ii-IMG_SIZE/2;
					int xx = x+jj-IMG_SIZE/2;
					if(xx>=0 && xx<width && yy>=0 && yy<height && zz>=0 && zz<depth)
						img[kk][ii][jj] = pnt[zz][yy][xx]+100;
					else
						img[kk][ii][jj] = 0;
				}
	}

	cvmSet(T,0,0,(double)x);
	cvmSet(T,1,0,(double)y);
	cvmSet(T,2,0,(double)z);

	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,-sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,-sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,-sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rz,Ryz);
	cvMatMul(Rx,Ryz,R); // R = RxRyRz

	dsum = 0;
	for (int i = 2;i<=POINT_NUM2-2;i++){
		for (int j = 0;j<config->num_ufrm[i];j++){
			cvmSet(V,0,0,a*config->coscos_ufrm[i][j]);
			cvmSet(V,1,0,b*config->cossin_ufrm[i][j]);
			cvmSet(V,2,0,c*config->sine_ufrm[i]);

			// ratation and translation
			cvMatMul(R,V,RV);
			cvScaleAdd(T, cvScalar(1), RV, V );

			nz = (int)(cvmGet(V,2,0)+0.5);
			ny = (int)(cvmGet(V,1,0)+0.5);
			nx = (int)(cvmGet(V,0,0)+0.5);
			if (nx>=0 && nx<width && ny>=0 && ny<height && nz>=0 && nz<depth)
			{
				mink = max(0,nz-NEAR_PNT_ZRANGE);
				maxk = min(depth-1,nz+NEAR_PNT_ZRANGE);
				mini = max(0,ny-NEAR_PNT_RANGE);
				maxi = min(height-1,ny+NEAR_PNT_RANGE);
				minj = max(0,nx-NEAR_PNT_RANGE);
				maxj = min(width-1,nx+NEAR_PNT_RANGE);
				minr = max_dist;
				for(int kk=mink; kk<=maxk ;kk++){
					for(int ii=mini; ii<=maxi ;ii++){
						for(int jj=minj; jj<=maxj ;jj++){
							if(pnt[kk][ii][jj]==id){
								tmp = Z_FACTOR2*(kk-nz)*(kk-nz)+(ii-ny)*(ii-ny)+(jj-nx)*(jj-nx);
								if(tmp<minr)
									minr = tmp;
							}
						}
					}
				}
				dtmp = exp(-(double)minr/var);
				dsum += dtmp;
				if(draw){
					minr = (int)(dtmp*255);
					if(minr>255) minr = 255;
					int xx = nx-x+IMG_SIZE/2;
					int yy = ny-y+IMG_SIZE/2;
					int zz = nz-z+IMG_SIZE/2;
					if(zz>=0 && zz<IMG_SIZE && yy>=0 && yy <IMG_SIZE && xx>=0 && xx<IMG_SIZE)
						img[zz][yy][xx] = minr;
				}
			}
		}
	}
	dsum /= (double)config->total_num_ufrm;

	(*d) = dsum;
	if (dsum<thredh)
		likely = 1-dsum/thredh;
	else
		likely = exp(-(dsum-thredh)/(0.5*dsum))-1;

	//FILE *fp;
	//if ((fp = fopen("data.txt", "wb")) == NULL ) {
	//	printf("Cannot open file \n");
	//	exit(1);
	//}
	//fprintf(fp,"x,y,z \r\n");
	//fprintf(fp,"%d,%d,%d \r\n", x,y,z);
	//fprintf(fp,"a,b,c \r\n");
	//fprintf(fp,"%d,%d,%d \r\n", a,b,c);
	//fprintf(fp,"alpha,beta,gamma \r\n");
	//fprintf(fp,"%f,%f,%f \r\n", alpha,beta,gamma);


//	fclose(fp);

	if(draw){
		struct TIFF_img inter_img;
		FILE *fp;
		char name_buff[2000];
		get_TIFF ( &inter_img, IMG_SIZE, IMG_SIZE, 'g' );

		for (int kk = 0; kk < IMG_SIZE; kk++ ){
			for (int ii = 0; ii < IMG_SIZE; ii++ )
				for (int jj = 0; jj < IMG_SIZE; jj++ ) {
				  inter_img.mono[ii][jj] = img[kk][ii][jj];
				}
			sprintf(name_buff, "point%03d.tif", kk);
			/* open image file */
			if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
			fprintf ( stderr, "cannot open file image.tif\n");
			exit ( 1 );
			}
			/* write image */
			if ( write_TIFF ( fp, &inter_img ) ) {
			fprintf ( stderr, "error writing TIFF file \n" );
			exit ( 1 );
			}
			/* close image file */
			fclose ( fp );
		}
		free_TIFF(&inter_img);
		free_volume((void ***)img);
	}

	cvReleaseMat(&V);
	cvReleaseMat(&RV);
	cvReleaseMat(&TV);
	cvReleaseMat(&T);
	cvReleaseMat(&Rz);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rx);
	cvReleaseMat(&Ryz);
	cvReleaseMat(&Ryx);
	cvReleaseMat(&R);
	return likely;
}

double ellipsoid_fitlikelyhood2(int ***pnt, int width, int height,int depth,int y,int x, int z,
						int a, int b, int c, double alpha, double beta, double gamma, double n, int id, 
						double *d, double thredh, ElipsoidConfig *config, int draw)
{
	int minx = width;
	int maxx = 0;
	int miny = height;
	int maxy = 0;
	int minz = depth;
	int maxz = 0;
	int NX,NY,NZ;
	double tmpx, tmpy, tmpz, dsum;
	double likely;
	CvMat *V = cvCreateMat(3,1,CV_64FC1);
	CvMat *RV = cvCreateMat(3,1,CV_64FC1);
	CvMat *TV = cvCreateMat(3,1,CV_64FC1);
	CvMat *T = cvCreateMat(3,1,CV_64FC1);
	CvMat *Rz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ry = cvCreateMat(3,3,CV_64FC1);
	CvMat *Rx = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryx = cvCreateMat(3,3,CV_64FC1);
	CvMat *R = cvCreateMat(3,3,CV_64FC1);
	CvMat *IR = cvCreateMat(3,3,CV_64FC1);
	unsigned char ***img;

	cvmSet(T,0,0,(double)x);
	cvmSet(T,1,0,(double)y);
	cvmSet(T,2,0,(double)z);

	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,-sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,-sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,-sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rz,Ryz);
	cvMatMul(Rx,Ryz,R); // R = RxRyRz

	for(int k = -c; k <= c;k+=2*c)
		for(int i = -b; i <= b;i+=2*b)
			for(int j = -a; j <= a;j+=2*a){
				cvmSet(V,0,0,(double)j);
				cvmSet(V,1,0,(double)i);
				cvmSet(V,2,0,(double)k);
				// ratation and translation
				cvMatMul(R,V,RV);
				cvScaleAdd(T, cvScalar(1), RV, V );

				NZ = (int)(cvmGet(V,2,0)+0.5);
				NY = (int)(cvmGet(V,1,0)+0.5);
				NX = (int)(cvmGet(V,0,0)+0.5);
				if(minx>NX)
					minx = NX;
				if(miny>NY)
					miny = NY;
				if(minz>NZ)
					minz = NZ;
				if(maxx<NX)
					maxx = NX;
				if(maxy<NY)
					maxy = NY;
				if(maxz<NZ)
					maxz = NZ;
			}
	minx -= FIT_MARGIN;
	miny -= FIT_MARGIN;
	minz -= FIT_MARGIN;
	maxx += FIT_MARGIN;
	maxy += FIT_MARGIN;
	maxz += FIT_MARGIN;
	if(minx<0)
		minx = 0;
	if(miny<0)
		miny = 0;
	if(minz<0)
		minz = 0;
	if(maxx>width-1)
		maxx = width-1;
	if(maxy>height-1)
		maxy = height-1;
	if(maxz>depth-1)
		maxz = depth-1;
	
	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,-sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,-sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,-sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rx,Ryx);
	cvMatMul(Rz,Ryx,IR); // R- = Rz-Ry-Rx-

	//FILE *fp;
	//if ((fp = fopen("data.txt", "wb")) == NULL ) {
	//	printf("Cannot open file \n");
	//	exit(1);
	//}
	//fprintf(fp,"x,y,z \r\n");
	//fprintf(fp,"%d,%d,%d \r\n", x,y,z);
	//fprintf(fp,"a,b,c \r\n");
	//fprintf(fp,"%d,%d,%d \r\n", a,b,c);
	//fprintf(fp,"alpha,beta,gamma \r\n");
	//fprintf(fp,"%f,%f,%f \r\n", alpha,beta,gamma);
	if(draw){
		img = (unsigned char ***)get_volume(maxx-minx, maxy-miny, maxz-minz, sizeof(unsigned char));
		for (int kk = 0; kk < maxz-minz; kk++ )
			for (int ii = 0; ii < maxy-miny; ii++ )
				for (int jj = 0; jj < maxx-minx; jj++ ) {
				  img[kk][ii][jj] = 0;
				}
	}

	dsum = 0.;
	for (int k = minz;k<=maxz;k++)
		for (int i = miny;i<=maxy;i++)
			for (int j =minx;j<=maxx;j++)
			{
				if(pnt[k][i][j]==id){
//	fprintf(fp,"%d,%d,%d \r\n", k,i,j);
					cvmSet(V,0,0,(double)j);
					cvmSet(V,1,0,(double)i);
					cvmSet(V,2,0,(double)k);
					// ratation and translation
					cvScaleAdd(T, cvScalar(-1), V, TV );
					cvMatMul(IR,TV,V);
					tmpx = cvmGet(V,0,0);
					tmpy = cvmGet(V,1,0);
					tmpz = cvmGet(V,2,0);
					tmpx = pow(fabs(tmpx/(double)a),n);
					tmpy = pow(fabs(tmpy/(double)b),n);
					tmpz = pow(fabs(tmpz/(double)c),n);
					tmpx = tmpx+tmpy+tmpz;
//					tmpx = pow(tmpx,0.333333333);
					if((tmpx<2.)&&(tmpx>=0)){
#if 1
						tmpx = tmpx-1;
						tmpx = exp(-tmpx*tmpx/FIT_VARI);
						if(draw) img[k-minz][i-miny][j-minx]=tmpx*200;
						dsum += tmpx-0.5;
#else
						if(tmpx>1.)
							tmpx = 2.-tmpx;
						//tmpx = exp(-tmpx/FIT_VARI);
//						img[k-minz][i-miny][j-minx]=tmpx*200;
						tmpx = tmpx-0.75;
						dsum += tmpx;
						//if(tmpx>0.1){
						//	dsum += tmpx;
						//}
#endif
					}
				}
			}
	dsum = dsum/((double)((a+b)*c));
	(*d) = dsum;
	if (dsum<thredh)
		likely = 1-dsum/thredh;
	else
		likely = exp(-(dsum-thredh)/(0.5*dsum))-1;
//	fclose(fp);

	if(draw){
		for (int i = 0;i<POINT_NUM;i++)
			for (int j = 0;j<POINT_NUM;j++)
			{
				cvmSet(V,0,0,a*config->coscos[i][j]);
				cvmSet(V,1,0,b*config->cossin[i][j]);
				cvmSet(V,2,0,c*config->sine[i]);

				// ratation and translation
				cvMatMul(R,V,RV);
				cvScaleAdd(T, cvScalar(1), RV, V );

				int NZ = (int)(cvmGet(V,2,0)+0.5);
				int NY = (int)(cvmGet(V,1,0)+0.5);
				int NX = (int)(cvmGet(V,0,0)+0.5);
				if (NX>=0 && NX<width&&NY>=0 && NY<height && NZ>=0 && NZ<depth)
				{
					img[NZ-minz][NY-miny][NX-minx] = 255;
				}
			}

		struct TIFF_img inter_img;
		FILE *fp;
		char name_buff[2000];
		get_TIFF ( &inter_img, maxy-miny, maxx-minx, 'g' );

		for (int kk = 0; kk < maxz-minz; kk++ ){
			for (int ii = 0; ii < maxy-miny; ii++ )
				for (int jj = 0; jj < maxx-minx; jj++ ) {
				  inter_img.mono[ii][jj] = img[kk][ii][jj];
				}
			sprintf(name_buff, "point%03d.tif", kk);
			/* open image file */
			if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
			fprintf ( stderr, "cannot open file image.tif\n");
			exit ( 1 );
			}
			/* write image */
			if ( write_TIFF ( fp, &inter_img ) ) {
			fprintf ( stderr, "error writing TIFF file \n" );
			exit ( 1 );
			}
			/* close image file */
			fclose ( fp );
		}
		free_TIFF(&inter_img);
		free_volume((void ***)img);
	}

	cvReleaseMat(&V);
	cvReleaseMat(&RV);
	cvReleaseMat(&TV);
	cvReleaseMat(&T);
	cvReleaseMat(&Rz);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rx);
	cvReleaseMat(&Ryz);
	cvReleaseMat(&Ryx);
	cvReleaseMat(&R);
	return likely;
}

double ellipsoid_fitlikelyhoodold(unsigned char ***pnt, int width, int height,int depth,int y,int x, int z,
						int a, int b, int c, double alpha, double beta, double gamma, double n, double thredh, ElipsoidConfig *config)
{
	int minx = width;
	int maxx = 0;
	int miny = height;
	int maxy = 0;
	int minz = depth;
	int maxz = 0;
	int NX,NY,NZ;
	double tmpx, tmpy, tmpz, dsum;
	double likely;
	CvMat *V = cvCreateMat(3,1,CV_64FC1);
	CvMat *RV = cvCreateMat(3,1,CV_64FC1);
	CvMat *TV = cvCreateMat(3,1,CV_64FC1);
	CvMat *T = cvCreateMat(3,1,CV_64FC1);
	CvMat *Rz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ry = cvCreateMat(3,3,CV_64FC1);
	CvMat *Rx = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryz = cvCreateMat(3,3,CV_64FC1);
	CvMat *Ryx = cvCreateMat(3,3,CV_64FC1);
	CvMat *R = cvCreateMat(3,3,CV_64FC1);
	CvMat *IR = cvCreateMat(3,3,CV_64FC1);
	unsigned char ***img;

	cvmSet(T,0,0,(double)x);
	cvmSet(T,1,0,(double)y);
	cvmSet(T,2,0,(double)z);

	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,-sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,-sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,-sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rz,Ryz);
	cvMatMul(Rx,Ryz,R); // R = RxRyRz

	for(int k = -c; k <= c;k+=2*c)
		for(int i = -b; i <= b;i+=2*b)
			for(int j = -a; j <= a;j+=2*a){
				cvmSet(V,0,0,(double)j);
				cvmSet(V,1,0,(double)i);
				cvmSet(V,2,0,(double)k);
				// ratation and translation
				cvMatMul(R,V,RV);
				cvScaleAdd(T, cvScalar(1), RV, V );

				NZ = (int)(cvmGet(V,2,0)+0.5);
				NY = (int)(cvmGet(V,1,0)+0.5);
				NX = (int)(cvmGet(V,0,0)+0.5);
				if(minx>NX)
					minx = NX;
				if(miny>NY)
					miny = NY;
				if(minz>NZ)
					minz = NZ;
				if(maxx<NX)
					maxx = NX;
				if(maxy<NY)
					maxy = NY;
				if(maxz<NZ)
					maxz = NZ;
			}
	minx -= FIT_MARGIN;
	miny -= FIT_MARGIN;
	minz -= FIT_MARGIN;
	maxx += FIT_MARGIN;
	maxy += FIT_MARGIN;
	maxz += FIT_MARGIN;
	if(minx<0)
		minx = 0;
	if(miny<0)
		miny = 0;
	if(minz<0)
		minz = 0;
	if(maxx>width-1)
		maxx = width-1;
	if(maxy>height-1)
		maxy = height-1;
	if(maxz>depth-1)
		maxz = depth-1;
	
	cvmSet(Rx,0,0,1);
	cvmSet(Rx,0,1,0);
	cvmSet(Rx,0,2,0);
	cvmSet(Rx,1,0,0);
	cvmSet(Rx,1,1,cos(alpha));
	cvmSet(Rx,1,2,sin(alpha));
	cvmSet(Rx,2,0,0);
	cvmSet(Rx,2,1,-sin(alpha));
	cvmSet(Rx,2,2,cos(alpha));

	cvmSet(Ry,0,0,cos(beta));
	cvmSet(Ry,0,1,0);
	cvmSet(Ry,0,2,-sin(beta));
	cvmSet(Ry,1,0,0);
	cvmSet(Ry,1,1,1);
	cvmSet(Ry,1,2,0);
	cvmSet(Ry,2,0,sin(beta));
	cvmSet(Ry,2,1,0);
	cvmSet(Ry,2,2,cos(beta));

	cvmSet(Rz,0,0,cos(gamma));
	cvmSet(Rz,0,1,sin(gamma));
	cvmSet(Rz,0,2,0);
	cvmSet(Rz,1,0,-sin(gamma));
	cvmSet(Rz,1,1,cos(gamma));
	cvmSet(Rz,1,2,0);
	cvmSet(Rz,2,0,0);
	cvmSet(Rz,2,1,0);
	cvmSet(Rz,2,2,1);

	cvMatMul(Ry,Rx,Ryx);
	cvMatMul(Rz,Ryx,IR); // R- = Rz-Ry-Rx-

	//FILE *fp;
	//if ((fp = fopen("data.txt", "wb")) == NULL ) {
	//	printf("Cannot open file \n");
	//	exit(1);
	//}
	//fprintf(fp,"x,y,z \r\n");
	//fprintf(fp,"%d,%d,%d \r\n", x,y,z);
	//fprintf(fp,"a,b,c \r\n");
	//fprintf(fp,"%d,%d,%d \r\n", a,b,c);
	//fprintf(fp,"alpha,beta,gamma \r\n");
	//fprintf(fp,"%f,%f,%f \r\n", alpha,beta,gamma);
	img = (unsigned char ***)get_volume(maxx-minx, maxy-miny, maxz-minz, sizeof(unsigned char));
	for (int kk = 0; kk < maxz-minz; kk++ )
		for (int ii = 0; ii < maxy-miny; ii++ )
			for (int jj = 0; jj < maxx-minx; jj++ ) {
			  img[kk][ii][jj] = 0;
			}

	dsum = 0.;
	for (int k = minz;k<=maxz;k++)
		for (int i = miny;i<=maxy;i++)
			for (int j =minx;j<=maxx;j++)
			{
				if(pnt[k][i][j]){
//	fprintf(fp,"%d,%d,%d \r\n", k,i,j);
					cvmSet(V,0,0,(double)j);
					cvmSet(V,1,0,(double)i);
					cvmSet(V,2,0,(double)k);
					// ratation and translation
					cvScaleAdd(T, cvScalar(-1), V, TV );
					cvMatMul(IR,TV,V);
					tmpx = cvmGet(V,0,0);
					tmpy = cvmGet(V,1,0);
					tmpz = cvmGet(V,2,0);
					tmpx = pow(fabs(tmpx/(double)a),n);
					tmpy = pow(fabs(tmpy/(double)b),n);
					tmpz = pow(fabs(tmpz/(double)c),n);
					tmpx = pow(tmpx,0.333333333);
					tmpx = tmpx+tmpy+tmpz-1;
					tmpx = exp(-tmpx*tmpx/FIT_VARI);
					img[k-minz][i-miny][j-minx]=(unsigned char)tmpx*200;
					dsum += tmpx;
					//if(tmpx>0.1){
					//	dsum += tmpx;
					//}
				}
			}
	dsum = dsum/((double)(a*b*c));
#if 1		
	for (int i = 0;i<POINT_NUM;i++)
		for (int j = 0;j<POINT_NUM;j++)
		{
			cvmSet(V,0,0,a*config->coscos[i][j]);
			cvmSet(V,1,0,b*config->cossin[i][j]);
			cvmSet(V,2,0,c*config->sine[i]);

			// ratation and translation
			cvMatMul(R,V,RV);
			cvScaleAdd(T, cvScalar(1), RV, V );

			int NZ = (int)(cvmGet(V,2,0)+0.5);
			int NY = (int)(cvmGet(V,1,0)+0.5);
			int NX = (int)(cvmGet(V,0,0)+0.5);
		//	if (NX>=0 && NX<width&&NY>=0 && NY<height && NZ>=0 && NZ<depth)
		//	{
			//	img[NZ-minz][NY-miny][NX-minx] = 255;
		//	}
		}

	struct TIFF_img inter_img;
	FILE *fp;
	char name_buff[2000];
	get_TIFF ( &inter_img, maxy-miny, maxx-minx, 'g' );

	for (int kk = 0; kk < maxz-minz; kk++ ){
		for (int ii = 0; ii < maxy-miny; ii++ )
			for (int jj = 0; jj < maxx-minx; jj++ ) {
			  inter_img.mono[ii][jj] = img[kk][ii][jj];
			}
		sprintf(name_buff, "point%03d.tif", kk);
		/* open image file */
		if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
		fprintf ( stderr, "cannot open file image.tif\n");
		exit ( 1 );
		}
		/* write image */
		if ( write_TIFF ( fp, &inter_img ) ) {
		fprintf ( stderr, "error writing TIFF file \n" );
		exit ( 1 );
		}
		/* close image file */
		fclose ( fp );
	}
	free_TIFF(&inter_img);
#endif
	if (dsum<thredh)
		likely = 1-dsum/thredh;
	else
		likely = exp(-(dsum-thredh)/(3*dsum))-1;
//	fclose(fp);

	cvReleaseMat(&V);
	cvReleaseMat(&RV);
	cvReleaseMat(&TV);
	cvReleaseMat(&T);
	cvReleaseMat(&Rz);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rx);
	cvReleaseMat(&Ryz);
	cvReleaseMat(&Ryx);
	cvReleaseMat(&R);
	free_volume((void ***)img);
	return dsum;
}

double MBC_ellipse3DMain3(unsigned char ***img, double ***Pmask, ElipsConfig *config, Config_Ellipses *cfg, int width, int height, int depth,
						 MPP_Parameters mpp, double *iter_time, IplImage *image, const char* win_name)
{
	int r_x = width-2*MARGIN;
	int r_y = height-2*MARGIN;
	int r_z = depth-1;
	int max_a  = mpp.max_a;
	int max_b = mpp.max_b;
	int min_a  = mpp.min_a;
	int min_b = mpp.min_b;
	int r_a = max_a-min_a;
	int r_b = max_a-min_b;

	double n = mpp.n;

	double pi = CV_PI;

	double max_theta1 = mpp.max_theta;
	double min_theta1 = mpp.min_theta;

	double max_theta2 = MAX_THETA2;
	double min_theta2 = MIN_THETA2;
	double r_theta = max_theta1-min_theta1;
	int kmax = height*width*depth/200;//mpp.aux; ///64

	double pen  = mpp.pen;		// overlap penalty
	double intra_pen = mpp.intra_pen;
	double inter_pen = mpp.inter_pen;
	double thredh = mpp.thredh; // Eq. (5) T
	int iter_num = mpp.iter_num;

	double b_zero = mpp.b_zero;
	double delta_mpp = mpp.delta;	// 0.9 sigma for death step 
	double beta_mpp = mpp.beta_mpp;	// 10  alpha for death step
	double F = mpp.F;			// sigma and alpha decreasing and increasing factor
	int perturb = 1;//mpp.perturbation;
	double prior = 0, inter_prior = 0;
	CvScalar green = CV_RGB(0,255,0);
	CvScalar white = CV_RGB(255,255,255);
	CvScalar orange = CV_RGB(255,128,0);
	CvScalar yellow = CV_RGB(255,255,0);
	Node *pobj;
	int id[MAX_IMAGE_NUM] = {0,};

	struct TIFF_img inter_img;
	int tag = 0;
	int select = 0;
	int select2 = 0;
	get_TIFF ( &inter_img, height, width, 'g' );
	LinkedList newlink = LinkedListInit();
	int newlink_num = 0;
	double overlap;
	int display_img_num = 26;
//	double dtmp;

	int tmpmask_len = 2*(mpp.max_a+9);
	double ***img_mark = (double ***)get_volume(width,height,depth,sizeof(double));   // mark exist
	double l_th = mpp.l_th;

	if(!mpp.inter_en) 
		inter_prior = 0;

    for(int k=0; k<depth; k++)   
		for (int i = 0;i<height;i++)
			for (int j = 0;j<width;j++){
				img_mark[k][i][j] = 0;
			}

	clock_t start_time=clock();
	for (int iter = 1;iter <= iter_num; iter++)//25
	{
		clock_t mid_time1=clock();
		if(!(iter%10)) printf("************* iter = %d/%d\n",iter,iter_num);

		// birth step
		double p1;
		if(perturb==0)
			p1 = 0; // MBC
		else
			p1 = rand()/(double)(RAND_MAX); // MBC with local perturbation

		if ((iter == 1)||(p1 < 0.3)){ // birth step
			newlink_num = 0;
			for(int k=0; (k<kmax)&&(newlink_num<N_MAX); k++) { 
				ellipseObj *elips = (ellipseObj *)malloc(sizeof(ellipseObj));
				elips->n = (int)n;
				elips->x = (int)floor((rand()/(double)(RAND_MAX)*(r_x)+MARGIN)+0.5);
				elips->y = (int)floor((rand()/(double)(RAND_MAX)*(r_y)+MARGIN)+0.5);
				elips->z = (int)floor((rand()/(double)(RAND_MAX)*(r_z))+0.5);
				if(Pmask[elips->z][elips->y][elips->x]){
					elips->a = (int)floor((rand()/(double)(RAND_MAX)*(r_a)+min_a)+0.5);
					elips->b = (int)floor((rand()/(double)(RAND_MAX)*(r_b)+min_b)+0.5);
					elips->theta = rand()/(double)(RAND_MAX)*(r_theta)+min_theta1;

					elips->dataterm = ellipse_likelyhood_sn(img[elips->z],width,height,elips->y,elips->x,
						elips->a,elips->b,elips->theta,elips->n,&(elips->dist),mpp.grad_dir,thredh,config);
					if(elips->dataterm<l_th){
						elips->gcut_weight = (int)(500.*(1+elips->dataterm));
						drawSuperEllipse3(elips, 1,tmpmask_len,config);
						if (iter == 1)
						{
							overlap = elips_overlap_check_all(elips, config->elink, config->mp_num);
							if(!overlap){
								config->mp_num++;
								LinkedListInsert(config->elink, config->mp_num, elips); 
								newlink_num++;
							}
							else{
								free_img((void**)elips->ptmpmask);
								elips->ptmpmask = NULL;
								free(elips);
							}
						}
						else
						{
							overlap = elips_overlap_check_all(elips, newlink, newlink_num);
							if(!overlap){
								newlink_num++;
								LinkedListInsert(newlink,newlink_num, elips); 
							}
							else{
								free_img((void**)elips->ptmpmask);
								elips->ptmpmask = NULL;
								free(elips);
//								printf("k=%d",k);
							}
						}
					}
					else
						free(elips);
				}
				else
					free(elips);
			}
		}
		else if((p1>=0.3)&&(p1<1.1)){	// BND in Neighborhood & perturbation
			pobj = config->elink;
			newlink_num = 0;
			int tmp;
			double dtmp = (double)(rand()/(double)(RAND_MAX));
//			if(mpp.aux==1)
				if(dtmp<0.3)
					tmp = 1;
				else if(dtmp<0.6)
					tmp = -1;
				else
					tmp = 0;
//			else
//				tmp = 0;

			for(int k=0; k<config->mp_num; k++) { 
				pobj = pobj->next;
				ellipseObj *elips = (ellipseObj *)malloc(sizeof(ellipseObj));
				elips->n = (int)n;
				elips->z = pobj->index->z + tmp;
				if ((elips->z<0)||(elips->z>depth-1))
					elips->z = pobj->index->z;
				int valid, try_num=0;
				do{
					try_num++;
					if(try_num>100)
						break;
					elips->x = (int)floor((rand()/(double)(RAND_MAX)*(2*EPSILON)-EPSILON)+0.5);
					elips->y = (int)floor((rand()/(double)(RAND_MAX)*(2*EPSILON)-EPSILON)+0.5);
					elips->x += pobj->index->x;
					elips->y += pobj->index->y;
					if((elips->x>=MARGIN)&&(elips->x<=width-MARGIN)
						&&(elips->y>=MARGIN)&&(elips->y<=height-MARGIN)){
						if(Pmask[elips->z][elips->y][elips->x]){
							elips->a = (int)floor((rand()/(double)(RAND_MAX)*(2*BND_NBR_DILATION)-BND_NBR_DILATION)+0.5);
							elips->b = (int)floor((rand()/(double)(RAND_MAX)*(2*BND_NBR_DILATION)-BND_NBR_DILATION)+0.5);
							elips->a += pobj->index->a;
							elips->b += pobj->index->b;
							if((elips->a>=mpp.min_a)&&(elips->a<=mpp.max_a)
								&&(elips->b>=mpp.min_b)&&(elips->b<=mpp.max_b)){
								elips->theta = rand()/(double)(RAND_MAX)*(2*BND_NBR_ANGLE)-BND_NBR_ANGLE;
								elips->theta += pobj->index->theta;
								if((elips->theta>=mpp.min_theta)&&(elips->theta<=mpp.max_theta))
									valid = 1;
								else
									valid = 0;
							}
							else
								valid = 0;
						}
						else
							valid = 0;
					}
					else
						valid = 0;
				}while(!valid);

				if(valid){
					elips->dataterm = ellipse_likelyhood_sn(img[elips->z],width,height,elips->y,elips->x,
							elips->a,elips->b,elips->theta,elips->n,&(elips->dist),mpp.grad_dir,thredh,config);

					if(elips->dataterm<l_th){
						elips->gcut_weight = (int)(500.*(1+elips->dataterm));

						drawSuperEllipse3(elips, 1,tmpmask_len,config);
						overlap = elips_overlap_check_all(elips,newlink,newlink_num);
						if(!overlap){
							newlink_num++;
							LinkedListInsert(newlink, newlink_num, elips); 
						}
						else{
							free_img((void**)elips->ptmpmask);
							elips->ptmpmask = NULL;
							free(elips);
						}
					}
					else
						free(elips);
				}
				else
					free(elips);
			}
		}
		else{	// local perturbation
			pobj = config->elink;
			newlink_num = 0;
			for(int k=0; k<config->mp_num; k++) { 
				pobj = pobj->next;
				ellipseObj *elips = (ellipseObj *)malloc(sizeof(ellipseObj));
				elips->n = (int)n;
				int valid;
				do{
					elips->x = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elips->y = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elips->x += pobj->index->x;
					elips->y += pobj->index->y;
					elips->z = pobj->index->z;
					if((elips->x>=MARGIN)&&(elips->x<=width-MARGIN)
						&&(elips->y>=MARGIN)&&(elips->y<=height-MARGIN)){
						if(Pmask[elips->z][elips->y][elips->x]){
							elips->a = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
							elips->b = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
							elips->a += pobj->index->a;
							elips->b += pobj->index->b;
							if((elips->a>=mpp.min_a)&&(elips->a<=mpp.max_a)
								&&(elips->b>=mpp.min_b)&&(elips->b<=mpp.max_b)){
								elips->theta = rand()/(double)(RAND_MAX)*(2*DELTA_ANGLE)-DELTA_ANGLE;
								elips->theta += pobj->index->theta;
								if((elips->theta>=mpp.min_theta)&&(elips->theta<=mpp.max_theta))
									valid = 1;
								else
									valid = 0;
							}
							else
								valid = 0;
						}
						else
							valid = 0;
					}
					else
						valid = 0;
				}while(!valid);

				elips->dataterm = ellipse_likelyhood_sn(img[elips->z],width,height,elips->y,elips->x,
						elips->a,elips->b,elips->theta,elips->n,&(elips->dist),mpp.grad_dir,thredh,config);

				if(elips->dataterm<l_th){
					elips->gcut_weight = (int)(500.*(1+elips->dataterm));

					drawSuperEllipse3(elips, 1,tmpmask_len,config);
					overlap = elips_overlap_check_all(elips,newlink,newlink_num);
					if(!overlap){
						newlink_num++;
						LinkedListInsert(newlink, newlink_num, elips); 
					}
					else{
						free_img((void**)elips->ptmpmask);
						elips->ptmpmask = NULL;
						free(elips);
					}
				}
				else
					free(elips);
			}
		}
		if(!(iter%10)) printf("obj num = %d\n",config->mp_num);

		// deadth step
		if(iter!=1){
//			printf("Graph cut start\n");
			graph_cut3(config,newlink,newlink_num);
//			printf("Graph cut end\n");

			pobj = config->elink->next;
			for(int i = 1;i<=config->mp_num;i++){
				if (pobj->index->my_type == 0){ // delete
					pobj = pobj->next;
					Qdelete_elips3(config, i);
					i--;
				}
				else
					pobj = pobj->next;
			}
			pobj = newlink->next;
			for(int i = 1;i<=newlink_num;i++){
				if (pobj->index->my_type == 0){ // delete
					pobj = pobj->next;
					Qdelete_elips3(newlink, i, &newlink_num);
					i--;
				}
				else{
					config->mp_num++;
					LinkedListInsert(config->elink,config->mp_num, pobj->index); 
					pobj = pobj->next;
					LinkedListDelete(newlink, i);//do not delete object itself
					newlink_num--;
					i--;
				}
			}
		}

		clock_t mid_time2=clock();
		if(!(iter%10))
			cout<< "Iteration running time is: "<<static_cast<double>(mid_time2-mid_time1)/CLOCKS_PER_SEC<<"s"<<endl;
		if(!(iter%config->record_step)){
			pobj = config->elink;
			int exe_time = static_cast<int>(mid_time2-start_time); //ms
			config->exe_time[(int)(iter/config->record_step)] = exe_time;
			config->num_obj[(int)(iter/config->record_step)] = config->mp_num;
			config->e_wo_inter_e[(int)(iter/config->record_step)] = 0;
			for (int k = 1;k<=config->mp_num;k++){
				pobj = pobj->next;
				config->e_wo_inter_e[(int)(iter/config->record_step)] += pobj->index->dataterm;
			}
		}
		if (mpp.show_running){
			if(kbhit()){
				getch();
				char ch = getch(); // ^ 72 v 80
				if(ch==72){
					if(display_img_num==0)
						display_img_num = depth-1;
					else
						display_img_num--;
				}
				else if(ch==80)
					display_img_num = (display_img_num+1)%depth;
				printf("%d : display_img_num = %d\n",ch,display_img_num+1);
				while(kbhit()) getch();
			}
			if((image != NULL)&&(iter>5)){
				LoadImageFromMemory(image, img[display_img_num]);
				CvFont font1;
				CvPoint pt1;
				char text[200];
				cvInitFont(&font1,CV_FONT_HERSHEY_DUPLEX,0.6,0.6,0.2,0,8);
				sprintf(text, "%d",display_img_num+1);
				pt1 = cvPoint(10, 20);
				cvPutText(image,text,pt1,&font1,white);
				pobj = config->elink;
				for (int k = 1;k<=config->mp_num;k++){
					pobj = pobj->next;
					prior = C_prior_elips3D3(config, k, pen);
					pobj->index->intralayer_e = prior;
					//ellipse.e0 = C_prior_interlayer_smooth (Qmark, tmpmaskarry,img_mark, tmpmask, tmpmask_len, width, height, depth, i, j, k, (int)a, (int)b,  theta, n, &m);
					if(pobj->index->z == display_img_num){
						//ellipse.obj_id = -1; // white
						//draw_ellipseObj( pobj->index, config, image, 0, white, TEXT_SINGLE_E_DIST, white, mpp);
						draw_ellipseObjwithsn( img[pobj->index->z],width,height,pobj->index, config, image, 0, white, TEXT_SINGLE_E_DIST, white, mpp);
					}
#if 0
					else if(pobj->index->z == display_img_num-1){
						//ellipse.obj_id = -1; // orange 
						draw_ellipseObj( pobj->index, config, image, 0, orange, TEXT_NONE, orange, mpp);
					}
					else if(pobj->index->z == display_img_num+1){
						//ellipse.obj_id = -1; // yellow
						draw_ellipseObj( pobj->index, config, image, 0, yellow, TEXT_NONE, yellow, mpp);
					}
#endif
				}
				cvShowImage(win_name, image);
				cvWaitKey((int)(1));
			}

		}
	}
	clock_t end_time=clock();
	double running_time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
	cout<< "Running time is: "<< running_time<<"s"<<endl;

	int z;
	pobj = config->elink;
	for (int k = 1;k<=config->mp_num;k++){
		pobj = pobj->next;
		z = pobj->index->z;
		prior = C_prior_elips3D3(config, k, pen);
		inter_prior = 0;
//		inter_prior = C_prior_interlayer_smooth0 (Qmark,img_mark, tmpmask, tmpmask_len, 
//			width, height, depth, i, j, k, (int)img_mark[k][i][j][1], (int)img_mark[k][i][j][2],
//			(int)img_mark[k][i][j][3], n, inter_pen);
		cfg[z].ellipse[id[z]].center.x = pobj->index->x;
		cfg[z].ellipse[id[z]].center.y = pobj->index->y;
		cfg[z].ellipse[id[z]].a = pobj->index->a;
		cfg[z].ellipse[id[z]].b = pobj->index->b;
		cfg[z].ellipse[id[z]].theta = pobj->index->theta;
		cfg[z].ellipse[id[z]].single_E = pobj->index->dataterm;
		cfg[z].ellipse[id[z]].multiple_E = prior;
		cfg[z].ellipse[id[z]].dist = pobj->index->dist;
		cfg[z].ellipse[id[z]].e0 = inter_prior;
		cfg[z].ellipse[id[z]].obj_id = id[z];
		cfg[z].ellipse[id[z]].n = pobj->index->n; // ELLIPSE_POWER
		id[z]++;
	}
	for (int k = 0;k<depth;k++){
		cfg[k].num_obj = id[k];
		printf("obj num = %d\n",cfg[k].num_obj);
	}

	char filename[100];
	sprintf(filename, "%stest0.cfg",mpp.seqprefix);
	save_config_data(filename, config, 7, 7);

	//free the memory   
	free_volume( (void***)img_mark );

	return running_time;
}

#define N_MAX3D (200)
double MBC_ellipsoidMain3(unsigned char ***img, double ***Pmask, ElipsoidConfig *config, Config_Ellipsoid *cfg,
			   int width, int height, int depth,MPP_Parameters mpp)
{
	int r_x = width-2*MARGIN;
	int r_y = height-2*MARGIN;
	int r_z = depth-2*MARGIN;
	int max_a = mpp.max_a;
	int max_b = mpp.max_b;
	int max_c = mpp.max_c;
	int min_a = mpp.min_a;
	int min_b = mpp.min_b;
	int min_c = mpp.min_c;
	int r_a = max_a-min_a;
	int r_b = max_b-min_b;
	int r_c = max_c-min_c;
	int kmax = height*width*depth/4;
	double n = mpp.n;
	double pi = CV_PI;
	double max_alpha = mpp.max_alpha;
	double min_alpha = mpp.min_alpha;
	double max_beta = mpp.max_beta;
	double min_beta = mpp.min_beta;
	double max_gamma = mpp.max_gamma;
	double min_gamma = mpp.min_gamma;
	double r_alpha = mpp.max_alpha-mpp.min_alpha;
	double r_beta  = mpp.max_beta-mpp.min_beta;
	double r_gamma = mpp.max_gamma-mpp.min_gamma;
	int b_alpha = (r_alpha>DELTA_ANGLE)? 1:0;
	int b_beta  = (r_beta >DELTA_ANGLE)? 1:0;
	int b_gamma = (r_gamma>DELTA_ANGLE)? 1:0;

	double pen  = mpp.pen; // 5 overlap penalty
	double thredh = mpp.thredh; // Eq. (5) T

	double b_zero = mpp.b_zero;
	double delta_mpp = mpp.delta; // 0.9 sigma for death step 
	double beta_mpp = mpp.beta_mpp;// 10  alpha for death step
	double F = mpp.F; // sigma and alpha decreasing and increasing factor
	int iter_num = mpp.iter_num;
	int perturb = 1;//mpp.perturbation;

	struct TIFF_img inter_img;
	int tag = 0;
	int select = 0;
	int select2 = 0;
	int k;
	double iter_runtime = 0, total_runtime = 0;
	double l_th = mpp.l_th;

	get_TIFF ( &inter_img, height, width, 'g' );

	CvScalar green = CV_RGB(0,255,0);
	CvScalar white = CV_RGB(255,255,255);
	CvScalar orange = CV_RGB(255,128,0);
	CvScalar yellow = CV_RGB(255,255,0);
	Node3D *pobj;
	int id[MAX_IMAGE_NUM] = {0,};

	LinkedList3D newlink = LinkedList3DInit();
	int newlink_num = 0;
	double overlap;
	int display_img_num = 10;

	int tmpmask_len = 2*(max_a+2+8);
	int tmpmask_depth = 2*(max_c+2+8);
	int ***tmpmask=(int ***)get_volume(tmpmask_len,tmpmask_len,tmpmask_depth,sizeof(int));

	clock_t start_time=clock();
	for (int iter = 1;iter <= iter_num; iter++)//25
	{
		clock_t mid_time1=clock();
		printf("************* iter = %d/%d\n",iter,iter_num);

		// birth step
		double p1;
		if(perturb==0)
			p1 = 0; // MBC
		else
			p1 = rand()/(double)(RAND_MAX); // MBC with local perturbation

		if ((iter == 1)||(p1 <0.5)){ // birth step
			newlink_num = 0;
			for(k=0; (k<kmax)&&(newlink_num<N_MAX3D); k++) { 
				ellipsoidObj *elipsoid = (ellipsoidObj *)malloc(sizeof(ellipsoidObj));
				elipsoid->type = mpp.type;
				elipsoid->n = (int)n;
				elipsoid->ptmpmask = NULL;
//				elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(11)+6)+0.5);
//				elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(9)+7)+0.5);
//				elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(8)+2)+0.5);
				elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(r_x)+MARGIN)+0.5);
				elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(r_y)+MARGIN)+0.5);
				elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(r_z)+MARGIN)+0.5);
				if(Pmask[elipsoid->z][elipsoid->y][elipsoid->x]){
					elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(r_a)+min_a)+0.5);
					elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(r_b)+min_b)+0.5);
					if(elipsoid->type != OBJ3D_DISK)
						elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(r_c)+min_c)+0.5);
					else
						elipsoid->c = 0;
					elipsoid->alpha = rand()/(double)(RAND_MAX)*(r_alpha)+min_alpha;
					elipsoid->beta = rand()/(double)(RAND_MAX)*(r_beta)+min_beta;
					elipsoid->gamma = rand()/(double)(RAND_MAX)*(r_gamma)+min_gamma;

					//for (int kk = 0;kk<tmpmask_depth;kk++)
					//	for (int ii = 0;ii<tmpmask_len;ii++)
					//		for (int jj = 0;jj<tmpmask_len;jj++)
					//		{
					//			tmpmask[kk][ii][jj] = 0;
					//		}					

//						elipsoid->dataterm = ellipsoid_likelyhood33(tmpmask,tmpmask_len,tmpmask_depth,img,width,height,depth,
//								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
//								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),thredh,config);
					if(elipsoid->type == OBJ3D_ELLIPSOID)
						elipsoid->dataterm = ellipsoid_likelyhood_sn(img,width,height,depth,
								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
					else
						elipsoid->dataterm = ellipsoid_likelyhood_cylinder(img,width,height,depth,
								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);

					if(elipsoid->dataterm<l_th){
						elipsoid->gcut_weight = (int)(500.*(1+elipsoid->dataterm));
						double dr = max(elipsoid->a,elipsoid->b);
						double dx = 2.*(fabs((double)elipsoid->c*sin(elipsoid->beta))+dr)+1.;
						double dy = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*sin(elipsoid->alpha))+dr)+1.;
						double dz = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*cos(elipsoid->alpha))+dr)+1.;
						if(elipsoid->type == OBJ3D_ELLIPSOID)
							drawSuperEllipsoid33(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
								elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						else if(elipsoid->type == OBJ3D_CYLINDER)
							drawSuperEllipsoid_Cylinder(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
								elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						if (iter == 1)
						{
							if(elipsoid->type != OBJ3D_DISK)
								overlap = elipsoid_overlap_check_all(elipsoid, config->elink, config->mp_num);
							else
								overlap = elipsoid_overlap_check_all_disk(elipsoid, config->elink, config->mp_num);
							if(!overlap){
								config->mp_num++;
								LinkedList3DInsert(config->elink, config->mp_num, elipsoid); 
								newlink_num++;
//			printf("k = %d, new_obj = %d\n",k,newlink_num);
							}
							else if(elipsoid->type != OBJ3D_DISK){
								free_volume((void***)elipsoid->ptmpmask);
								elipsoid->ptmpmask = NULL;
								free(elipsoid);
							}
							else {
								free(elipsoid);
							}
						}
						else
						{
							if(elipsoid->type != OBJ3D_DISK)
								overlap = elipsoid_overlap_check_all(elipsoid, newlink, newlink_num);
							else
								overlap = elipsoid_overlap_check_all_disk(elipsoid, newlink, newlink_num);
							if(!overlap){
								newlink_num++;
								LinkedList3DInsert(newlink,newlink_num, elipsoid); 
//			printf("k = %d, new_obj = %d\n",k,newlink_num);
							}
							else if(elipsoid->type != OBJ3D_DISK){
								free_volume((void***)elipsoid->ptmpmask);
								elipsoid->ptmpmask = NULL;
								free(elipsoid);
							}
							else {
								free(elipsoid);
							}
						}
					}
					else{
						free(elipsoid);
					}
				}
				else{
					free(elipsoid);
				}
			}
			printf("k = %d, new_obj = %d\n",k,newlink_num);
		}
		else {	// local perturbation
			pobj = config->elink;
			newlink_num = 0;
			for(k=0; k<config->mp_num; k++) { 
				pobj = pobj->next;
				ellipsoidObj *elipsoid = (ellipsoidObj *)malloc(sizeof(ellipsoidObj));
				elipsoid->type = pobj->index->type;
				elipsoid->n = pobj->index->n;
				elipsoid->ptmpmask = NULL;
				int valid;
				do{
					elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_TRANSLATION)-DELTAZ_TRANSLATION)+0.5);
					elipsoid->x += pobj->index->x;
					elipsoid->y += pobj->index->y;
					elipsoid->z += pobj->index->z;
					if((elipsoid->x>=MARGIN)&&(elipsoid->x<=width-MARGIN)
						&&(elipsoid->y>=MARGIN)&&(elipsoid->y<=height-MARGIN)
						&&(elipsoid->z>=MARGIN)&&(elipsoid->z<=depth-MARGIN)){
						if(Pmask[elipsoid->z][elipsoid->y][elipsoid->x]){
							elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
							elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
							elipsoid->a += pobj->index->a;
							elipsoid->b += pobj->index->b;
							if((elipsoid->a>=mpp.min_a)&&(elipsoid->a<=mpp.max_a)
								&&(elipsoid->b>=mpp.min_b)&&(elipsoid->b<=mpp.max_b)){
								if(b_alpha)
									elipsoid->alpha = rand()/(double)(RAND_MAX)*(2*DELTA_ANGLE)-DELTA_ANGLE;
								else
									elipsoid->alpha = rand()/(double)(RAND_MAX)*(2*r_alpha)-r_alpha;
								if(b_beta)
									elipsoid->beta = rand()/(double)(RAND_MAX)*(2*DELTA_ANGLE)-DELTA_ANGLE;
								else
									elipsoid->beta = rand()/(double)(RAND_MAX)*(2*r_beta)-r_beta;
								if(b_gamma)
									elipsoid->gamma = rand()/(double)(RAND_MAX)*(2*DELTA_ANGLE)-DELTA_ANGLE;
								else
									elipsoid->gamma = rand()/(double)(RAND_MAX)*(2*r_gamma)-r_gamma;
								elipsoid->alpha += pobj->index->alpha;
								elipsoid->beta  += pobj->index->beta;
								elipsoid->gamma += pobj->index->gamma;
								if((elipsoid->alpha>=mpp.min_alpha)&&(elipsoid->alpha<=mpp.max_alpha)
									&&(elipsoid->beta>=mpp.min_beta)&&(elipsoid->beta<=mpp.max_beta)
									&&(elipsoid->gamma>=mpp.min_gamma)&&(elipsoid->gamma<=mpp.max_gamma)){
										if(elipsoid->type != OBJ3D_DISK){
											elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_DILATION)-DELTAZ_DILATION)+0.5);
											elipsoid->c += pobj->index->c;
											if((elipsoid->c>=mpp.min_c)&&(elipsoid->c<=mpp.max_c))
												valid = 1;
											else
												valid = 0;
										}
										else
											valid = 1;
								}
								else
									valid = 0;
							}
							else
								valid = 0;
						}
						else
							valid = 0;
					}
					else
						valid = 0;
				}while(!valid);

				for (int kk = 0;kk<tmpmask_depth;kk++)
					for (int ii = 0;ii<tmpmask_len;ii++)
						for (int jj = 0;jj<tmpmask_len;jj++)
						{
							tmpmask[kk][ii][jj] = 0;
						}					

//						elipsoid->dataterm = ellipsoid_likelyhood33(tmpmask,tmpmask_len,tmpmask_depth,img,width,height,depth,
//								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
//								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),thredh,config);
				if(elipsoid->type == OBJ3D_ELLIPSOID)
					elipsoid->dataterm = ellipsoid_likelyhood_sn(img,width,height,depth,
							elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
							elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
				else if(elipsoid->type == OBJ3D_CYLINDER)
					elipsoid->dataterm = ellipsoid_likelyhood_cylinder(img,width,height,depth,
							elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
							elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
				if(elipsoid->dataterm<l_th){
					elipsoid->gcut_weight = (int)(500.*(1+elipsoid->dataterm));
					double dr = max(elipsoid->a,elipsoid->b);
					double dx = 2.*(fabs((double)elipsoid->c*sin(elipsoid->beta))+dr)+1.;
					double dy = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*sin(elipsoid->alpha))+dr)+1.;
					double dz = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*cos(elipsoid->alpha))+dr)+1.;
					if(elipsoid->type == OBJ3D_ELLIPSOID){
						drawSuperEllipsoid33(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
							elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						overlap = elipsoid_overlap_check_all(elipsoid,newlink,newlink_num);
					}
					else if(elipsoid->type == OBJ3D_CYLINDER){
						drawSuperEllipsoid_Cylinder(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
							elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						overlap = elipsoid_overlap_check_all(elipsoid,newlink,newlink_num);
					}
					else
						overlap = elipsoid_overlap_check_all_disk(elipsoid,newlink,newlink_num);
					if(!overlap){
						newlink_num++;
						LinkedList3DInsert(newlink, newlink_num, elipsoid); 
					}
					else if(elipsoid->type != OBJ3D_DISK){
						free_volume((void***)elipsoid->ptmpmask);
						elipsoid->ptmpmask = NULL;
						free(elipsoid);
					}
					else {
						free(elipsoid);
					}
				}
				else{
					free(elipsoid);
				}
			}
			printf("k = %d, new_obj = %d\n",k,newlink_num);
		}


		// deadth step
		if((iter!=1)&&(newlink_num)){
//			printf("Graph cut start\n");
			graph_cut33(config,newlink,newlink_num);
//			printf("Graph cut end\n");

			pobj = config->elink->next;
			for(int i = 1;i<=config->mp_num;i++){
				if (pobj->index->my_type == 0){ // delete
					pobj = pobj->next;
					Qdelete_elipsoid3(config, i);
					i--;
				}
				else
					pobj = pobj->next;
			}
			pobj = newlink->next;
			for(int i = 1;i<=newlink_num;i++){
				if (pobj->index->my_type == 0){ // delete
					pobj = pobj->next;
					Qdelete_elipsoid3(newlink, i, &newlink_num);
					i--;
				}
				else{
					config->mp_num++;
					LinkedList3DInsert(config->elink,config->mp_num, pobj->index); 
					pobj = pobj->next;
					LinkedList3DDelete(newlink, i);//do not delete object itself
					newlink_num--;
					i--;
				}
			}
		}
		printf("obj num = %d\n",config->mp_num);

		clock_t mid_time2=clock();
		cout<< "Iteration running time is: "<<static_cast<double>(mid_time2-mid_time1)/CLOCKS_PER_SEC<<"s"<<endl;
		pobj = config->elink;
		int exe_time = static_cast<int>(mid_time2-start_time); //ms
		config->exe_time[iter] = exe_time;
		config->num_obj[iter] = config->mp_num;
		config->energy[iter] = 0;
		for (k = 1;k<=config->mp_num;k++){
			pobj = pobj->next;
//			prior = C_prior_elipsoid3(config, k, pen);
//			dtmp = pobj->index->dataterm
//				+pobj->index->inter_e*mpp.intra_pen/2.;
			config->energy[iter] += pobj->index->dataterm;
		}
		printf("E = %1.5f\n",config->energy[iter]);
	}
	clock_t end_time=clock();
	double running_time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
	cout<< "Running time is: "<< running_time<<"s"<<endl;

	int z;
	pobj = config->elink;
	cfg->num_obj=0;
	cfg->total_E = 0;
	for (k = 1;k<=config->mp_num;k++){
		pobj = pobj->next;
		z = pobj->index->z;
				double dataterm = ellipsoid_likelyhood33(tmpmask,tmpmask_len,tmpmask_depth,img,width,height,depth,
					pobj->index->y,pobj->index->x,pobj->index->z,pobj->index->a,pobj->index->b,pobj->index->c,
					pobj->index->alpha,pobj->index->beta,pobj->index->gamma,pobj->index->n,&(pobj->index->dist),thredh,config);
//		prior = C_prior_elipsoid3(config, k, pen);
		cfg->ellipsoid[cfg->num_obj].center.x = pobj->index->x;
		cfg->ellipsoid[cfg->num_obj].center.y = pobj->index->y;
		cfg->ellipsoid[cfg->num_obj].center.z = pobj->index->z;
		cfg->ellipsoid[cfg->num_obj].a = pobj->index->a;
		cfg->ellipsoid[cfg->num_obj].b = pobj->index->b;
		cfg->ellipsoid[cfg->num_obj].c = pobj->index->c;
		cfg->ellipsoid[cfg->num_obj].alpha = pobj->index->alpha;
		cfg->ellipsoid[cfg->num_obj].beta = pobj->index->beta;
		cfg->ellipsoid[cfg->num_obj].gamma = pobj->index->gamma;
		cfg->ellipsoid[cfg->num_obj].single_E = pobj->index->dataterm;
		cfg->ellipsoid[cfg->num_obj].e0 = pobj->index->dist;
		cfg->ellipsoid[cfg->num_obj].multiple_E = 0;
		cfg->ellipsoid[cfg->num_obj].obj_id = cfg->num_obj;
		cfg->ellipsoid[cfg->num_obj].n = pobj->index->n; // ELLIPSE_POWER
		cfg->total_E += cfg->ellipsoid[cfg->num_obj].single_E;
		cfg->num_obj++;
	}
	printf("obj num = %d\n",cfg->num_obj);

	//free the memory   
	free_volume( (void***)tmpmask );

	return running_time;
}

typedef struct ellipse0
{
	int x;
	int y;
	int z;
	int a;
	int b;
}Ellipse0;

int compare_zposition (const void * a, const void * b)
{
	if(((Ellipse0 *)a)->z > ((Ellipse0 *)b)->z) 
		return -1;
	else if(((Ellipse0 *)a)->z == ((Ellipse0 *)b)->z) 
		return 0;
	else //if(((Ellipse0 *)a)->z < ((Ellipse0 *)b)->z) 
		return 1;
}

void drawline(int **img, int x1, int y1, int x2, int y2, int color)
{
    int dx, dy, inx, iny, e;

	dx = x2 - x1;
	dy = y2 - y1;
	inx = dx > 0 ? 1 : -1;
	iny = dy > 0 ? 1 : -1;

	dx = ABS(dx);
	dy = ABS(dy);

	if(dx >= dy) {
	dy <<= 1;
	e = dy - dx;
	dx <<= 1;
		while (x1 != x2) {
			img[y1][x1] = color;
			if(e >= 0) {
				y1 += iny;
					e-= dx;
			}
			e += dy; x1 += inx;
		}
	} 
	else {
		dx <<= 1;
		e = dx - dy;
		dy <<= 1;
		while (y1 != y2) {
			img[y1][x1] = color;
			if(e >= 0) {
				x1 += inx;
				e -= dy;
			}
			e += dx; y1 += iny;
		}
	}
	img[y1][x1] = color;
}


void MBC_GeneratePointCloud(Config_Ellipses *cfg, int ***pnt, int width, int height, int depth, ElipsConfig *config)
{
	int jj0, jj1=0, jj2, ii0, ii1=0, ii2;
	for (int k = 0;k<depth;k++)
		for (int i = 0;i<height;i++)
			for (int j =0;j<width;j++)
			{
				pnt[k][i][j] = -1;
			}
	
	for (int i = 0;i<depth;i++){
		for (int j = 0;j<cfg[i].num_obj;j++){
			for (int k = 0;k<=POINT_NUM;k++)
			{
				if(k<POINT_NUM){
					jj0 = jj1;
					ii0 = ii1;
					double x = cfg[i].ellipse[j].a*config->cosine[k];
					double y = cfg[i].ellipse[j].b*config->sine[k];
					double rx = rotatex(x,y,cfg[i].ellipse[j].theta);
					double ry = rotatey(x,y,cfg[i].ellipse[j].theta);
					jj1 = (int)floor(rx + cfg[i].ellipse[j].center.x+0.5);
					ii1 = (int)floor(ry + cfg[i].ellipse[j].center.y+0.5);
					if(k==0){
						jj2 = jj1;
						ii2 = ii1;
					}
					else
						if((jj0>=0)&&(jj0<width)&&(ii0>=0)&&(ii0<height)&&
							(jj1>=0)&&(jj1<width)&&(ii1>=0)&&(ii1<height))
							drawline(pnt[i],jj0,ii0,jj1,ii1,cfg[i].ellipse[j].obj_id);
				}
				else
					if((jj2>=0)&&(jj2<width)&&(ii2>=0)&&(ii2<height)&&
						(jj1>=0)&&(jj1<width)&&(ii1>=0)&&(ii1<height))
						drawline(pnt[i],jj2,ii2,jj1,ii1,cfg[i].ellipse[j].obj_id);
			}
		}
	}

#if 0		
	struct TIFF_img inter_img;
	FILE *fp;
	char name_buff[2000];
	get_TIFF ( &inter_img, height, width, 'g' );

	for (int kk = 0; kk < depth; kk++ ){
		for (int ii = 0; ii < height; ii++ )
		for (int jj = 0; jj < width; jj++ ) {
			if(pnt[kk][ii][jj]==-1)
				inter_img.mono[ii][jj] = 0;
			else
				inter_img.mono[ii][jj] = 255;
		}
		sprintf(name_buff, "pointcloud%03d.tif", kk);
		/* open image file */
		if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
		fprintf ( stderr, "cannot open file image.tif\n");
		exit ( 1 );
		}
		/* write image */
		if ( write_TIFF ( fp, &inter_img ) ) {
		fprintf ( stderr, "error writing TIFF file \n" );
		exit ( 1 );
		}
		/* close image file */
		fclose ( fp );
	}
	free_TIFF(&inter_img);
#endif
}

int MBC_GenerateBirhMap(Config_Ellipses *clustered_cfg, Birthmap ***bm, int num_clust, int width, int height, int depth, MPP_Parameters mpp)
{
	Ellipse0 e[MAX_IMAGE_NUM];
	int num, bm_num;
	int min, max;
	CvPoint3D64f centroid, direct;
	CvPoint3D64f center[2*LINE_FIT_NUM+1];

	for (int k = 0;k<depth;k++)
		for (int i = 0;i<height;i++)
			for (int j =0;j<width;j++)
			{
				bm[k][i][j].max_c = 255;
			}
	
	bm_num = 0;
	for (int k = 0;k<num_clust;k++){
		num = 0;
		for (int i = 0;i<depth;i++){
			for (int j = 0;j<clustered_cfg[i].num_obj;j++){
				if(clustered_cfg[i].ellipse[j].obj_id == k){
					e[num].x = (int)clustered_cfg[i].ellipse[j].center.x;
					e[num].y = (int)clustered_cfg[i].ellipse[j].center.y;
					e[num].z = i;
					e[num].a = (int)clustered_cfg[i].ellipse[j].a;
					e[num].b = (int)clustered_cfg[i].ellipse[j].b;
					num++;
				}
			}
		}
		qsort(e, num, sizeof(Ellipse0), compare_zposition);
		max = e[0].z;
		if(max==depth-1)
			max = 2*depth;
		min = e[num-1].z;
		if(min==0)
			min = -depth;
		for (int ii = -CENTER_LINE_MARGIN;ii<=CENTER_LINE_MARGIN;ii++){
			for (int jj = -CENTER_LINE_MARGIN;jj<=CENTER_LINE_MARGIN;jj++){
				if(ii*ii+jj*jj<=CENTER_LINE_MARGIN2){
					if(num>=2*LINE_FIT_NUM+1){
						for (int i = LINE_FIT_NUM;i<num-LINE_FIT_NUM;i++){
							int x = jj+e[i].x;
							int y = ii+e[i].y;
							if((x>=0)&&(x<width)&&(y>=0)&&(y<height)){
								bm[e[i].z][y][x].min_c = min(e[i].z-min, max-e[i].z);
								bm[e[i].z][y][x].max_c = bm[e[i].z][y][x].min_c+1;
								bm[e[i].z][y][x].id = k;
								bm[e[i].z][y][x].min_ab = min(e[i].a, e[i].b);
								bm[e[i].z][y][x].max_ab = max(e[i].a, e[i].b);
								// beta, gamma : angle of line which include centers of superellipses
								for(int z=-LINE_FIT_NUM;z<=LINE_FIT_NUM;z++){
									center[z+LINE_FIT_NUM].x = (double)e[i+z].x;
									center[z+LINE_FIT_NUM].y = (double)e[i+z].y;
									center[z+LINE_FIT_NUM].z = (double)e[i+z].z;
								}
								//if((k==79)&&(e[i].z==35))
								//	k=k;
								line3D_fit(center, &centroid, &direct, 2*LINE_FIT_NUM+1);
								double alpha = get_alpha(direct.y,direct.z);
								bm[e[i].z][y][x].alpha = alpha;
								double dz = sin(-alpha)*direct.y+cos(-alpha)*direct.z;
								bm[e[i].z][y][x].beta = get_beta(direct.x,dz);
								bm_num++;
							}
						}
					}
					else{
						for (int i = 0;i<num;i++){
							int x = jj+e[i].x;
							int y = ii+e[i].y;
							if((x>=0)&&(x<width)&&(y>=0)&&(y<height)){
								bm[e[i].z][y][x].min_c = min(e[i].z-min, max-e[i].z);
								bm[e[i].z][y][x].max_c = bm[e[i].z][y][x].min_c;
								bm[e[i].z][y][x].id = k;
								bm[e[i].z][y][x].min_ab = min(e[i].a, e[i].b);
								bm[e[i].z][y][x].max_ab = max(e[i].a, e[i].b);
								// beta, gamma : angle of line which include centers of superellipses
								bm[e[i].z][y][x].alpha = 0;
								bm[e[i].z][y][x].beta = 0;
								bm_num++;
							}
						}
					}
				}
			}
		}
	}
#if 0		
	struct TIFF_img inter_img;
	FILE *fp;
	char name_buff[2000];
	get_TIFF ( &inter_img, height, width, 'g' );

	for (int kk = 0; kk < depth; kk++ ){
		for (int ii = 0; ii < height; ii++ )
		for (int jj = 0; jj < width; jj++ ) {
			if(bm[kk][ii][jj].max_c!=255){
				if(bm[kk][ii][jj].max_c>35)
					ii = ii;
				inter_img.mono[ii][jj] = bm[kk][ii][jj].max_c*6;
			}
			else
				inter_img.mono[ii][jj] = 0;
		}
		sprintf(name_buff, "birthmap%03d.tif", kk);
		/* open image file */
		if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
		fprintf ( stderr, "cannot open file image.tif\n");
		exit ( 1 );
		}
		/* write image */
		if ( write_TIFF ( fp, &inter_img ) ) {
		fprintf ( stderr, "error writing TIFF file \n" );
		exit ( 1 );
		}
		/* close image file */
		fclose ( fp );
	}
	free_TIFF(&inter_img);
#endif
	return bm_num;
}

typedef struct ellipse1
{
	int x;
	int y;
	int z;
	int a;
	int b;
	double  theta;
	double	min_ab;
	double	max_ab;
	double	c_pz;
	double	alpha_mean;
	double	beta_mean;
	double	gamma_mean;
	double  fit_error;
}Ellipse1;

int compare_zposition1 (const void * a, const void * b)
{
	if(((Ellipse1 *)a)->z > ((Ellipse1 *)b)->z) 
		return -1;
	else if(((Ellipse1 *)a)->z == ((Ellipse1 *)b)->z) 
		return 0;
	else //if(((Ellipse1 *)a)->z < ((Ellipse1 *)b)->z) 
		return 1;
}

typedef struct ellipsoid1
{
	int x;
	int y;
	int z;
	int a;
	int b;
	int c;
	double  alpha;
	double  beta;
	double  gamma;
	double	min_ab;
	double	max_ab;
	double	c_pz;
	double	alpha_mean;
	double	beta_mean;
	double	gamma_mean;
	double  fit_error;
}Ellipsoid1;

#define BM_STDEV1  1.5
#define BM_2VARI   4.5 //2.*BM_STDEV1*BM_STDEV1
#define BM_STDEV2  0.4
#define BM_STDEV3  0.5 // 0.4~1.5 , 0.5(final good but slow)
#define WZ		  10
//#define BM_STDEV1  1.5
//#define BM_2VARI   4.5 //2.*BM_STDEV1*BM_STDEV1
//#define BM_STDEV2  0.4
//#define BM_STDEV3  0.5 // 0.4~1.5 , 0.5(final good but slow)
//#define WZ		  10
#define min3(a,b,c)		min(min((a),(b)),(c))
#define max3(a,b,c)		max(max((a),(b)),(c))
int MBC_GenerateBirhMap_Normal(Config_Ellipses *clustered_cfg, Birthmap2 ***bm, int num_clust, int width, int height, int depth, MPP_Parameters mpp)
{
	int num;
	int minz, maxz;
	CvPoint3D64f centroid, direct;
	CvPoint3D64f center[2*LINE_FIT_NUM+1];
	int id;

    int *N_k = (int *)malloc(num_clust*sizeof(int));
	Ellipse1 **e=(Ellipse1 **)get_img(mpp.max_c*4,num_clust,sizeof(Ellipse1));

	clock_t start_time=clock();
	for (int k = 0;k<num_clust;k++)
		N_k[k] = 0;
	for (int i = 0;i<depth;i++){
		for (int j = 0;j<clustered_cfg[i].num_obj;j++){
			id = clustered_cfg[i].ellipse[j].obj_id;
			e[id][N_k[id]].x = (int)clustered_cfg[i].ellipse[j].center.x;
			e[id][N_k[id]].y = (int)clustered_cfg[i].ellipse[j].center.y;
			e[id][N_k[id]].z = i;
			e[id][N_k[id]].a = (int)clustered_cfg[i].ellipse[j].a;
			e[id][N_k[id]].b = (int)clustered_cfg[i].ellipse[j].b;
			e[id][N_k[id]].theta = clustered_cfg[i].ellipse[j].theta;
			N_k[id]++;
		}
	}
	
	for (int k = 0;k<num_clust;k++){
		qsort(e[k], N_k[k], sizeof(Ellipse1), compare_zposition1);
		maxz = e[k][0].z;
		if(maxz==depth-1)
			maxz = 2*depth;
		minz = e[k][N_k[k]-1].z;
		if(minz==0)
			minz = -depth;
		for (int i = 0;i<N_k[k];i++){
			e[k][i].min_ab = min(e[k][i].a,e[k][i].b);
			e[k][i].max_ab = max(e[k][i].a,e[k][i].b);
			e[k][i].c_pz   = min(e[k][i].z-minz, maxz-e[k][i].z)+1;
			if ((i>=LINE_FIT_NUM)&&(i<N_k[k]-LINE_FIT_NUM)){
#if 1
				double max1, max2, max3; 
				double min1, min2, min3; 
				min1 = min(e[k][i-1].a,e[k][i-1].b);
				max1 = max(e[k][i-1].a,e[k][i-1].b);
				min2 = min(e[k][i].a,e[k][i].b);
				max2 = max(e[k][i].a,e[k][i].b);
				min3 = min(e[k][i+1].a,e[k][i+1].b);
				max3 = max(e[k][i+1].a,e[k][i+1].b);
				e[k][i].min_ab = min3(min1,min2,min3);
				e[k][i].max_ab = max3(max1,max2,max3);
#endif
			// beta, gamma : angle of line which include centers of superellipses
				for(int z=-LINE_FIT_NUM;z<=LINE_FIT_NUM;z++){
					center[z+LINE_FIT_NUM].x = (double)e[k][i+z].x;
					center[z+LINE_FIT_NUM].y = (double)e[k][i+z].y;
					center[z+LINE_FIT_NUM].z = (double)e[k][i+z].z;
				}
				//if((k==57)&&(e[k][i].z==27))
				//	k=k;
				e[k][i].fit_error = line3D_fit(center, &centroid, &direct, 2*LINE_FIT_NUM+1);
				double alpha = get_alpha(direct.y,direct.z);
				e[k][i].alpha_mean = alpha;
				double dz = sin(-alpha)*direct.y+cos(-alpha)*direct.z;
				e[k][i].beta_mean = get_beta(direct.x,dz);
				e[k][i].gamma_mean = e[k][i].theta;
			}
			else{
				e[k][i].min_ab = 0; // just uniformly draw parameters
				e[k][i].alpha_mean = 0;
				e[k][i].beta_mean  = 0;
				e[k][i].gamma_mean = 0;
				e[k][i].fit_error = 0;
			}
		}
	}

	for (int z = 0;z<depth;z++)
		for (int y = 0;y<height;y++)
			for (int x =0;x<width;x++)
			{
				int mind = INF;
				int kk,ii;
				for (int k = 0;k<num_clust;k++)
					for (int i =0;i<N_k[k];i++)
					{
						int dx, dy, dz;
						dx = x-e[k][i].x;
						dx = dx*dx;
						if(dx<mind){
							dy = y-e[k][i].y;
							dy = dy*dy;
							if(dy<mind){
								dz = z-e[k][i].z;
								dz = WZ*dz*dz;
								if(dz<mind){
									int d =	dx+dy+dz;
									if (d<mind){
										mind = d;
										kk = k;
										ii = i;
									}
								}
							}
						}
					}
				bm[z][y][x].id = kk;
				bm[z][y][x].p = exp(-(double)mind/BM_2VARI);
				
//				if((mind<20)&&(ii>mpp.min_c)&&(ii<N_k[kk]-mpp.min_c)){
				if(mind<20){
					bm[z][y][x].min_ab = e[kk][ii].min_ab;
					bm[z][y][x].max_ab = e[kk][ii].max_ab;
					bm[z][y][x].c_pz   = e[kk][ii].c_pz  ;
					bm[z][y][x].alpha_mean = e[kk][ii].alpha_mean;
					bm[z][y][x].beta_mean  = e[kk][ii].beta_mean ;
					bm[z][y][x].gamma_mean = e[kk][ii].gamma_mean;
					bm[z][y][x].fit_error = e[kk][ii].fit_error*0.57+0.5;
				}
				else
					bm[z][y][x].min_ab = 0; // just uniformly draw parameters
			}

#if 0		
	struct TIFF_img inter_img;
	FILE *fp;
	char name_buff[2000];
	get_TIFF ( &inter_img, height, width, 'g' );

	for (int kk = 0; kk < depth; kk++ ){
		for (int ii = 0; ii < height; ii++ )
		for (int jj = 0; jj < width; jj++ ) {
			inter_img.mono[ii][jj] = (int)(bm[kk][ii][jj].p*255);
		}
		sprintf(name_buff, "birthmap%03d.tif", kk);
		/* open image file */
		if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
		fprintf ( stderr, "cannot open file image.tif\n");
		exit ( 1 );
		}
		/* write image */
		if ( write_TIFF ( fp, &inter_img ) ) {
		fprintf ( stderr, "error writing TIFF file \n" );
		exit ( 1 );
		}
		/* close image file */
		fclose ( fp );
	}
	free_TIFF(&inter_img);
#endif
	free_img((void**)e );
	free(N_k);

	clock_t end_time=clock();
	int running_time = static_cast<int>(end_time-start_time)/CLOCKS_PER_SEC;
	cout<< "Birth map exe time: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

	return running_time;
}

double MBC_GenerateBirhMap_DiskNormal(Config_Ellipsoid *clustered_cfg, Birthmap2 ***bm, int num_clust, int width, int height, int depth, MPP_Parameters mpp)
{
	int num;
	int minz, maxz;
	CvPoint3D64f centroid, direct;
	CvPoint3D64f center[2*LINE_FIT_NUM+1];
	int id;
	int idcnt[MAX_2DOBJ_NUM];

 	for(int i=0;i<num_clust;i++){
		idcnt[i] = 0;
	}
	for(int i=0;i<clustered_cfg->num_obj;i++){
		idcnt[clustered_cfg->ellipsoid[i].obj_id]++;
	}
	int mx = 0;
	for(int i=1;i<num_clust;i++){
		if(mx<idcnt[i])
			mx = idcnt[i];
	}
   int *N_k = (int *)malloc(num_clust*sizeof(int));
	Ellipsoid1 **e=(Ellipsoid1 **)get_img(mx,num_clust,sizeof(Ellipsoid1));

	clock_t start_time=clock();
	for (int k = 1;k<num_clust;k++){
		N_k[k] = 0;
	}
	for (int j = 0;j<clustered_cfg->num_obj;j++){
		id = clustered_cfg->ellipsoid[j].obj_id;
		if(id){
			e[id][N_k[id]].x = (int)clustered_cfg->ellipsoid[j].center.x;
			e[id][N_k[id]].y = (int)clustered_cfg->ellipsoid[j].center.y;
			e[id][N_k[id]].z = (int)clustered_cfg->ellipsoid[j].center.z;
			e[id][N_k[id]].a = (int)clustered_cfg->ellipsoid[j].a;
			e[id][N_k[id]].b = (int)clustered_cfg->ellipsoid[j].b;
			e[id][N_k[id]].alpha = clustered_cfg->ellipsoid[j].alpha;
			e[id][N_k[id]].beta = clustered_cfg->ellipsoid[j].beta;
			e[id][N_k[id]].gamma = clustered_cfg->ellipsoid[j].gamma;
			N_k[id]++;
		}
	}
	
#if 0
	for (int k = 0;k<num_clust;k++){
		qsort(e[k], N_k[k], sizeof(Ellipse1), compare_zposition1);
		maxz = e[k][0].z;
		if(maxz==depth-1)
			maxz = 2*depth;
		minz = e[k][N_k[k]-1].z;
		if(minz==0)
			minz = -depth;
		for (int i = 0;i<N_k[k];i++){
			e[k][i].min_ab = min(e[k][i].a,e[k][i].b);
			e[k][i].max_ab = max(e[k][i].a,e[k][i].b);
			e[k][i].c_pz   = min(e[k][i].z-minz, maxz-e[k][i].z)+1;
			if ((i>=LINE_FIT_NUM)&&(i<N_k[k]-LINE_FIT_NUM)){
				double max1, max2, max3; 
				double min1, min2, min3; 
				min1 = min(e[k][i-1].a,e[k][i-1].b);
				max1 = max(e[k][i-1].a,e[k][i-1].b);
				min2 = min(e[k][i].a,e[k][i].b);
				max2 = max(e[k][i].a,e[k][i].b);
				min3 = min(e[k][i+1].a,e[k][i+1].b);
				max3 = max(e[k][i+1].a,e[k][i+1].b);
				e[k][i].min_ab = min3(min1,min2,min3);
				e[k][i].max_ab = max3(max1,max2,max3);
			// beta, gamma : angle of line which include centers of superellipses
				for(int z=-LINE_FIT_NUM;z<=LINE_FIT_NUM;z++){
					center[z+LINE_FIT_NUM].x = (double)e[k][i+z].x;
					center[z+LINE_FIT_NUM].y = (double)e[k][i+z].y;
					center[z+LINE_FIT_NUM].z = (double)e[k][i+z].z;
				}
				//if((k==57)&&(e[k][i].z==27))
				//	k=k;
				e[k][i].fit_error = line3D_fit(center, &centroid, &direct, 2*LINE_FIT_NUM+1);
				double alpha = get_alpha(direct.y,direct.z);
				e[k][i].alpha_mean = alpha;
				double dz = sin(-alpha)*direct.y+cos(-alpha)*direct.z;
				e[k][i].beta_mean = get_beta(direct.x,dz);
				e[k][i].gamma_mean = e[k][i].theta;
			}
			else{
				e[k][i].min_ab = 0; // just uniformly draw parameters
				e[k][i].alpha_mean = 0;
				e[k][i].beta_mean  = 0;
				e[k][i].gamma_mean = 0;
				e[k][i].fit_error = 0;
			}
		}
	}
#endif

	for (int z = 0;z<depth;z++)
		for (int y = 0;y<height;y++)
			for (int x =0;x<width;x++)
			{
				int mind = INF;
				int kk,ii;
				for (int k = 1;k<num_clust;k++)
					for (int i =0;i<N_k[k];i++)
					{
						int dx, dy, dz;
						dx = x-e[k][i].x;
						dx = dx*dx;
						if(dx<mind){
							dy = y-e[k][i].y;
							dy = dy*dy;
							if(dy<mind){
								dz = z-e[k][i].z;
								dz = WZ*dz*dz;
								if(dz<mind){
									int d =	dx+dy+dz;
									if (d<mind){
										mind = d;
										kk = k;
										ii = i;
									}
								}
							}
						}
					}
				bm[z][y][x].id = kk;
				bm[z][y][x].p = exp(-(double)mind/BM_2VARI);
				
//				if(mind<20){
				if(0){
					bm[z][y][x].min_ab = e[kk][ii].min_ab;
					bm[z][y][x].max_ab = e[kk][ii].max_ab;
					bm[z][y][x].c_pz   = e[kk][ii].c_pz  ;
					bm[z][y][x].alpha_mean = e[kk][ii].alpha_mean;
					bm[z][y][x].beta_mean  = e[kk][ii].beta_mean ;
					bm[z][y][x].gamma_mean = e[kk][ii].gamma_mean;
					bm[z][y][x].fit_error = e[kk][ii].fit_error*0.57+0.5;
				}
				else
					bm[z][y][x].min_ab = 0; // just uniformly draw parameters
			}

#if 0		
	struct TIFF_img inter_img;
	FILE *fp;
	char name_buff[2000];
	get_TIFF ( &inter_img, height, width, 'g' );

	for (int kk = 0; kk < depth; kk++ ){
		for (int ii = 0; ii < height; ii++ )
		for (int jj = 0; jj < width; jj++ ) {
			inter_img.mono[ii][jj] = (int)(bm[kk][ii][jj].p*255);
		}
		sprintf(name_buff, "birthmap%03d.tif", kk);
		/* open image file */
		if ( ( fp = fopen ( name_buff, "wb" ) ) == NULL ) {
		fprintf ( stderr, "cannot open file image.tif\n");
		exit ( 1 );
		}
		/* write image */
		if ( write_TIFF ( fp, &inter_img ) ) {
		fprintf ( stderr, "error writing TIFF file \n" );
		exit ( 1 );
		}
		/* close image file */
		fclose ( fp );
	}
	free_TIFF(&inter_img);
#endif
	free_img((void**)e );
	free(N_k);

	clock_t end_time=clock();
	double running_time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
	cout<< "Birth map exe time: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

	return running_time;
}

double MBC_ellipsoidFitMain3(unsigned char ***img, int ***pnt, Birthmap ***bm, int bm_num, ElipsoidConfig *config, Config_Ellipsoid *cfg,
			   int width, int height, int depth,MPP_Parameters mpp)
{
	int r_x = width-2*MARGIN;
	int r_y = height-2*MARGIN;
	int r_z = depth-2*MARGIN;
	int max_a = mpp.max_a;
	int max_b = mpp.max_b;
	int max_c = mpp.max_c;
	int min_a = mpp.min_a;
	int min_b = mpp.min_b;
	int min_c = mpp.min_c;
	int r_a = mpp.max_a-mpp.min_a;
	int r_b = mpp.max_b-mpp.min_b;
	int r_c = mpp.max_c-mpp.min_c;
	int kmax = height*width*depth/15;//mpp.aux; //40
	double n = mpp.n;
	double pi = CV_PI;
	double max_alpha = mpp.max_alpha;
	double min_alpha = mpp.min_alpha;
	double max_beta = mpp.max_beta;
	double min_beta = mpp.min_beta;
	double max_gamma = mpp.max_gamma;
	double min_gamma = mpp.min_gamma;
	double r_alpha = mpp.max_alpha-mpp.min_alpha;
	double r_beta  = mpp.max_beta-mpp.min_beta;
	double r_gamma = mpp.max_gamma-mpp.min_gamma;
	int b_alpha = (r_alpha>DELTA_AB_ANGLE)? 1:0;
	int b_beta  = (r_beta >DELTA_AB_ANGLE)? 1:0;
	int b_gamma = (r_gamma>DELTA_G_ANGLE)? 1:0;

	double pen  = mpp.pen; // 5 overlap penalty
	double thredh = mpp.thredh; // Eq. (5) T

	double b_zero = mpp.b_zero;
	double delta_mpp = mpp.delta; // 0.9 sigma for death step 
	double beta_mpp = mpp.beta_mpp;// 10  alpha for death step
	double F = mpp.F; // sigma and alpha decreasing and increasing factor
	int iter_num = mpp.iter_num;
	int perturb = 1;//mpp.perturbation;

	struct TIFF_img inter_img;
	int tag = 0;
	int select = 0;
	int select2 = 0;
	double iter_runtime = 0, total_runtime = 0;
	double dtmp;
	double l_th = mpp.l_th;

	get_TIFF ( &inter_img, height, width, 'g' );

	CvScalar green = CV_RGB(0,255,0);
	CvScalar white = CV_RGB(255,255,255);
	CvScalar orange = CV_RGB(255,128,0);
	CvScalar yellow = CV_RGB(255,255,0);
	Node3D *pobj;
//	int id[MAX_IMAGE_NUM] = {0,};

	LinkedList3D newlink = LinkedList3DInit();
	int newlink_num = 0;
	double overlap;
	int display_img_num = 10;

	int tmpmask_len = 2*(max_a+2+8);
	int tmpmask_depth = 2*(max_c+2+8);
	double pc_thredh = 0.4;

	clock_t start_time=clock();
	for (int iter = 1;iter <= iter_num; iter++)//25
	{
		clock_t mid_time1=clock();
		printf("************* iter = %d/%d\n",iter,iter_num);

		// birth step
		double p1;
		if(perturb==0)
			p1 = 0; // MBC
		else
			p1 = rand()/(double)(RAND_MAX); // MBC with local perturbation
		if ((iter == 1)||(p1 < 0.5)){ // birth step
			newlink_num = 0;
			int k;
			for(k=0; (k<kmax)&&(newlink_num<N_MAX3D); k++) { 
				ellipsoidObj *elipsoid = (ellipsoidObj *)malloc(sizeof(ellipsoidObj));
				elipsoid->n = (int)n;
				//elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(4)+132)+0.5);
				//elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(5)+81)+0.5);
				//elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(4)+33)+0.5);
				elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(r_x)+MARGIN)+0.5);
				elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(r_y)+MARGIN)+0.5);
				elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(r_z)+MARGIN)+0.5);
				if(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c!=255){
#if 1
					max_alpha = bm[elipsoid->z][elipsoid->y][elipsoid->x].alpha+ALPHA_MARGIN;
					min_alpha = bm[elipsoid->z][elipsoid->y][elipsoid->x].alpha-ALPHA_MARGIN;
					max_alpha = min(mpp.max_alpha,max_alpha);
					min_alpha = max(mpp.min_alpha,min_alpha);
					if(max_alpha<=min_alpha){
						max_alpha = mpp.max_alpha;
						min_alpha = mpp.min_alpha;
					}
					max_beta = bm[elipsoid->z][elipsoid->y][elipsoid->x].beta+BETA_MARGIN;
					min_beta = bm[elipsoid->z][elipsoid->y][elipsoid->x].beta-BETA_MARGIN;
					max_beta = min(mpp.max_beta,max_beta);
					min_beta = max(mpp.min_beta,min_beta);
					if(max_beta<=min_beta){
						max_beta = mpp.max_beta;
						min_beta = mpp.min_beta;
					}
#endif
					double r1_alpha = max_alpha-min_alpha;
					double r1_beta = max_beta-min_beta;
					elipsoid->alpha = rand()/(double)(RAND_MAX)*(r1_alpha)+min_alpha;
					elipsoid->beta = rand()/(double)(RAND_MAX)*(r1_beta)+min_beta;
					elipsoid->gamma = rand()/(double)(RAND_MAX)*(r_gamma)+min_gamma;
					int r_a,r_b,r_c;
#if 1
					dtmp = cos(elipsoid->alpha)*cos(elipsoid->beta);
					max_a = bm[elipsoid->z][elipsoid->y][elipsoid->x].max_ab+AB_MARGIN;
					min_a = (int)(bm[elipsoid->z][elipsoid->y][elipsoid->x].min_ab*dtmp)-AB_MARGIN;
					max_b = max_a;
					min_b = min_a;
					max_a = min(mpp.max_a,max_a);
					min_a = max(mpp.min_a,min_a);
					max_b = min(mpp.max_b,max_b);
					min_b = max(mpp.min_b,min_b);
					max_c = (int)(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c/dtmp)+C_MARGIN;
					min_c = bm[elipsoid->z][elipsoid->y][elipsoid->x].min_c-C_MARGIN;
					max_c = min(mpp.max_c,max_c);
					min_c = max(mpp.min_c,min_c);
#endif
					r_a = max_a-min_a;
					r_b = max_b-min_b;
					r_c = max_c-min_c;
					if((r_a>=0)&&(r_b>=0)&&(r_c>=0)){
						elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(r_a)+min_a)+0.5);
						elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(r_b)+min_b)+0.5);
						elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(r_c)+min_c)+0.5);
						elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
						//elipsoid->x = 135;
						//elipsoid->y = 83;
						//elipsoid->z = 35;
						//elipsoid->a = 6;
						//elipsoid->b = 6;
						//elipsoid->c = 12;
						//elipsoid->alpha = -0.52;
						//elipsoid->beta  = -0.52;
						//elipsoid->gamma = 0.52;

//						elipsoid->dataterm = ellipsoid_likelyhood33(tmpmask,tmpmask_len,tmpmask_depth,img,width,height,depth,
//								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
//								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),thredh,config);
						elipsoid->dataterm = ellipsoid_likelyhood_sn(img,width,height,depth,
								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
//						elipsoid->dataterm =  ellipsoid_fitlikelyhood(pnt, width, height,depth,elipsoid->y,elipsoid->x,elipsoid->z,
//								elipsoid->a,elipsoid->b,elipsoid->c, elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,
//								elipsoid->obj_id, &(elipsoid->dist),pc_thredh,config,0);
//							printf("%1.3f, %1.3f, %1.3f\n",elipsoid->dataterm,elipsoid->dist,dtmp);
						if(elipsoid->dataterm<=l_th){
//				printf("k = %d ",k);
							elipsoid->gcut_weight = (int)(500.*(1+elipsoid->dataterm));
							drawSuperEllipsoid33(elipsoid, tmpmask_len, tmpmask_len, tmpmask_depth, elipsoid->a, elipsoid->b, elipsoid->c,
								elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
							if (iter == 1)
							{
								overlap = elipsoid_overlap_check_all(elipsoid, config->elink, config->mp_num);
								if(!overlap){
									config->mp_num++;
									LinkedList3DInsert(config->elink, config->mp_num, elipsoid); 
									newlink_num++;
	//			printf("k = %d, new_obj = %d\n",k,newlink_num);
								}
								else{
									free_volume((void***)elipsoid->ptmpmask);
									elipsoid->ptmpmask = NULL;
									free(elipsoid);
//				printf("k = %d\n",k);
								}
							}
							else
							{
								overlap = elipsoid_overlap_check_all(elipsoid, newlink, newlink_num);
								if(!overlap){
									newlink_num++;
									LinkedList3DInsert(newlink,newlink_num, elipsoid); 
	//			printf("k = %d, new_obj = %d\n",k,newlink_num);
								}
								else{
									free_volume((void***)elipsoid->ptmpmask);
									elipsoid->ptmpmask = NULL;
									free(elipsoid);
//				printf("k = %d\n",k);
								}
							}
						}
						else{
							free(elipsoid);
						}
					}
					else{
						free(elipsoid);
					}
				}
				else{
					free(elipsoid);
				}
			}
			printf("k = %d, new_obj = %d\n",k,newlink_num);
		}
		else {	// local perturbation
			pobj = config->elink;
			newlink_num = 0;
			int k;
			for(k=0; k<config->mp_num; k++) { 
				pobj = pobj->next;
				ellipsoidObj *elipsoid = (ellipsoidObj *)malloc(sizeof(ellipsoidObj));
				elipsoid->n = (int)n;
				int valid = 0;
#if 1
				elipsoid->x = pobj->index->x;
				elipsoid->y = pobj->index->y;
				elipsoid->z = pobj->index->z;
				elipsoid->a = pobj->index->a;
				elipsoid->b = pobj->index->b;
				elipsoid->c = pobj->index->c;
				elipsoid->alpha = pobj->index->alpha;
				elipsoid->beta  = pobj->index->beta;
				elipsoid->gamma = pobj->index->gamma;
				elipsoid->obj_id = pobj->index->obj_id;
				p1 = rand()/(double)(RAND_MAX); // MBC with local perturbation
				if(p1<0.1){
					do{
						elipsoid->x = pobj->index->x;
						elipsoid->x += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
						if((elipsoid->x>=MARGIN)&&(elipsoid->x<=width-MARGIN)){
							if(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c!=0){
								valid = 1;
								elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
							}
						}
					}while(!valid);
				}
				else if(p1<0.2){
					do{
						elipsoid->y = pobj->index->y;
						elipsoid->y += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
						if((elipsoid->y>=MARGIN)&&(elipsoid->y<=height-MARGIN)){
						if(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c!=0){
								valid = 1;
								elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
							}
						}
					}while(!valid);
				}
				else if(p1<0.3){
					do{
						elipsoid->z = pobj->index->z;
						elipsoid->z += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_TRANSLATION)-DELTAZ_TRANSLATION)+0.5);
						if((elipsoid->z>=MARGIN)&&(elipsoid->z<=depth-MARGIN)){
							if(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c!=0){
								valid = 1;
								elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
							}
						}
					}while(!valid);
				}
				else if(p1<0.4){
					do{
						elipsoid->a = pobj->index->a;
						if(r_a>DELTA_DILATION)
							elipsoid->a += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
						else
							elipsoid->a += (int)floor((rand()/(double)(RAND_MAX)*(2*r_a)-r_a)+0.5);
						if((elipsoid->a>=mpp.min_a)&&(elipsoid->a<=mpp.max_a))
							valid = 1;
					}while(!valid);
				}
				else if(p1<0.5){
					do{
						elipsoid->b = pobj->index->b;
						if(r_b>DELTA_DILATION)
							elipsoid->b += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
						else
							elipsoid->b += (int)floor((rand()/(double)(RAND_MAX)*(2*r_b)-r_b)+0.5);
						if((elipsoid->b>=mpp.min_b)&&(elipsoid->b<=mpp.max_b))
							valid = 1;
					}while(!valid);
				}
				else if(p1<0.6){
					do{
						elipsoid->c = pobj->index->c;
						if(r_c>DELTAZ_DILATION)
							elipsoid->c += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_DILATION)-DELTAZ_DILATION)+0.5);
						else
							elipsoid->c += (int)floor((rand()/(double)(RAND_MAX)*(2*r_c)-r_c)+0.5);
						if((elipsoid->c>=mpp.min_c)&&(elipsoid->c<=mpp.max_c))
							valid = 1;
					}while(!valid);
				}
				else if(p1<0.7){
					do{
						elipsoid->alpha = pobj->index->alpha;
						if(b_alpha)
							elipsoid->alpha += rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
						else
							elipsoid->alpha += rand()/(double)(RAND_MAX)*(2*r_alpha)-r_alpha;
						if((elipsoid->alpha>=mpp.min_alpha)&&(elipsoid->alpha<=mpp.max_alpha))
							valid = 1;

					}while(!valid);
				}
				else if(p1<0.8){
					do{
						elipsoid->beta = pobj->index->beta;
						if(b_beta)
							elipsoid->beta += rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
						else
							elipsoid->beta += rand()/(double)(RAND_MAX)*(2*r_beta)-r_beta;
						if((elipsoid->beta>=mpp.min_beta)&&(elipsoid->beta<=mpp.max_beta))
							valid = 1;
					}while(!valid);
				}
				else{
					do{
						elipsoid->gamma = pobj->index->gamma;
						if(b_gamma)
							elipsoid->gamma += rand()/(double)(RAND_MAX)*(2*DELTA_G_ANGLE)-DELTA_G_ANGLE;
						else
							elipsoid->gamma += rand()/(double)(RAND_MAX)*(2*r_gamma)-r_gamma;
						if((elipsoid->gamma>=mpp.min_gamma)&&(elipsoid->gamma<=mpp.max_gamma))
							valid = 1;
					}while(!valid);
				}
#else
				do{
					elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_TRANSLATION)-DELTAZ_TRANSLATION)+0.5);
					elipsoid->x += pobj->index->x;
					elipsoid->y += pobj->index->y;
					elipsoid->z += pobj->index->z;
					if((elipsoid->x>=MARGIN)&&(elipsoid->x<=width-MARGIN)
						&&(elipsoid->y>=MARGIN)&&(elipsoid->y<=height-MARGIN)
						&&(elipsoid->z>=MARGIN)&&(elipsoid->z<=depth-MARGIN)){
						if(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c!=0){
							if(b_alpha)
								elipsoid->alpha = rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
							else
								elipsoid->alpha = rand()/(double)(RAND_MAX)*(2*r_alpha)-r_alpha;
							if(b_beta)
								elipsoid->beta = rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
							else
								elipsoid->beta = rand()/(double)(RAND_MAX)*(2*r_beta)-r_beta;
							if(b_gamma)
								elipsoid->gamma = rand()/(double)(RAND_MAX)*(2*DELTA_G_ANGLE)-DELTA_G_ANGLE;
							else
								elipsoid->gamma = rand()/(double)(RAND_MAX)*(2*r_gamma)-r_gamma;
							elipsoid->alpha += pobj->index->alpha;
							elipsoid->beta  += pobj->index->beta;
							elipsoid->gamma += pobj->index->gamma;
							if((elipsoid->alpha>=mpp.min_alpha)&&(elipsoid->alpha<=mpp.max_alpha)
								&&(elipsoid->beta>=mpp.min_beta)&&(elipsoid->beta<=mpp.max_beta)
								&&(elipsoid->gamma>=mpp.min_gamma)&&(elipsoid->gamma<=mpp.max_gamma)){
#if 1
								dtmp = cos(elipsoid->alpha)*cos(elipsoid->beta);
								max_a = bm[elipsoid->z][elipsoid->y][elipsoid->x].max_ab+AB_MARGIN;
								min_a = (int)(bm[elipsoid->z][elipsoid->y][elipsoid->x].min_ab*dtmp)-AB_MARGIN;
								max_a = min(mpp.max_a,max_a);
								min_a = max(mpp.min_a,min_a);
								max_b = max_a;
								min_b = min_a;
								max_c = (int)(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c/dtmp)+C_MARGIN;
								min_c = bm[elipsoid->z][elipsoid->y][elipsoid->x].min_c-C_MARGIN;
								max_c = min(mpp.max_c,max_c);
								min_c = max(mpp.min_c,min_c);
#endif
								tmp = max_a-min_a;
								if(tmp>DELTA_DILATION)
									elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
								else
									elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(2*tmp)-tmp)+0.5);
								tmp = max_b-min_b;
								if(tmp>DELTA_DILATION)
									elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
								else
									elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(2*tmp)-tmp)+0.5);
								tmp = max_c-min_c;
								if(tmp>DELTAZ_DILATION)
									elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_DILATION)-DELTAZ_DILATION)+0.5);
								else
									elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(2*tmp)-tmp)+0.5);
								elipsoid->a += pobj->index->a;
								elipsoid->b += pobj->index->b;
								elipsoid->c += pobj->index->c;
								if((elipsoid->a>=mpp.min_a)&&(elipsoid->a<=mpp.max_a)
									&&(elipsoid->b>=mpp.min_b)&&(elipsoid->b<=mpp.max_b)
									&&(elipsoid->c>=mpp.min_c)&&(elipsoid->c<=mpp.max_c)){
									valid = 1;
									elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
								}
								else
									valid = 0;
							}
							else
								valid = 0;
						}
						else
							valid = 0;
					}
					else
						valid = 0;
				}while(!valid);
#endif

//						elipsoid->dataterm = ellipsoid_likelyhood33(tmpmask,tmpmask_len,tmpmask_depth,img,width,height,depth,
//								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
//								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),thredh,config);
						elipsoid->dataterm = ellipsoid_likelyhood_sn(img,width,height,depth,
								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
//						elipsoid->dataterm =  ellipsoid_fitlikelyhood(pnt, width, height,depth,elipsoid->y,elipsoid->x,elipsoid->z,
//								elipsoid->a,elipsoid->b,elipsoid->c, elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,
//								elipsoid->obj_id, &(elipsoid->dist),pc_thredh,config,0);
//				printf("%1.3f, %1.3f, %1.3f, %d\n",data,dist,elipsoid->dataterm,elipsoid->c);
				if(elipsoid->dataterm<=l_th){
					elipsoid->gcut_weight = (int)(500.*(1+elipsoid->dataterm));

					drawSuperEllipsoid33(elipsoid, tmpmask_len, tmpmask_len, tmpmask_depth,	elipsoid->a, elipsoid->b, elipsoid->c, 
						elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
					overlap = elipsoid_overlap_check_all(elipsoid,newlink,newlink_num);
					if(!overlap){
						newlink_num++;
						LinkedList3DInsert(newlink, newlink_num, elipsoid); 
					}
					else{
						free_volume((void***)elipsoid->ptmpmask);
						elipsoid->ptmpmask = NULL;
						free(elipsoid);
					}
				}
				else{
					free(elipsoid);
				}
			}
			printf("k = %d, new_obj = %d\n",k,newlink_num);
		}


		// deadth step
		if((iter!=1)&&(newlink_num)){
//			printf("Graph cut start\n");
			graph_cut33(config,newlink,newlink_num);
//			printf("Graph cut end\n");

			pobj = config->elink->next;
			for(int i = 1;i<=config->mp_num;i++){
				if (pobj->index->my_type == 0){ // delete from exist config
					pobj = pobj->next;
					Qdelete_elipsoid3(config, i);
					i--;
				}
				else							// keep in config
					pobj = pobj->next;
			}
			pobj = newlink->next;
			for(int i = 1;i<=newlink_num;i++){
				if (pobj->index->my_type == 0){ // delete from newlink
					pobj = pobj->next;
					Qdelete_elipsoid3(newlink, i, &newlink_num);
					i--;
				}
				else{							// move from newlink to config
					config->mp_num++;
					LinkedList3DInsert(config->elink,config->mp_num, pobj->index); 
					pobj = pobj->next;
					LinkedList3DDelete(newlink, i);//do not delete object itself
					newlink_num--;
					i--;
				}
			}
		}
		printf("obj num = %d\n",config->mp_num);

		clock_t mid_time2=clock();
		cout<< "Iteration running time is: "<<static_cast<double>(mid_time2-mid_time1)/CLOCKS_PER_SEC<<"s"<<endl;
		pobj = config->elink;
		int exe_time = static_cast<int>(mid_time2-start_time); //ms
		config->exe_time[iter] = exe_time;
		config->num_obj[iter] = config->mp_num;
		config->energy[iter] = 0;
		for (int k = 1;k<=config->mp_num;k++){
			pobj = pobj->next;
//			prior = C_prior_elipsoid3(config, k, pen);
//			dtmp = pobj->index->dataterm
//				+pobj->index->inter_e*mpp.intra_pen/2.;
			config->energy[iter] += pobj->index->dataterm;
		}
		printf("E = %1.5f\n",config->energy[iter]);
	}
	clock_t end_time=clock();
	double running_time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
	cout<< "Running time is: "<< running_time<<"s"<<endl;

	int z;
	pobj = config->elink;
	cfg->num_obj=0;
	cfg->total_E = 0;
	for (int k = 1;k<=config->mp_num;k++){
		pobj = pobj->next;
		z = pobj->index->z;
//		prior = C_prior_elipsoid3(config, k, pen);
		cfg->ellipsoid[cfg->num_obj].center.x = pobj->index->x;
		cfg->ellipsoid[cfg->num_obj].center.y = pobj->index->y;
		cfg->ellipsoid[cfg->num_obj].center.z = pobj->index->z;
		cfg->ellipsoid[cfg->num_obj].a = pobj->index->a;
		cfg->ellipsoid[cfg->num_obj].b = pobj->index->b;
		cfg->ellipsoid[cfg->num_obj].c = pobj->index->c;
		cfg->ellipsoid[cfg->num_obj].alpha = pobj->index->alpha;
		cfg->ellipsoid[cfg->num_obj].beta = pobj->index->beta;
		cfg->ellipsoid[cfg->num_obj].gamma = pobj->index->gamma;
		cfg->ellipsoid[cfg->num_obj].single_E = pobj->index->dataterm;
		cfg->ellipsoid[cfg->num_obj].e0 = pobj->index->dist;
		cfg->ellipsoid[cfg->num_obj].multiple_E = 0;
		cfg->ellipsoid[cfg->num_obj].obj_id = pobj->index->obj_id;
		cfg->ellipsoid[cfg->num_obj].n = pobj->index->n; // ELLIPSE_POWER
		cfg->total_E += pobj->index->dataterm;
		cfg->num_obj++;
	}
	printf("obj num = %d\n",cfg->num_obj);

	return running_time;
}
int sample_normal(double *sample, double mean, double stdev, double min, double max)
{
	int k = 0;
	do{
		*sample = normal()*stdev + mean;
		if((*sample<=max)&&(*sample>=min))
			break;
		else
			k++;
	}while(k<1000);
	if(k==1000)
		return 0;
	else
		return 1;
}

#define EPSILON_R 2

double MBC_ellipsoidFitMain3_Normal(unsigned char ***img, int ***pnt, Birthmap2 ***bm, ElipsoidConfig *config, Config_Ellipsoid *cfg,
			   int width, int height, int depth,MPP_Parameters mpp)
{
	int r_x = width-2*MARGIN;
	int r_y = height-2*MARGIN;
	int r_z = depth-2*MARGIN;
	int max_a = mpp.max_a;
	int max_b = mpp.max_b;
	int max_c = mpp.max_c;
	int min_a = mpp.min_a;
	int min_b = mpp.min_b;
	int min_c = mpp.min_c;
	int r_a = mpp.max_a-mpp.min_a;
	int r_b = mpp.max_b-mpp.min_b;
	int r_c = mpp.max_c-mpp.min_c;
	double max_r = sqrt((double)(max_a*max_a+max_b*max_b));
	int kmax = height*width*depth;//height*width*depth/15;//mpp.aux; //40
	double n = mpp.n;
	double pi = CV_PI;
	double max_alpha = mpp.max_alpha;
	double min_alpha = mpp.min_alpha;
	double max_beta = mpp.max_beta;
	double min_beta = mpp.min_beta;
	double max_gamma = mpp.max_gamma;
	double min_gamma = mpp.min_gamma;
	double r_alpha = mpp.max_alpha-mpp.min_alpha;
	double r_beta  = mpp.max_beta-mpp.min_beta;
	double r_gamma = mpp.max_gamma-mpp.min_gamma;
	int b_alpha = (r_alpha>DELTA_AB_ANGLE)? 1:0;
	int b_beta  = (r_beta >DELTA_AB_ANGLE)? 1:0;
	int b_gamma = (r_gamma>DELTA_G_ANGLE)? 1:0;

	double pen  = mpp.pen; // 5 overlap penalty
	double thredh = mpp.thredh; // Eq. (5) T

	double b_zero = mpp.b_zero;
	double delta_mpp = mpp.delta; // 0.9 sigma for death step 
	double beta_mpp = mpp.beta_mpp;// 10  alpha for death step
	double F = mpp.F; // sigma and alpha decreasing and increasing factor
	int iter_num = mpp.iter_num;
	int perturb = 1;//mpp.perturbation;
	int test[60] = {0,};

	struct TIFF_img inter_img;
	int tag = 0;
	int select = 0;
	int select2 = 0;
	double iter_runtime = 0, total_runtime = 0;
	int tmp;
	double dtmp;
	double l_th = mpp.l_th;
	int num_overlap;

	get_TIFF ( &inter_img, height, width, 'g' );

	CvScalar green = CV_RGB(0,255,0);
	CvScalar white = CV_RGB(255,255,255);
	CvScalar orange = CV_RGB(255,128,0);
	CvScalar yellow = CV_RGB(255,255,0);
	Node3D *pobj;
//	int id[MAX_IMAGE_NUM] = {0,};

	LinkedList3D newlink = LinkedList3DInit();
	int newlink_num = 0;
	double overlap;
	int display_img_num = 10;

	int tmpmask_len = 2*(max_a+2+8);
	int tmpmask_depth = 2*(max_c+2+8);
	double pc_thredh = 0.4;
	int rst;
	double fit_error_max = 0, fit_error_min = 10000, fit_error_sum = 0;
	int fit_error_num = 0;
	double **cos_phi = (double **)get_img((int)(CV_PI*1000)+1,(int)(CV_PI*1000)+1,sizeof(double));
	double *tan_phi = (double *)malloc(sizeof(double)*1001);
	double **r_p_M = (double **)get_img(mpp.max_a-mpp.min_a+1,mpp.max_b-mpp.min_b+1,sizeof(double));

	for(int i = 0;i<(int)(CV_PI*1000)+1;i++)
		for(int j = 0;j<(int)(CV_PI*1000)+1;j++)
			cos_phi[i][j] = cos((double)i/1000-PI_2)*cos((double)j/1000-PI_2);
	for(int i = 0;i<mpp.max_a-mpp.min_a+1;i++)
		for(int j = 0;j<mpp.max_b-mpp.min_b+1;j++)
			r_p_M[i][j] = hypot((double)(i+mpp.min_a),(double)(j+mpp.min_b));
	for(int i = 0;i<1001;i++){
		double cos_phi = (double)i/1000.;
		tan_phi[i] = sqrt(1/(cos_phi*cos_phi)-1);
	}

	clock_t start_time=clock();
	for (int iter = 1;iter <= iter_num; iter++)//25
	{
		clock_t mid_time1=clock();
		printf("************* iter = %d/%d\n",iter,iter_num);

		// birth step
		double p1;
		if(perturb==0)
			p1 = 0; // MBC
		else
			p1 = rand()/(double)(RAND_MAX); // MBC with local perturbation
		if ((iter == 1)||(p1 < 0.5)){ // birth step
			newlink_num = 0;
			int k;
			num_overlap = 0;
			for(k=0; (k<kmax)&&(num_overlap<5); k++) { 
				ellipsoidObj *elipsoid = (ellipsoidObj *)malloc(sizeof(ellipsoidObj));
				elipsoid->type = mpp.type;
				elipsoid->n = (int)n;
				//elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(7)+40)+0.5);
				//elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(7)+29)+0.5);
				//elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(25)+22)+0.5);
				//elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(7)+40)+0.5);
				//elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(7)+29)+0.5);
				//elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(7)+30)+0.5);
				elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(r_x)+MARGIN)+0.5);
				elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(r_y)+MARGIN)+0.5);
				elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(r_z)+MARGIN)+0.5);
				p1 = rand()/(double)(RAND_MAX);
				if(bm[elipsoid->z][elipsoid->y][elipsoid->x].p >= p1){
					Birthmap2 bmap = bm[elipsoid->z][elipsoid->y][elipsoid->x];
					if((mpp.aux==1)||(mpp.aux==3)){
						if(bmap.min_ab!=0){
							rst = sample_normal(&(elipsoid->alpha),bmap.alpha_mean,/*bmap.fit_error**/BM_STDEV2,min_alpha,max_alpha);
							//if(fit_error_min>bmap.fit_error)
							//	fit_error_min=bmap.fit_error;
							//if(fit_error_max<bmap.fit_error)
							//	fit_error_max=bmap.fit_error;
							//fit_error_num++;
							//fit_error_sum += bmap.beta_mean;
							//printf("fit_error = %1.3f ",fit_error_sum/(double)fit_error_num);
							if(!rst){
								//printf("sample alpha error!!\n");
								free(elipsoid);
								continue;
							}
							rst = sample_normal(&(elipsoid->beta),bmap.beta_mean,/*bmap.fit_error**/BM_STDEV2,min_beta,max_beta);
							if(!rst){
								//printf("sample beta error!!\n");
								free(elipsoid);
								continue;
							}
							//rst = sample_normal(&(elipsoid->gamma),bmap.gamma_mean,BM_STDEV2,min_gamma,max_gamma);
							//if(!rst){
							//	//printf("sample gamma error!!\n");
							//	free(elipsoid);
							//	continue;
							//}
							elipsoid->gamma = rand()/(double)(RAND_MAX)*(r_gamma)+min_gamma;
							elipsoid->alpha_mean = bmap.alpha_mean;
							elipsoid->beta_mean = bmap.beta_mean;
							elipsoid->fit_error = bmap.fit_error*BM_STDEV2;
						}
						else{
							elipsoid->alpha = rand()/(double)(RAND_MAX)*(r_alpha)+min_alpha;
							elipsoid->beta = rand()/(double)(RAND_MAX)*(r_beta)+min_beta;
							elipsoid->gamma = rand()/(double)(RAND_MAX)*(r_gamma)+min_gamma;
							elipsoid->alpha_mean = 0;
							elipsoid->beta_mean = 0;
							elipsoid->fit_error = 0;
						}
					}
					else{
						elipsoid->alpha = rand()/(double)(RAND_MAX)*(r_alpha)+min_alpha;
						elipsoid->beta = rand()/(double)(RAND_MAX)*(r_beta)+min_beta;
						elipsoid->gamma = rand()/(double)(RAND_MAX)*(r_gamma)+min_gamma;
					}

					if((mpp.aux==2)||(mpp.aux==3)){
						if(bmap.min_ab!=0){
							double cosphi = cos_phi[(int)((elipsoid->alpha+PI_2)*1000)][(int)((elipsoid->beta+PI_2)*1000)];
							double pm = bmap.min_ab*cosphi;
							double pM = bmap.max_ab;
							double mean = (pm+pM)/2.0;
							double stdev = stdev = (pM-pm+EPSILON_R)/2.0*BM_STDEV3;
							double sample;
							rst = sample_normal(&sample,mean,stdev,min_a,max_a);
							if(!rst){
								//printf("sample a error!!\n");
								free(elipsoid);
								continue;
							}
							elipsoid->a = (int)(sample+0.5);
							rst = sample_normal(&sample,mean,stdev,min_b,max_b);
							if(!rst){
								//printf("sample b error!!\n");
								free(elipsoid);
								continue;
							}
							elipsoid->b = (int)(sample+0.5);
							pM = bmap.c_pz/cosphi;
							pm = pM-r_p_M[elipsoid->a-mpp.min_a][elipsoid->b-mpp.min_b]*tan_phi[(int)(cosphi*1000)];
							//pm = pM-max_r*tan_phi;
							//double pm = bmap.c_pz*cos_phi;
							mean = (pm+pM)/2.0;
							stdev = (pM-pm+EPSILON_R)/2.0*BM_STDEV3;
							rst = sample_normal(&sample,mean,stdev,min_c,max_c);
							if(!rst){
								//printf("sample c error!!\n");
								free(elipsoid);
								continue;
							}
							elipsoid->c = (int)(sample+0.5);
						}
						else{
							elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(r_a)+min_a)+0.5);
							elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(r_b)+min_b)+0.5);
							elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(r_c)+min_c)+0.5);
						}
					}
					else{
						elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(r_a)+min_a)+0.5);
						elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(r_b)+min_b)+0.5);
						elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(r_c)+min_c)+0.5);
					}
					elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
					//elipsoid->x = 35;
					//elipsoid->y = 127;
					//elipsoid->z = 15;
					//elipsoid->a = 5;
					//elipsoid->b = 9;
					//elipsoid->c = 13;
					//elipsoid->alpha = -0.52;
					//elipsoid->beta  = 0.52;
					//elipsoid->gamma = 1.2566;
					if(elipsoid->type == OBJ3D_ELLIPSOID)
						elipsoid->dataterm = ellipsoid_likelyhood_sn(img,width,height,depth,
								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
					else
						elipsoid->dataterm = ellipsoid_likelyhood_cylinder(img,width,height,depth,
								elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
								elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
					if(elipsoid->dataterm<=l_th){
						test[elipsoid->z]++;
						elipsoid->gcut_weight = (int)(500.*(1+elipsoid->dataterm));
						double dr = max(elipsoid->a,elipsoid->b);
						double dx = 2.*(fabs((double)elipsoid->c*sin(elipsoid->beta))+dr)+1.;
						double dy = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*sin(elipsoid->alpha))+dr)+1.;
						double dz = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*cos(elipsoid->alpha))+dr)+1.;
						if(elipsoid->type == OBJ3D_ELLIPSOID)
							drawSuperEllipsoid33(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
								elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						else if(elipsoid->type == OBJ3D_CYLINDER)
							drawSuperEllipsoid_Cylinder(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
								elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						if (iter == 1)
						{
							if(elipsoid->type != OBJ3D_DISK)
								overlap = elipsoid_overlap_check_all(elipsoid, config->elink, config->mp_num);
							else
								overlap = elipsoid_overlap_check_all_disk(elipsoid, config->elink, config->mp_num);
							if(!overlap){
								config->mp_num++;
								LinkedList3DInsert(config->elink, config->mp_num, elipsoid); 
								newlink_num++;
//			printf("k = %d, new_obj = %d\n",k,newlink_num);
							}
							else if(elipsoid->type != OBJ3D_DISK){
								free_volume((void***)elipsoid->ptmpmask);
								elipsoid->ptmpmask = NULL;
								free(elipsoid);
							}
							else {
								free(elipsoid);
							}
						}
						else
						{
							if(elipsoid->type != OBJ3D_DISK)
								overlap = elipsoid_overlap_check_all(elipsoid, newlink, newlink_num);
							else
								overlap = elipsoid_overlap_check_all_disk(elipsoid, newlink, newlink_num);
							if(!overlap){
								newlink_num++;
								LinkedList3DInsert(newlink,newlink_num, elipsoid); 
//			printf("k = %d, new_obj = %d\n",k,newlink_num);
							}
							else if(elipsoid->type != OBJ3D_DISK){
								free_volume((void***)elipsoid->ptmpmask);
								elipsoid->ptmpmask = NULL;
								free(elipsoid);
							}
							else {
								free(elipsoid);
							}
						}
					}
					else{
						free(elipsoid);
					}
				}
				else{
					free(elipsoid);
				}
			}
			printf("k = %d, new_obj = %d\n",k,newlink_num);
		}
		else {	// local perturbation
			pobj = config->elink;
			newlink_num = 0;
			int k;
			for(k=0; k<config->mp_num; k++) { 
				pobj = pobj->next;
				ellipsoidObj *elipsoid = (ellipsoidObj *)malloc(sizeof(ellipsoidObj));
				elipsoid->type = pobj->index->type;
				elipsoid->n = (int)n;
				elipsoid->ptmpmask = NULL;
				int valid = 0;
#if 1
				elipsoid->x = pobj->index->x;
				elipsoid->y = pobj->index->y;
				elipsoid->z = pobj->index->z;
				elipsoid->a = pobj->index->a;
				elipsoid->b = pobj->index->b;
				elipsoid->c = pobj->index->c;
				elipsoid->alpha = pobj->index->alpha;
				elipsoid->beta  = pobj->index->beta;
				elipsoid->gamma = pobj->index->gamma;
				elipsoid->obj_id = pobj->index->obj_id;
				elipsoid->alpha_mean = pobj->index->alpha_mean;
				elipsoid->beta_mean = pobj->index->beta_mean;
				elipsoid->fit_error = pobj->index->fit_error;
				p1 = rand()/(double)(RAND_MAX); // MBC with local perturbation
				if(p1<0.1){
					do{
						elipsoid->x = pobj->index->x;
						elipsoid->x += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
						if((elipsoid->x>=MARGIN)&&(elipsoid->x<=width-MARGIN)){
							valid = 1;
							elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
						}
					}while(!valid);
				}
				else if(p1<0.2){
					do{
						elipsoid->y = pobj->index->y;
						elipsoid->y += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
						if((elipsoid->y>=MARGIN)&&(elipsoid->y<=height-MARGIN)){
							valid = 1;
							elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
						}
					}while(!valid);
				}
				else if(p1<0.3){
					do{
						elipsoid->z = pobj->index->z;
						elipsoid->z += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_TRANSLATION)-DELTAZ_TRANSLATION)+0.5);
						if((elipsoid->z>=MARGIN)&&(elipsoid->z<=depth-MARGIN)){
							valid = 1;
							elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
						}
					}while(!valid);
				}
				else if(p1<0.4){
					do{
						elipsoid->a = pobj->index->a;
						if(r_a>DELTA_DILATION)
							elipsoid->a += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
						else
							elipsoid->a += (int)floor((rand()/(double)(RAND_MAX)*(2*r_a)-r_a)+0.5);
						if((elipsoid->a>=mpp.min_a)&&(elipsoid->a<=mpp.max_a))
							valid = 1;
					}while(!valid);
				}
				else if(p1<0.5){
					do{
						elipsoid->b = pobj->index->b;
						if(r_b>DELTA_DILATION)
							elipsoid->b += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
						else
							elipsoid->b += (int)floor((rand()/(double)(RAND_MAX)*(2*r_b)-r_b)+0.5);
						if((elipsoid->b>=mpp.min_b)&&(elipsoid->b<=mpp.max_b))
							valid = 1;
					}while(!valid);
				}
				else if(p1<0.6){
					do{
						elipsoid->c = pobj->index->c;
						if(r_c>DELTAZ_DILATION)
							elipsoid->c += (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_DILATION)-DELTAZ_DILATION)+0.5);
						else
							elipsoid->c += (int)floor((rand()/(double)(RAND_MAX)*(2*r_c)-r_c)+0.5);
						if((elipsoid->c>=mpp.min_c)&&(elipsoid->c<=mpp.max_c))
							valid = 1;
					}while(!valid);
				}
				else if(p1<0.7){
					do{
						elipsoid->alpha = pobj->index->alpha;
						if(b_alpha)
							elipsoid->alpha += rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
						else
							elipsoid->alpha += rand()/(double)(RAND_MAX)*(2*r_alpha)-r_alpha;
						if((elipsoid->alpha>=mpp.min_alpha)&&(elipsoid->alpha<=mpp.max_alpha))
							valid = 1;

					}while(!valid);
				}
				else if(p1<0.8){
					do{
						elipsoid->beta = pobj->index->beta;
						if(b_beta)
							elipsoid->beta += rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
						else
							elipsoid->beta += rand()/(double)(RAND_MAX)*(2*r_beta)-r_beta;
						if((elipsoid->beta>=mpp.min_beta)&&(elipsoid->beta<=mpp.max_beta))
							valid = 1;
					}while(!valid);
				}
				else{
					do{
						elipsoid->gamma = pobj->index->gamma;
						if(b_gamma)
							elipsoid->gamma += rand()/(double)(RAND_MAX)*(2*DELTA_G_ANGLE)-DELTA_G_ANGLE;
						else
							elipsoid->gamma += rand()/(double)(RAND_MAX)*(2*r_gamma)-r_gamma;
						if((elipsoid->gamma>=mpp.min_gamma)&&(elipsoid->gamma<=mpp.max_gamma))
							valid = 1;
					}while(!valid);
				}
#else
				do{
					elipsoid->x = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elipsoid->y = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_TRANSLATION)-DELTA_TRANSLATION)+0.5);
					elipsoid->z = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_TRANSLATION)-DELTAZ_TRANSLATION)+0.5);
					elipsoid->x += pobj->index->x;
					elipsoid->y += pobj->index->y;
					elipsoid->z += pobj->index->z;
					if((elipsoid->x>=MARGIN)&&(elipsoid->x<=width-MARGIN)
						&&(elipsoid->y>=MARGIN)&&(elipsoid->y<=height-MARGIN)
						&&(elipsoid->z>=MARGIN)&&(elipsoid->z<=depth-MARGIN)){
						if(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c!=0){
							if(b_alpha)
								elipsoid->alpha = rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
							else
								elipsoid->alpha = rand()/(double)(RAND_MAX)*(2*r_alpha)-r_alpha;
							if(b_beta)
								elipsoid->beta = rand()/(double)(RAND_MAX)*(2*DELTA_AB_ANGLE)-DELTA_AB_ANGLE;
							else
								elipsoid->beta = rand()/(double)(RAND_MAX)*(2*r_beta)-r_beta;
							if(b_gamma)
								elipsoid->gamma = rand()/(double)(RAND_MAX)*(2*DELTA_G_ANGLE)-DELTA_G_ANGLE;
							else
								elipsoid->gamma = rand()/(double)(RAND_MAX)*(2*r_gamma)-r_gamma;
							elipsoid->alpha += pobj->index->alpha;
							elipsoid->beta  += pobj->index->beta;
							elipsoid->gamma += pobj->index->gamma;
							if((elipsoid->alpha>=mpp.min_alpha)&&(elipsoid->alpha<=mpp.max_alpha)
								&&(elipsoid->beta>=mpp.min_beta)&&(elipsoid->beta<=mpp.max_beta)
								&&(elipsoid->gamma>=mpp.min_gamma)&&(elipsoid->gamma<=mpp.max_gamma)){
#if 1
								dtmp = cos(elipsoid->alpha)*cos(elipsoid->beta);
								max_a = bm[elipsoid->z][elipsoid->y][elipsoid->x].max_ab+AB_MARGIN;
								min_a = (int)(bm[elipsoid->z][elipsoid->y][elipsoid->x].min_ab*dtmp)-AB_MARGIN;
								max_a = min(mpp.max_a,max_a);
								min_a = max(mpp.min_a,min_a);
								max_b = max_a;
								min_b = min_a;
								max_c = (int)(bm[elipsoid->z][elipsoid->y][elipsoid->x].max_c/dtmp)+C_MARGIN;
								min_c = bm[elipsoid->z][elipsoid->y][elipsoid->x].min_c-C_MARGIN;
								max_c = min(mpp.max_c,max_c);
								min_c = max(mpp.min_c,min_c);
#endif
								tmp = max_a-min_a;
								if(tmp>DELTA_DILATION)
									elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
								else
									elipsoid->a = (int)floor((rand()/(double)(RAND_MAX)*(2*tmp)-tmp)+0.5);
								tmp = max_b-min_b;
								if(tmp>DELTA_DILATION)
									elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTA_DILATION)-DELTA_DILATION)+0.5);
								else
									elipsoid->b = (int)floor((rand()/(double)(RAND_MAX)*(2*tmp)-tmp)+0.5);
								tmp = max_c-min_c;
								if(tmp>DELTAZ_DILATION)
									elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(2*DELTAZ_DILATION)-DELTAZ_DILATION)+0.5);
								else
									elipsoid->c = (int)floor((rand()/(double)(RAND_MAX)*(2*tmp)-tmp)+0.5);
								elipsoid->a += pobj->index->a;
								elipsoid->b += pobj->index->b;
								elipsoid->c += pobj->index->c;
								if((elipsoid->a>=mpp.min_a)&&(elipsoid->a<=mpp.max_a)
									&&(elipsoid->b>=mpp.min_b)&&(elipsoid->b<=mpp.max_b)
									&&(elipsoid->c>=mpp.min_c)&&(elipsoid->c<=mpp.max_c)){
									valid = 1;
									elipsoid->obj_id = bm[elipsoid->z][elipsoid->y][elipsoid->x].id;
								}
								else
									valid = 0;
							}
							else
								valid = 0;
						}
						else
							valid = 0;
					}
					else
						valid = 0;
				}while(!valid);
#endif

				if(elipsoid->type == OBJ3D_ELLIPSOID)
					elipsoid->dataterm = ellipsoid_likelyhood_sn(img,width,height,depth,
							elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
							elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
				else if(elipsoid->type == OBJ3D_CYLINDER)
					elipsoid->dataterm = ellipsoid_likelyhood_cylinder(img,width,height,depth,
							elipsoid->y,elipsoid->x,elipsoid->z,elipsoid->a,elipsoid->b,elipsoid->c,
							elipsoid->alpha,elipsoid->beta,elipsoid->gamma,elipsoid->n,&(elipsoid->dist),mpp.grad_dir,thredh,config);
				if(elipsoid->dataterm<l_th){
					elipsoid->gcut_weight = (int)(500.*(1+elipsoid->dataterm));
					double dr = max(elipsoid->a,elipsoid->b);
					double dx = 2.*(fabs((double)elipsoid->c*sin(elipsoid->beta))+dr)+1.;
					double dy = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*sin(elipsoid->alpha))+dr)+1.;
					double dz = 2.*(fabs((double)elipsoid->c*cos(elipsoid->beta)*cos(elipsoid->alpha))+dr)+1.;
					if(elipsoid->type == OBJ3D_ELLIPSOID){
						drawSuperEllipsoid33(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
							elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						overlap = elipsoid_overlap_check_all(elipsoid,newlink,newlink_num);
					}
					else if(elipsoid->type == OBJ3D_CYLINDER){
						drawSuperEllipsoid_Cylinder(elipsoid, (int)dx, (int)dy, (int)dz, elipsoid->a, elipsoid->b, elipsoid->c,
							elipsoid->alpha, elipsoid->beta, elipsoid->gamma, elipsoid->n,1,config);
						overlap = elipsoid_overlap_check_all(elipsoid,newlink,newlink_num);
					}
					else
						overlap = elipsoid_overlap_check_all_disk(elipsoid,newlink,newlink_num);
					if(!overlap){
						newlink_num++;
						LinkedList3DInsert(newlink, newlink_num, elipsoid); 
					}
					else{
						free_volume((void***)elipsoid->ptmpmask);
						elipsoid->ptmpmask = NULL;
						free(elipsoid);
					}
				}
				else{
					free(elipsoid);
				}
			}
			printf("k = %d, new_obj = %d\n",k,newlink_num);
		}

		// deadth step
		if((iter!=1)&&(newlink_num)){
//			printf("Graph cut start\n");
			graph_cut33(config,newlink,newlink_num);
//			printf("Graph cut end\n");

			pobj = config->elink->next;
			for(int i = 1;i<=config->mp_num;i++){
				if (pobj->index->my_type == 0){ // delete from exist config
					pobj = pobj->next;
					Qdelete_elipsoid3(config, i);
					i--;
				}
				else							// keep in config
					pobj = pobj->next;
			}
			pobj = newlink->next;
			for(int i = 1;i<=newlink_num;i++){
				if (pobj->index->my_type == 0){ // delete from newlink
					pobj = pobj->next;
					Qdelete_elipsoid3(newlink, i, &newlink_num);
					i--;
				}
				else{							// move from newlink to config
					config->mp_num++;
					LinkedList3DInsert(config->elink,config->mp_num, pobj->index); 
					pobj = pobj->next;
					LinkedList3DDelete(newlink, i);//do not delete object itself
					newlink_num--;
					i--;
				}
			}
		}
		printf("obj num = %d\n",config->mp_num);

		clock_t mid_time2=clock();
		cout<< "Iteration running time is: "<<static_cast<double>(mid_time2-mid_time1)/CLOCKS_PER_SEC<<"s"<<endl;
		cout<< "Total running time is    : "<<static_cast<double>(mid_time2-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
		pobj = config->elink;
		int exe_time = static_cast<int>(mid_time2-start_time); //ms
		config->exe_time[iter] = exe_time;
		config->num_obj[iter] = config->mp_num;
		config->energy[iter] = 0;
		for (int k = 1;k<=config->mp_num;k++){
			pobj = pobj->next;
//			prior = C_prior_elipsoid3(config, k, pen);
//			dtmp = pobj->index->dataterm
//				+pobj->index->inter_e*mpp.intra_pen/2.;
			config->energy[iter] += pobj->index->dataterm;
		}
		printf("E = %1.5f\n",config->energy[iter]);
	}
	clock_t end_time=clock();
	double running_time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
	cout<< "Running time is: "<< running_time<<"s"<<endl;

	int z;
	pobj = config->elink;
	cfg->num_obj=0;
	cfg->total_E = 0;
	for (int k = 1;k<=config->mp_num;k++){
		pobj = pobj->next;
		z = pobj->index->z;
//		prior = C_prior_elipsoid3(config, k, pen);
		cfg->ellipsoid[cfg->num_obj].center.x = pobj->index->x;
		cfg->ellipsoid[cfg->num_obj].center.y = pobj->index->y;
		cfg->ellipsoid[cfg->num_obj].center.z = pobj->index->z;
		cfg->ellipsoid[cfg->num_obj].a = pobj->index->a;
		cfg->ellipsoid[cfg->num_obj].b = pobj->index->b;
		cfg->ellipsoid[cfg->num_obj].c = pobj->index->c;
		cfg->ellipsoid[cfg->num_obj].alpha = pobj->index->alpha;
		cfg->ellipsoid[cfg->num_obj].beta = pobj->index->beta;
		cfg->ellipsoid[cfg->num_obj].gamma = pobj->index->gamma;
		cfg->ellipsoid[cfg->num_obj].single_E = pobj->index->dataterm;
		cfg->ellipsoid[cfg->num_obj].e0 = pobj->index->dist;
		cfg->ellipsoid[cfg->num_obj].multiple_E = 0;
		cfg->ellipsoid[cfg->num_obj].obj_id = pobj->index->obj_id;
		cfg->ellipsoid[cfg->num_obj].n = pobj->index->n; // ELLIPSE_POWER
		cfg->total_E += pobj->index->dataterm;
		cfg->ellipsoid[cfg->num_obj].alpha_mean = pobj->index->alpha_mean;
		cfg->ellipsoid[cfg->num_obj].beta_mean = pobj->index->beta_mean;
		cfg->ellipsoid[cfg->num_obj].fit_error = pobj->index->fit_error;
		cfg->num_obj++;
	}
	printf("obj num = %d\n",cfg->num_obj);
	free_img((void **)r_p_M);
	free(tan_phi);

	return running_time;
}