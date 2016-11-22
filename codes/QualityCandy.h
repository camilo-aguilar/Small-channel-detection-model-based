#ifndef _QUALITYCANDY_H_
#define _QUALITYCANDY_H_

#include <stdlib.h>
#include <stdio.h>
#include "randlib.h"
#include "math.h"
#include "util.h" //KDW
#include "ndmulti_e.h" //KDW

#define QUALITY_CANDY
#define EN_NEW_SEG_AT_NON_FREE_END
#define MAX_CONNECTION_NUM  30
#define MAX_NEAR_NUM		                30
#define _PI					3.141592654


#define L_MAX			10 //11KDW 30 		//10
#define L_MIN			6 //KDW 10.0//8.0
#define W_MIN			6 //KDW 4
#define W_MAX			13 //KDW 9

#define THETA_MIN			 (0.0)
#define THETA_MAX			(_PI)
#define DELTA_MAX		0//0.4 //KDW 0.0
#define PARRALLEL_RANGE  -1.0 //_PI/40
#define TAU_MAX			_PI/4.0
#define TAU_MAX_FREESEG	 _PI//KDW _PI/6.0
#define DATATERM		exp(double(2))
#define NEIGHBOORHOOD	2//5
#define SEARCH_R 10

#define I_WIDTH width
#define I_HEIGHT height

#define THIN_W	4
#define THICK_W	6
//#define THIN_W	6
//#define THICK_W	4
#define PEN_1 2.6
#define PEN_2 1.8
#define PEN_3 3.5

#define RADIUS_SEGS 30.0

#define SMALL_BETA	1.5//0.5
#define BIG_BETA	1.5
#define SEG_OFFSET  0.7
#define MAX_K    20

#define SIGN(a) ((a)>=0?1:(-1))
#define ABS(A) ((A)>0?(A):-(A))
#ifndef MAX
	#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
	#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#define DIST(x,y,xx,yy) ((x-xx)*(x-xx)+(y-yy)*(y-yy))
#define sigma(tau_ij2, tau_m2)  (((1+(tau_m2))/(1+(tau_ij2))-1)/(tau_m2))

typedef struct 
{
	int x;
	int y;
}site;

typedef struct lineObj
{
	int x;
	int y;
	int len;
	int width;
	double theta;
	site enda; // enda_x < endb_x, if enda_x = endb_x, then enda_y < endb_y 
	site endb;

	int img_num; // img_num where object theta belongs

	int type;  /*0: free segment 1: single segment  2: double segment*/
	int enda_L_Num;
	int endb_L_Num;
	struct lineObj *enda_L[MAX_NEAR_NUM];   /*objects near to enda within the distance of NEIGHBOORHOOD*/
	struct lineObj *endb_L[MAX_NEAR_NUM];	/*objects near to endb within the distance of NEIGHBOORHOOD*/

	int enda_C_Num;
	int endb_C_Num;
	struct lineObj *enda_C[MAX_CONNECTION_NUM]; /*objects connecting to enda*/
	struct lineObj *endb_C[MAX_CONNECTION_NUM];	/*objects connecting to endb*/

	double dataterm;
	double engergy_for_transition;

}lineObj;

typedef struct Node
{
	lineObj *index;
	struct Node *next;
}Node,*LinkedList;

typedef struct NClinks
{
	lineObj *object1;
	lineObj *object2;
	struct NClinks *next;

}NClinks,*Neighbor_links, *Connection_links;

typedef struct 
{
	int n_f;    /*number of free segments*/
	int n_s;    /*number of single segments*/
	int n_d;    /*number of double segments*/

	double lambda;  //poison intensity

	double p_b_d;
	double p_t;
	double p_c;  //connect free ends

	double p_c_FtoC;
	double p_c_CtoF;

	double p_t_f;
	double p_t_s;
	double p_t_d;

	double p_f;
	double p_s;
	double p_d;
	

	double p_f_d;
	double p_f_b;
	double p_s_d;
	double p_s_b;
	double p_d_d;
	double p_d_b;

	double w_f;
	double w_s;
	double w_d;

	double w_io;
	double w_eo;
	double beta;
	double gamma_d; //KDW

	LinkedList link_f;   /*list of free segments*/
	LinkedList link_s;   /*list of single segments*/
	LinkedList link_d;   /*list of double segments*/

 
	int neighbor_n;
	int connection_n;
	Neighbor_links neighbor_l;		/*list of neighbors*/
	Connection_links connection_l;	/*list of connections*/

	double energy;
}Candy;

Candy* CandyInit(MPP_Parameters mpp);
void freeCandy(Candy *C);
void DrawLine(double **img,Candy *C);
int SaveCandy2MP(Candy *C, NeckDent *mp);
//void Candy_Model(double **input_img,double ***img_mpp_l,double ***img_seg_l,double **output_img,int height,int width,double T,int iter,double **patch, int patch_len);
int Candy_Model(double **input_img, double **lm, double ****img_mpp_l, double ****img_seg_l, double **output_img, int **output_seg, int height,int width,double T_origin,int iter,double **patch, int patch_len,double *mu, double *cov, int label,
				 NeckDent *mp, MPP_Parameters mpp); 

LinkedList LinkedListInit();
void LinkedListInsert(LinkedList L,int i, lineObj* x);
void LinkedListDelete(LinkedList L,lineObj* x);
lineObj* ReturnObj (LinkedList L, int i);

//void AddFreeSeg (Candy *M, double **img, double ***mid_img,int img_height, int img_width, double T,int test,double **patch, int patch_len);
//void AddSingleSeg (Candy *M, double **img,double ***mid_img, int img_height, int img_width, double T,double **patch, int patch_len);
//void AdddoubleSeg (Candy *M, double **img, double ***mid_img, int img_height, int img_width, double T,double **patch, int patch_len, int test);

void AddFreeSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,int test,double **patch, int patch_len, double ***Matrix, double prior_pen1, double prior_pen2,int **dilated_img);
void AddSingleSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len, double ***Matrix, double prior_pen1, double prior_pen2);
void AdddoubleSeg (Candy *M, double **img, double **lm, int **img_seg, double ****img_mpp_l, double ****img_seg_l,int img_height, int img_width, double T,double **patch, int patch_len , double ***Matrix, double prior_pen1, double prior_pen2, int test);
void KillFreeSeg(Candy *M, double **img, int img_height, int img_width, double T);
void KillSingleSeg(Candy *M, double **img, int img_height, int img_width, double T);
void KillDoubleSeg(Candy *M, double **img, int img_height, int img_width, double T);

//int FreeSeg_length_move(Candy *M, double **img, double ***mid_img, int img_height, int img_width, double T,double **patch, int patch_len);
//int FreeSeg_theta_move(Candy *M, double **img, double ***mid_img, int img_height, int img_width, double T,double **patch, int patch_len);
//int SingleSeg_freeEnd_move(Candy *M, double **img, double ***mid_img, int img_height, int img_width, double T, double **patch, int patch_len);

int FreeSeg_length_move(Candy *M, double **img, int **seg_img, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);
int FreeSeg_theta_move(Candy *M, double **img, int **seg_img, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);

//KDW
int FreeSeg_width_move(Candy *M, double **img, int **seg_img, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T,double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);
int FreeSeg_freeEnd_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);
int FreeSeg_center_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);

int SingleSeg_freeEnd_move(Candy *M, double **img, int **seg_img, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);
int SingleDoubleSeg_Connection_move(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, int img_height, int img_width, double T, double **patch, int patch_len,double ***Matrix, double prior_pen1, double prior_pen2);

double theta_from_two_ends(site enda, site endb);

int Bad_EO_freeSeg(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc);
int Bad_EO(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc);
int Bad_EO_from_Link(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, double *g_Rc);
int Bad_EO_freeSeg_death(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc);
int Bad_EO_death(int neighborR,double searchRatio, lineObj *line, Candy *M, double *g_Rc);
int Bad_EO_from_Link_death(lineObj *line,int neighboorR,double searchRatio, LinkedList link, int n, double *g_Rc);
int Bad_EO_transition(int neighborR,double searchRatio, lineObj *line, Candy *M, int choose_num, int end_type, double *g_Rc);
int Bad_EO_connection_move(int neighborR,double searchRatio, lineObj *line, Candy *M, lineObj **obj, int obj_num, double *g_Rc);
int Bad_EO_objects(int neighboorR,double searchRatio, lineObj **obj, int obj_num, double *g_Rc);
int Bad_EO_freeEndsC(int neighborR,double searchRatio, lineObj *line, Candy *M, lineObj *obj1, double *g_Rc);

int Bad_IO(lineObj *Seg, Candy *M, double *g_Rio);
int Bad_IO_from_Link(lineObj *Seg, LinkedList link, int n, double *g_Rio);
int Bad_IO_transition(lineObj *Seg, Candy *M, int choose_num, int end_type, double *g_Rio);
int Bad_IO_connection_move(lineObj *Seg, Candy *M, lineObj **obj, int obj_num, double *g_Rio);
int Bad_IO_objects(lineObj **obj, int obj_num, double *g_Rio);
int Bad_IO_freeEndsC(lineObj *Seg, lineObj *obj1, Candy *M, double *g_Rio);


site GenerateEndb(site enda,int len,double theta,int img_height,int img_width);
site SelectEndfromSingleLink(LinkedList link,int num, int choose_num);
site SelectEndfromFreeLink(LinkedList link,int num, int choose_num);

int DistSquare(site a, site b);
//int CheckConnection(site a, site b);
double CheckConnection(site a, site b, site aa, site bb, double searchRatio, int neighboorR,int *if_eo, int *if_connect, int *which_end_is_neighboor);
void UpdateItsNeighboors_born (lineObj* seg, Candy* M);
void updateNeighboorLink_born(Candy *M, lineObj* target,lineObj* seg,site end);
void UpdateItsNeighboors_death (lineObj* seg, Candy* M);
void updateNeighboorLink_death(Candy *M, lineObj* target,lineObj* seg,site end);

void UpdateItsNeighboors_born_freeSeg (lineObj* seg, Candy* M);
void UpdateItsNeighboors_born_SingleSeg (lineObj* seg, Candy* M);
void UpdateItsNeighboors_born_DoubleSeg (lineObj* seg, Candy* M);

void UpdateItsNeighboors_death_freeSeg (lineObj* seg, Candy* M);
void UpdateItsNeighboors_death_SingleSeg (lineObj* seg, Candy* M);
void UpdateItsNeighboors_death_DoubleSeg (lineObj* seg, Candy* M);


void line(int x1, int y1, int x2, int y2, int color, double **img);
bool doIntersect(site p1, site q1, site p2, site q2);

#if 0
void AddSingleSeg_allends (Candy *M, double **img, int img_height, int img_width, double T);
void KillSingleSeg_allends(Candy *M, double **img, int img_height, int img_width, double T);
void AdddoubleSeg_allends (Candy *M, double **img, int img_height, int img_width, double T);
void KillDoubleSeg_allends(Candy *M, double **img, int img_height, int img_width, double T);
#endif
double dataterm(double **input_img,site enda, site endb);
double dataterm_double(double **input_img,site enda, site endb);
//double dataterm(double **input_img,site enda, site endb, double theta,int len,int theta_num,double **patch, int patch_len);
void likelyhoodMap (double **input_img, double **out_img_mpp, double **out_img_seg, double *mu, double *cov,double **patch, int L, int upper, int lower, double theta, int height, int width,int patch_len,int flag);
double GetLikelyhood(int y, int x, double **input_img, int L, double theta, int upper_w, int lower_w,double **patch,int patch_len, int height, int width);
double GetLikelyhood_seg(int y, int x, double **input_img, int L, double theta, int upper_w, int lower_w,  double **patch,int patch_len,double *mu, double *cov, int **indexmask,int flag);
site drawRect(double **patch, int patch_len, int y, int x, int L, double theta, int upper_w, int lower_w);
void subfill (double **tmpmask, int x, int y,int new_c );

double Echange_from_neighboors_born(Candy *M, lineObj *obj);
double Echange_from_neighboors_death(Candy *M, lineObj *obj);

NClinks *NClinksInit();
void AddNClinks(NClinks *l, int index, lineObj *obj1, lineObj *obj2);
void KillNClinks(NClinks *l, lineObj *obj1, lineObj *obj2);
//int Conncet_freeEnds(Candy *M, double **img, double ***mid_img, double **patch, int patch_len,int img_height, int img_width, double T);
//int Seperate_connectedEnds(Candy *M, double **img, double ***mid_img, double **patch, int patch_len, int img_height, int img_width, double T);
int Seperate_connectedEnds(Candy *M, double **img, int **img_seg, double ****img_mpp_l, double ****img_seg_l, double **patch, int patch_len,int img_height, int img_width, double ***Matrix, double prior_pen1, double prior_pen2 ,double T);
int Conncet_freeEnds(Candy *M, double **img, int **img_seg, double ****img_mpp_l,double ****img_seg_l, double **patch, int patch_len,int img_height, int img_width, double ***Matrix, double prior_pen1, double prior_pen2, double T);
void NC_pairs(lineObj *obj1, lineObj *obj2,site *end1, site *end2);
void killlineSeg(lineObj *obj, Candy *M);
void AddlineSeg(lineObj *obj, Candy *M);

void K_means(double **img1, double *mu,double *cov, int height,int width,int k);
//double dataterm_rec(lineObj *line, double **img,double **patch, int patch_len, int line_w);
//double dataterm_rec(lineObj *line, double **img_mpp_l,double **img_seg_l,double **patch, int patch_len, int line_w);
double dataterm_rec(int **img_seg, lineObj *line, double **patch, int patch_len, int line_w,double ***Matrix, double prior_pen1, double prior_pen2, int height,int width);

//KDW DrawCandyLine
//void DrawCandyLine(IplImage *image, lineObj *mp, int i, CvScalar color, int text_type, int line_thickness);
int count_all_ends(Candy *M, site *end);
int count_all_ends(Candy *M, site *end, lineObj **obj);
int SelectNearestEnd(site a,lineObj *belong,Candy *M,int searchR, site *candid);

#endif /* _QUALITYCANDY_H_ */
