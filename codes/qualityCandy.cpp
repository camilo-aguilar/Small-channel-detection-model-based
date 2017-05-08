#include "QualityCandy.h"
#include "tiff.h"
#include "allocate.h"
#include <iostream> //KDW
#include <time.h>



#define DEBUG_CAMILO 0
#define INCLUDE_OPENCV_END 0

#if INCLUDE_OPENCV
	#include "gui_functions.h"
#endif

#if INCLUDE_OPENCV_END
	#include "gui_functions.h"
#endif


using namespace std;




void drawSegfromfile(Candy *C)
{
	FILE *fp;

 if ( ( fp = fopen ( "C_output99000000.txt", "rb" ) ) == NULL ) {
    fprintf ( stderr, "cannot open file MBD_img.tif\n");
    exit ( 1 );
  }
		int enda_x;
		int enda_y;
		int endb_x;
		int endb_y;
		double theta;
		int len;
		int flag;


	while (!feof(fp))
	{
		fscanf(fp,"%d,%d,%d,%d,%d,%lf,%d\n",&flag,&enda_x,&enda_y,&endb_x,&endb_y,&theta,&len);

		if (flag != 4)
		{
			lineObj *freeSeg = (lineObj *)malloc(sizeof(lineObj));

			freeSeg->width = W_MIN;
			freeSeg->len = len;
	
			freeSeg->theta = theta;
	
			freeSeg->enda.x  = enda_x;
			freeSeg->enda.y = enda_y;
			freeSeg->endb.x  = endb_x;
			freeSeg->endb.y = endb_y;
						
			C->n_f++;
			freeSeg->type = 0;
		    LinkedListInsert( C->link_f,C->n_f, freeSeg); 
		}

	}
}

/* Draws a line between (x1,y1) and (x2,y2) in img[y][x] = color */
void line_back(int x1, int y1, int x2, int y2, int color, double **img)
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
	img[y1][x1] += color;
	if(e >= 0) {
	y1 += iny;
		e-= dx;
	}
	e += dy; x1 += inx;
	}
	} else {
	dx <<= 1;
	e = dx - dy;
	dy <<= 1;
	while (y1 != y2) {
	img[y1][x1] += color;
	if(e >= 0) {
	x1 += inx;
	e -= dy;
	}
	e += dx; y1 += iny;
	}
	}
	img[y1][x1] += color;
}


/*Draws all the lines in the Candy so far*/
void DrawLine_back(double **img,Candy *C)
{
	int n_f = C->n_f;
	int n_s = C->n_s;
	int n_d = C->n_d;

	site enda, endb;

	Node *p = C->link_f;

	for (int i = 0;i<n_f;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;

		line_back(enda.x,enda.y,endb.x,endb.y,160,img);
	}

	p = C->link_s;
	for (int i = 0;i<n_s;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;
		line_back(enda.x,enda.y,endb.x,endb.y,160,img);
	}

	p = C->link_d;
	for (int i = 0;i<n_d;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;
		line_back(enda.x,enda.y,endb.x,endb.y,160,img);
	}
}



Candy* CandyInit(MPP_Parameters mpp)
{
	double Constant_For_Normalizing_Lambda = (double)((L_MAX-L_MIN)*(W_MAX-W_MIN)*(THETA_MAX-THETA_MIN));
	Candy *C;
	C = (Candy *)malloc(sizeof(Candy));

	C->beta = 1;					//Beta
	C->lambda = 1.0/Constant_For_Normalizing_Lambda;					//Lambda



	C->link_d = LinkedListInit();	//Double Seg Links
	C->link_f = LinkedListInit();	//Free Seg Links
	C->link_s = LinkedListInit();	//Single Seg Links

	C->connection_l = NClinksInit();//
	C->neighbor_l = NClinksInit();	//
	
	C->connection_n = 0;			//
	C->neighbor_n = 0;				//
	
	C->n_d = 0;						//Number of Double Segs
	C->n_f = 0;						//Number of Free Segments
	C->n_s = 0;						//Number of Single Segments

	C->p_f_b = P_BIRTH_F; 				/* Free Segment Birth */
	C->p_f_d = 1-C->p_f_b; 				/* Free Segment Death */
	
	C->p_s_b = P_BIRTH_S; 				/* Single Segment Birth */
	C->p_s_d = 1 - C->p_s_b; 				/* Single Segment Death */
	
	C->p_d_b = P_BIRTH_D; 				/* Double Segment Birth */
	C->p_d_d = 1- C->p_d_b; 				/* Double Segment Death */


	C->gamma_d = mpp.gamma_d; 		//
	
	C->w_eo = mpp.w_eo; 			/* Weight External Orientation */
	C->w_f = mpp.w_f;   			/* Weight free Seg */
	C->w_s = mpp.w_s;   			/* Weight Single Seg */
	C->w_d = mpp.w_d;   			/* Weight Double Seg */
	C->w_io = mpp.w_io; 			/* Weight Internal Bad Orientation */

	C->p_b_d =  P_BIRTH_DEATH_STEP;  				/* Birth and Death Step */
	C->p_t =    P_TRANSLATION_STEP;     				/* Translation Step */
	C->p_c =    P_CONNECTION_STEP;     				/* Connection Step */

	C->p_c_CtoF = CONNETECTED_TO_FREE; 				/* Connected to Free */
	C->p_c_FtoC = FREE_TO_CONNECTED; 				/* Free to Connected */

	C->p_t_f = TRANSITION_FREE_SEGMENT; 				/* Transition Free Seg*/
	C->p_t_s = TRANSITION_SINGLE_SEGMENT; 				/* Transition Single Seg*/
	C->p_t_d = TRANSITION_DOUBLE_SEGMENT; 				/* Transition Double Seg*/
 
	C->p_f = P_PICK_F; 					/* Free (B&D) */
	C->p_s = P_PICK_S ;  					/* Single */
	C->p_d = P_PICK_D ;  					/*Double */

	C->energy = 0;
	C->Vo = 0;
	C->VRio = 0;
	C->VReo = 0;

	return C;
}


void freeCandy(Candy *C)
{
	free(C->link_d);
	free(C->link_f);
	free(C->link_s);
	free(C);
}



int DrawEach_Seg(int **input_img,Node *obj, int height, int width, double *mu, double *cov, double **patch, int patch_len)
{
	site enda, endb, offset;
	int lower_w = 2;
	int upper_w = obj->index->width;

	for (int i =0;i<patch_len;i++)
	for (int j = 0;j<patch_len;j++)
		patch[i][j] = 0;

	enda = obj->index->enda;
	endb = obj->index->endb;
	double theta = obj->index->theta;
	int len = obj->index->len;
	site start_end;

	int xx = enda.x + (int)((double)len*cos(theta));
	int yy = enda.y - (int)((double)len*sin(theta));

	if(DIST(xx,yy,endb.x,endb.y)< 3*3)
		start_end = enda;
	else
	{
		xx = endb.x + (int)((double)len*cos(theta));
		yy = endb.y - (int)((double)len*sin(theta));
		if(DIST(xx,yy,enda.x,enda.y)< 3*3)
			start_end = endb;
		else
			{
				printf("find start_end wrong\n");
				return 0;
		}
	}

	offset = drawRect(patch,patch_len,start_end.y,start_end.x,len,theta,upper_w,lower_w);


	for (int i = 0;i < patch_len;i++)
		for (int j = 0; j< patch_len;j++)
	{
		if (patch[i][j]== 128)
		{ 
			if ((i + offset.y > 0 )&& (i + offset.y < height) && (j+offset.x > 0) &&(j+offset.x< width))
			{
				int yyy =i+offset.y;
				int xxx = j+offset.x;

				input_img[yyy][xxx] = 1;

	
			}
		}
	}
	return 1;  
}


void DrawLine_Seg(int **tmp_img, Candy *C, int height,int width,double *mu,double *cov,double **patch,int patch_len)
{
//KDW	site offset;
		 
	int n_f = C->n_f;
	int n_s = C->n_s;
	int n_d = C->n_d;
	 

	Node *p = C->link_f;

	for (int i = 0;i<n_f;i++)
	{
		p= p->next; 
		DrawEach_Seg(tmp_img,p,height,width,mu,cov,patch,patch_len);
	}

	p = C->link_s;
	for (int i = 0;i<n_s;i++)
	{
		p= p->next;
		DrawEach_Seg(tmp_img,p,height,width,mu,cov,patch,patch_len);
	}

	p = C->link_d;
	for (int i = 0;i<n_d;i++)
	{
		p= p->next;
		DrawEach_Seg(tmp_img,p,height,width,mu,cov,patch,patch_len);
	}
}



void writeC(Candy *C,int i)
{
	FILE *fp;
	
	char a[100];
	sprintf(a,"C_output%d.txt",i);

	 if ( ( fp = fopen ( a, "wb" ) ) == NULL ) {
    exit ( 1 );
  }

	int n_f = C->n_f;
	int n_s = C->n_s;
	int n_d = C->n_d;
	 

	Node *p = C->link_f;

	for (int i = 0;i<n_f;i++)
	{
		p= p->next; 
		int enda_x = p->index->enda.x;
		int enda_y = p->index->enda.y;
		int endb_x = p->index->endb.x;
		int endb_y = p->index->endb.y;
		double theta = p->index->theta;
		int len = p->index->len;

		 fprintf(fp,"0,%d,%d,%d,%d,%f,%d\n",enda_x,enda_y,endb_x,endb_y,theta,len);
	}

	p = C->link_s;
	for (int i = 0;i<n_s;i++)
	{
		p= p->next;
		int enda_x = p->index->enda.x;
		int enda_y = p->index->enda.y;
		int endb_x = p->index->endb.x;
		int endb_y = p->index->endb.y;
		double theta = p->index->theta;
		int len = p->index->len;

		 fprintf(fp,"1,%d,%d,%d,%d,%f,%d\n",enda_x,enda_y,endb_x,endb_y,theta,len);
	}

	p = C->link_d;
	for (int i = 0;i<n_d;i++)
	{
		p= p->next;
		int enda_x = p->index->enda.x;
		int enda_y = p->index->enda.y;
		int endb_x = p->index->endb.x;
		int endb_y = p->index->endb.y;
		double theta = p->index->theta;
		int len = p->index->len;

		 fprintf(fp,"2,%d,%d,%d,%d,%f,%d\n",enda_x,enda_y,endb_x,endb_y,theta,len);
	}
	fclose(fp);
}



void dilation(int **seg_img,int **dilated_img,int width,int height)
{
	for(int i = 1;i<height-1;i++)
		for(int j = 1;j<width-1;j++)
		{
			int flag = 0;
			for(int k = i-1;k<=i+1;k++)
				for(int m = j-1;m<=j+1;m++)
				{
					if(seg_img[k][m] == 2)
					{
						flag = 1;
						break;
					}
				}
			if(flag == 1)
				dilated_img[i][j] = 1;
			else
				dilated_img[i][j] = 0;
		}
}





int Candy_Model(double **input_img, double **lm, double ****img_mpp_l, double ****img_seg_l, double **output_img, int **output_seg, int height,int width,double T_origin,int iter,double **patch, int patch_len,double *mu, double *cov, int label,
				 NeckDent **mp, MPP_Parameters mpp) 
{

	int mp_num = 0;
	double Echange = 0;

	FILE *energy_log_file;

 	if ( ( energy_log_file = fopen ( "Energy.txt", "wb") ) == NULL ) {
    fprintf ( stderr, "cannot open file Energy.txt\n");
    exit ( 1 );
  	}

	struct TIFF_img inter_img_mrf;

	get_TIFF ( &inter_img_mrf, height, width, 'g' );

	Candy *C;
	C = CandyInit(mpp);

	double T = 0;

	int **index_matrix = (int **)get_img(2,patch_len*patch_len,sizeof(int));
	int **indexmask=(int **)get_img(patch_len,patch_len,sizeof(int));

	int **tmp_img = (int **)get_img(width,height,sizeof(int));
	int **dilated_img = (int **)get_img(width,height,sizeof(int));

	for(int ii = 0; ii<height;ii++)
		for(int jj = 0; jj<width;jj++)
			{
				tmp_img[ii][jj] = 0;
				dilated_img[ii][jj] = 0;
			}	

	double prior_pen1 = mpp.error_th;//KDW
	double prior_pen2 = PEN_2;//0.5;


	double **Matrix[MAX_K];
	
	for(int i = 0; i<label;i++)
	{
		Matrix[i] = (double **)get_img(width,height,sizeof(double));
	}


	for(int k = 0; k<label;k++)
		for (int i = 0;i<height;i++)
			for(int j = 0; j< width;j++)
			{
				Matrix[k][i][j] = (input_img[i][j]-mu[k])*(input_img[i][j]-mu[k])/(2*cov[k])+log(sqrt(2*_PI*cov[k]));
			}

	for(int ii = 0; ii<height;ii++)
		for(int jj = 0; jj<width;jj++)
		{
			output_seg[ii][jj] = (int)input_img[ii][jj];
		}


	for (int ii = 0; ii < height; ii++ )
	{
		for (int jj = 0; jj < width; jj++ ) 
		{
			int pixel = (int32_t)(output_seg[ii][jj]);				
			if(pixel>0) 
			{
				
				inter_img_mrf.mono[ii][jj] = MIN(pixel*256/3,255);
			}
			else 
			{
				inter_img_mrf.mono[ii][jj] = 0;
			}
		}
		
	 }
	
	clock_t start_time=clock();

	T = mpp.T0;
	printf("Start RJMCMC\n");
	for (int i = 0;i<mpp.iter_num;i++)
	{
		double E1 = C->energy;
		double E2 = -1* ((C->Vo*C->gamma_d) + (C->VReo*C->w_eo) + (C->VRio*C->w_io) + (C->n_s*C->w_s) + (C->n_f*C->w_f) + (C->n_d * C->w_d));
		#if INCLUDE_OPENCV
				if(Echange)//i%(mpp.iter_num/NUM_WINDOWS)==0)
				{						
					
					printf("\n1:%.2f\n", E1);
					printf("2:%.2f\n", E2);
					
					//printf("2:%.2f\n", (-C->Vo*C->gamma_d) - ((double)C->n_f*C->w_f));
					display_image_double(input_img, height, width, C);
					
				}
		#else
				if(i%(mpp.iter_num/NUM_WINDOWS)==0)
				{						

				}
		#endif

		if(fabs(E1 - E2) > 0.0001)
		{
			printf("Error in Energies");
			printf("\n1:%.2f\n", E1);
			printf("2:%.2f\n", E2);
			exit(1);
		}
		Echange = 0;
		if(i > 0 && i % (mpp.iter_num/10) == 0) 
		{
			printf(" %2.2f Percent of Iterations Completed \n", (float)i/(float)mpp.iter_num * 100.00);
		}

		double r = random2();
		
		
		if((C->n_f + C->n_s + C->n_d)==0)
		{
			r = 0;
		}

		if (r < C->p_b_d)   //birth and death step
		{
			double rr = random2();
			double rrr = random2();

				
			if (rr < C->p_f)  // free segment
			{


				if (rrr < C->p_f_b)
				{
					Echange = AddFreeSeg(C,input_img,lm,output_seg,height,width,T,patch,patch_len, Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}
				else
				{

					Echange = KillFreeSeg(C,input_img,height,width,T);
					C->energy = C->energy + Echange; 

				}
			}
			else if (rr < C->p_s+C->p_f)  //single segment
			{
				if (rrr < C->p_s_b)
				{
					Echange =AddSingleSeg(C,input_img,lm,output_seg, img_mpp_l,img_seg_l,height,width,T,patch,patch_len, Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}
				else
				{
					Echange = KillSingleSeg(C,input_img,height,width,T);
					C->energy = C->energy + Echange;
				}
			}
			else
			{
				if (rrr < C->p_d_b)
				{
					Echange = AdddoubleSeg(C,input_img,lm,output_seg,img_mpp_l,img_seg_l,height,width,T,patch,patch_len, Matrix, prior_pen1, prior_pen2,0);
					C->energy = C->energy + Echange;
				}
				else
				{

					Echange = KillDoubleSeg(C,input_img,height,width,T);
					C->energy = C->energy + Echange;
				}
			}
			/* End Birth/Death Kernel */
		}
		else if (r <C->p_b_d +C->p_t)
		{
			/* Transition Kernel */
			
			double rr = random2();
			double rrr = random2();
			if (rr < C->p_t_f)
			{
				if (rrr < 0.2)
				{
					Echange = FreeSeg_length_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height,width,T,patch,patch_len,Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}
				else if (rrr< 0.4)
				{
					Echange = FreeSeg_theta_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height,width,T,patch,patch_len,Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}
				else  if (rrr< 0.6) //KDW 
				{
					Echange = FreeSeg_width_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height,width,T,patch,patch_len,Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}
				else  if (rrr< 0.8)  //KDW
				{
					Echange = FreeSeg_freeEnd_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}
				else //KDW
				{
					Echange = FreeSeg_center_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);
					C->energy = C->energy + Echange;
				}

				/* End of Transition Kernel */

			}
			else if (rr < C->p_t_s+C->p_t_f)
			{
				Echange = SingleSeg_freeEnd_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);
				C->energy = C->energy + Echange;

			}
			else  // no disconnect & move
			{
				Echange = SingleDoubleSeg_Connection_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);
				C->energy = C->energy + Echange;
			}
		}
		else   //connetion and separation step
		{

			double rr = random2();
			if (rr < C->p_c_FtoC)
			{
				Conncet_freeEnds(C,input_img,output_seg, img_mpp_l,img_seg_l,patch,patch_len,height,width,Matrix, prior_pen1,prior_pen2,1.0/T);
			}
			else
			{
				Seperate_connectedEnds(C,input_img,output_seg, img_mpp_l,img_seg_l,patch,patch_len,height,width,Matrix, prior_pen1,prior_pen2,1.0/T);
			}

		}

		T = T*(mpp.de_coeff); 
		if(T<0.00000000001)
		{
			printf("Warning, Temperature Reached 0\n");
			break;
		}

		if(i%(mpp.iter_num/NUM_WINDOWS)==0)
		{
			//fprintf(energy_log_file,"%.5f\n",C->energy);
		}	

	}
	fprintf(energy_log_file,"%.5f\r\n %.5f\r\n %.5f\r\n %d\r\n %d\r\n %d\r\n",C->Vo, C->VRio, C->VReo, C->n_f,C->n_s,C->n_s);
	fclose(energy_log_file);

	#if INCLUDE_OPENCV_END
	display_image_double(input_img, height, width, C);
	#endif


	clock_t end_time=clock();
	cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<endl;
	
	printf("free=%d,single=%d,double=%d\n",C->n_f,C->n_s,C->n_d);
	DrawLine(output_img,C);
	//int mp_num = SaveCandy2MP(C, &mp);
	*mp = SaveCandy2MP(C, &mp_num);
	freeCandy(C);
	free_img( (void**)index_matrix );
	free_img( (void**)indexmask );
	free_img((void**)dilated_img);
	printf("MP_NUM = %d ", mp_num);
	printf("\n Finished Quality Candy \n \n");
	return mp_num;

}





site SelectEndfromFreeLink(LinkedList link,int num, int choose_num)
{
	Node *p = link;
	for (int i = 0;i< choose_num;i++)
		p = p->next;
	if (random2()<0.5)
		return ((p->index)->enda);
	else
		return ((p->index)->endb);
}


site SelectEndfromSingleLink(LinkedList link,int num, int choose_num)
{
	site no_end;
	no_end.x =0;
	no_end.y = 0;
	Node *p = link;
	for (int i = 0;i< choose_num;i++)
		p = p->next;
	if ((p->index)->enda_C_Num == 0)
		return ((p->index)->enda);
	else if ((p->index)->endb_C_Num == 0)
		return ((p->index)->endb);
	else
		return(no_end);
}



void UpdateItsNeighboors_born (lineObj* seg, Candy* M)
{
	lineObj* target;

	if (seg->enda_L_Num != 0)
	{
		for (int i = 0; i < seg->enda_L_Num;i++)
		{
			target = seg->enda_L[i];
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
				if (target->enda_L_Num != 0 && target->endb_L_Num != 0)   //single change to double seg, otherwise, remain in single
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
	if (seg->endb_L_Num != 0)
	{
		for (int i = 0; i < seg->endb_L_Num;i++)
		{
			target = seg->endb_L[i];
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
				if (target->enda_L_Num != 0 && target->endb_L_Num != 0)   //single change to double seg, otherwise, remain in single
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

void UpdateItsNeighboors_death (lineObj* seg, Candy* M)
{
	lineObj* target;

	if (seg->enda_L_Num != 0)
	{
		for (int i = 0; i < seg->enda_L_Num;i++)
		{
			target = seg->enda_L[i];
			 if (target->type == 2)
			 {
				updateNeighboorLink_death(M,target,seg,seg->enda);
				if (target->enda_L_Num == 0 || target->endb_L_Num == 0)   //double change to single seg, otherwise, remain in double
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
				 if(target->enda_L_Num == 0 && target->endb_L_Num == 0)//single seg change to free seg, otherwise remain single
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
	if (seg->endb_L_Num != 0)
	{
		for (int i = 0; i < seg->endb_L_Num;i++)
		{
			target = seg->endb_L[i];
			 if (target->type == 2)
			 {
				updateNeighboorLink_death(M,target,seg,seg->endb);
				if (target->enda_L_Num == 0 || target->endb_L_Num == 0)   //double change to single seg, otherwise, remain in double
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
				 if(target->enda_L_Num == 0 && target->endb_L_Num == 0)//single seg change to free seg, otherwise remain single
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
void updateNeighboorLink_death(Candy *M, lineObj* target,lineObj* seg,site end)
{
	int if_find = 0;
	if (DIST(target->enda.x,target->enda.y,end.x,end.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
		{
			for(int i = 0; i < target->enda_C_Num;i++)
			{
				if(target->enda_C[i] == seg)
				{
					if_find = 1;
					if (i == target->enda_C_Num-1)
					{
						target->enda_C_Num--;
						
						M->connection_n--;
						KillNClinks(M->connection_l,target,seg);
						break;
					}
					for (int ii = i;ii<target->enda_C_Num-1;ii++)
					{
						target->enda_C[ii] = target->enda_C[ii+1];
					}
					target->enda_C_Num--;
					M->connection_n--;
					KillNClinks(M->connection_l,target,seg);
					break;
				}
			}
			for(int i = 0;i<target->enda_L_Num;i++)
			{
				if (target->enda_L[i] == seg)
				{
					if_find = 1;
					if(i == target->enda_L_Num-1)
					{
						target->enda_L_Num--;
						M->neighbor_n--;
						KillNClinks(M->neighbor_l,target,seg);
						break;
					}
					for(int ii = i;ii<target->enda_L_Num-1;ii++)
					{
						target->enda_L[ii] = target->enda_L[ii+1];
					}
					target->enda_L_Num--;
					M->neighbor_n--;
					KillNClinks(M->neighbor_l,target,seg);
					break;
				}
	
			}
		}
		 if (DIST(target->endb.x,target->endb.y,end.x,end.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
		{
			for(int i = 0; i < target->endb_C_Num;i++)
			{
				if(target->endb_C[i] == seg)
				{
					if_find = 1;
					if (i == target->endb_C_Num-1)
					{
						target->endb_C_Num--;
						M->connection_n--;
						KillNClinks(M->connection_l,target,seg);
						break;
					}
					for (int ii = i;ii<target->endb_C_Num-1;ii++)
					{
						target->endb_C[ii] = target->endb_C[ii+1];
					}
					target->endb_C_Num--;
					M->connection_n--;
					KillNClinks(M->connection_l,target,seg);
					break;
				}
			}
			for(int i = 0;i<target->endb_L_Num;i++)
			{
				if (target->endb_L[i] == seg)
				{
					if_find = 1;
					if(i == target->endb_L_Num-1)
					{
						target->endb_L_Num--;
						M->neighbor_n--;
						KillNClinks(M->neighbor_l,target,seg);
						break;
					}
					for(int ii = i;ii<target->endb_L_Num-1;ii++)
					{
						target->endb_L[ii] = target->endb_L[ii+1];
					}
					target->endb_L_Num--;
					M->neighbor_n--;
					KillNClinks(M->neighbor_l,target,seg);
					break;
				}
			}
		}
}

void updateNeighboorLink_born(Candy *M, lineObj* target,lineObj* seg,site end)
{
	if (DIST(target->enda.x,target->enda.y,end.x,end.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
		{
			if (target->enda.x == end.x && target->enda.y == end.y)
			{
				target->enda_C[target->enda_C_Num] = seg;
				target->enda_C_Num++;

				M->connection_n++;
				AddNClinks(M->connection_l,M->connection_n,target, seg);
			}
			else
			{
				target->enda_L[target->enda_L_Num] = seg;
				target->enda_L_Num++;

				M->neighbor_n++;
				AddNClinks(M->neighbor_l,M->neighbor_n,target, seg);
			}
		}
		if (DIST(target->endb.x,target->endb.y,end.x,end.y) < NEIGHBOORHOOD*NEIGHBOORHOOD)
		{
			if (target->endb.x == end.x && target->endb.y == end.y)
			{
				target->endb_C[target->endb_C_Num] = seg;
				target->endb_C_Num++;

				M->connection_n++;
				AddNClinks(M->connection_l,M->connection_n,target, seg);
			}
			else
			{
				target->endb_L[target->endb_L_Num] = seg;
				target->endb_L_Num++;

				M->neighbor_n++;
				AddNClinks(M->neighbor_l,M->neighbor_n,target, seg);
			}
		}
}




LinkedList LinkedListInit()
{
    Node *L;
    L = (Node *)malloc(sizeof(Node));  
    if(L == NULL)                     
        printf("fail to allocate node\n");
    L->next = NULL;   

	return L;
}

void LinkedListInsert(LinkedList L,int i, lineObj* x)
{
    Node *pre;                     
    pre = L;
    int tempi = 0;
    for (tempi = 1; tempi < i; tempi++)
        pre = pre->next;                
    Node *p;                                
    p = (Node *)malloc(sizeof(Node));
    p->index = x; 
    p->next = pre->next;
    pre->next = p;
                             
} 

void LinkedListDelete(LinkedList L,lineObj* x)
{
    Node *p,*pre;
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




lineObj* ReturnObj (LinkedList L, int i)
{
	Node *pre;
	pre = L;
	int tempi = 0;
	for (tempi = 1;tempi<i;tempi++)
		pre=pre->next;
	return ((pre->next)->index);
}

/*******************************************************
GenerateEndb
Give Enda, it generates another end at location enda + length
in theta direction

Inputs: 													
enda       : First end
len        : Length for new segment
theta      : Theta for new segment
img_height : Input Image Height
img_width  : Input Image Width
*********************************************************/
site GenerateEndb(site enda,int len,double theta,int img_height,int img_width)
{
	site end;

	if(random2()<0.5)
	{
		end.x = (int)floor(enda.x + len*cos(theta)+0.5);
		end.y = (int)floor(enda.y - len*sin(theta)+0.5);
		if(end.x < 0)
		{
			end.x = 0;
			end.y = (int)floor(enda.y+enda.x * tan(theta)+0.5);
			if (end.y < 0)
			{
				end.y = 0;
				end.x = (int)floor (enda.x+enda.y/tan(theta)+0.5);
			}
		}
		if (end.x >= img_width)
		{
			end.x = img_width-1;
			end.y = (int)floor(enda.y -(end.x-enda.x)*tan(theta)+0.5);
			if (end.y < 0)
			{
				end.y = 0;
				end.x = (int)floor(enda.x + enda.y/tan(theta)+0.5);
			}
		}
		if (end.y < 0)
		{
			end.y = 0;
			end.x = (int)floor (enda.x+enda.y/tan(theta)+0.5);
		}
	}
	else
	{
		end.x = (int)floor(enda.x + len*cos(theta+_PI)+0.5);
		end.y = (int)floor(enda.y - len*sin(theta+_PI)+0.5);
		if (end.x < 0)
		{
			end.x = 0;
			end.y = (int)floor(enda.y+enda.x * tan(theta)+0.5);
			if (end.y >= img_height)
			{
				end.y = img_height-1;
				end.x = (int)floor (enda.x-(end.y-enda.y)/tan(theta)+0.5);
			}
		}
		if (end.x >= img_width)
		{
			end.x = img_width-1;
			end.y = (int)floor(enda.y -(end.x-enda.x)*tan(theta)+0.5);
			if (end.y >= img_height)
			{
				end.y = img_height-1;
				end.x = (int)floor(enda.x - (end.y-enda.y)/tan(theta)+0.5);
			}
		}
		if (end.y >= img_height)
		{
			end.y = img_height-1;
			end.x = (int)floor (enda.x-(end.y-enda.y)/tan(theta)+0.5);
		}
	}

	return end;
}

/* 
	Return Distance Squared  between points a and b
*/
int DistSquare(site a, site b)
{
	return((a.x - b.x)*(a.x - b.x)+(a.y - b.y)*(a.y - b.y));
}


/*******************************************************
Check Connection:
Given 2 segments seg1:(a-b) and seg2:(aa-bb), checks if dist(a,aa)/dist(a,bb)/dist(b,aa)/dist(b,bb)
are within |neighboorR| of each other.
If so, then checks if seg2 is included in the ball of radious |neighboorR| of the node a/b of seg1 (or if seg1 
is included in the ball of radious |neighboorR| of node aa/bb of seg 2)

Inputs:
	site a 		: Coordinates x,y of point a of seg1.
	site b 		: Coordinates x,y of point b of seg1.
	site aa 	: Coordinates x,y of point aa of seg2.
	site bb 	: Coordinates x,y of point bb of seg2.
	searchRatio : Radious of search as a proportion to min(|seg1|,|seg2|).  
	neighboorR 	: Radious of search as a fixed quantity.

Outputs:
	if_eo 				  : 
							Returns 1 if dist(a,aa)/dist(a,bb)/dis(b,aa)/dist(b,bb) is small
	if_connect 			  : 
							Returns 1 if dist(a,aa) or dist(a,bb) is small 
							Returns 2 if dist(b,aa) or dist(b,bb) is small
							Returns 0 otherwise 

	which_end_is_neighboor: 
							Returns 1 if ONLY dist(a,aa)/dist(a,bb) is small
							Returns 2 if ONLY dist(b,aa)/dist(b,bb) is small
							Returns 3 if more than one dis(x,y) is small. 
Returns:
	Minimun distance between (a,b) and (aa,bb)
*********************************************************/

double CheckConnection(site a, site b, site aa, site bb, double searchRatio, int neighboorR,int *if_eo,
					 int *if_connect, int *which_end_is_neighboor)
{

	*which_end_is_neighboor = 0;
	
	double len_a = sqrt(double(DistSquare(a,b)));
	double len_aa = sqrt(double(DistSquare(aa,bb)));
	double dist;
	int searchR = int(MIN(len_a,len_aa)*searchRatio);
	
	if (searchR == 0) 
		searchR = 1; 

	searchR = neighboorR;

	//1 Check if distance between node a and aa is less than searchR
	if ((dist = DistSquare(a,aa)) < searchR*searchR )
	{
		*if_eo = 1;

		//Check if they are connected
		if ((a.x == aa.x) && (a.y == aa.y))
			*if_connect = 1;

		//Check if distance between node a and bb is more than searchR
		//Otherwise segment (aa-bb) is included in ball searchR around (a-b)
		if (DistSquare(a,bb) > searchR*searchR)
		{
			//Check if distance between node b and aa is more than searchR
			//Otherwise segment (a-b) is included in ball searchR around (aa-bb)
			if (DistSquare(b,aa) > searchR*searchR)
				*which_end_is_neighboor = 1;
			else
				*which_end_is_neighboor = 3; 
		}
		else
		{
			*which_end_is_neighboor = 3; 
		}
	}
	//2 Check if distance between node a and bb is less than searchR
	else if ((dist = DistSquare(a,bb)) < searchR*searchR)
	{
		*if_eo = 1;
		if((a.x == bb.x)&&(a.y == bb.y))
			*if_connect = 1; 

		//Check if distance between node a and aa is more than searchR
		//Otherwise segment (aa-bb) is included in ball searchR around (a-b)
		if (DistSquare(a,aa) > searchR*searchR)
			if (DistSquare(b,bb) > searchR*searchR)
				*which_end_is_neighboor = 1;
			else
				*which_end_is_neighboor = 3; 
		else
			*which_end_is_neighboor = 3; 
	}

	//3 Check if distance between node b and aa is less than searchR
	else if ((dist = DistSquare(b,aa)) < searchR*searchR)
	{
		*if_eo = 1;

		//Check if they are connected
		if((b.x == aa.x) && (b.y == aa.y))
			(*if_connect) = 2; 

		if (DistSquare(b,bb) > searchR*searchR)
			if (DistSquare(a,aa) > searchR*searchR)
				*which_end_is_neighboor = 2;
			else
				*which_end_is_neighboor = 3; //KDW 3 5
		else
			*which_end_is_neighboor = 3; //KDW 3 5
	}
	else if ((dist = DistSquare(b,bb)) < searchR*searchR)// && DistSquare(b,aa) > searchR*searchR )
	{
		*if_eo = 1;
		if ((b.x == bb.x) && (b.y == bb.y))
			(*if_connect) = 2; //KDW 2 4

		if (DistSquare(b,aa) > searchR*searchR)
			if (DistSquare(a,bb) > searchR*searchR)
				*which_end_is_neighboor = 2;
			else
				*which_end_is_neighboor = 3; //KDW 3 5
		else
			*which_end_is_neighboor = 3; //KDW 3 5
	}
	return dist;
}




void line(int x1, int y1, int x2, int y2, int color, double **img)
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
	} else {
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


void DrawLine(double **img,Candy *C)
{
	int n_f = C->n_f;
	int n_s = C->n_s;
	int n_d = C->n_d;

	site enda, endb;

	Node *p = C->link_f;

	for (int i = 0;i<n_f;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;

		line(enda.x,enda.y,endb.x,endb.y,255,img);
	}

	p = C->link_s;
	for (int i = 0;i<n_s;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;
		line(enda.x,enda.y,endb.x,endb.y,255,img);
	}

	p = C->link_d;
	for (int i = 0;i<n_d;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;
		line(enda.x,enda.y,endb.x,endb.y,255,img);
	}
}


NeckDent *SaveCandy2MP(Candy *C, int *num)
//int SaveCandy2MP(Candy *C, NeckDent **mp)
{
	int n_f = C->n_f;
	int n_s = C->n_s;
	int n_d = C->n_d;

	site enda, endb;

	NeckDent *mp=NULL;

	mp = new NeckDent [n_f + n_s + n_d];
	//mp = (NeckDent *)malloc( (n_f + n_s + n_d) *sizeof(NeckDent));
	Node *p = C->link_f;
	int j=0;
	for (int i = 0;i<n_f;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;

		mp[j].num			= j;
		mp[j].state			= STATE_EXIST;
		mp[j].type			= NDTYPE_CANDY;
		mp[j].center.x		= p->index->x;
		mp[j].center.y		= p->index->y;
		mp[j].enda.x		= p->index->enda.x;
		mp[j].enda.y		= p->index->enda.y;
		mp[j].endb.x		= p->index->endb.x;
		mp[j].endb.y		= p->index->endb.y;
		mp[j].width			= p->index->width;
		mp[j].length		= p->index->len;
		mp[j].theta			= -p->index->theta;
		mp[j].single_E		= p->index->dataterm;
		mp[j].multiple_E	= 0;
		mp[j].sort_idx		= 0;
		mp[j].e[0]			= 0;
		mp[j].e[1]			= 0;
		mp[j].e[2]			= 0;
		mp[j].e[3]			= 0;
		j++;
	}



	p = C->link_s;
	for (int i = 0;i<n_s;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;

		mp[j].num			= j;
		mp[j].state			= STATE_EXIST;
		mp[j].type			= NDTYPE_CANDY;
		mp[j].center.x		= p->index->x;
		mp[j].center.y		= p->index->y;
		mp[j].enda.x		= p->index->enda.x;
		mp[j].enda.y		= p->index->enda.y;
		mp[j].endb.x		= p->index->endb.x;
		mp[j].endb.y		= p->index->endb.y;
		mp[j].width			= p->index->width;
		mp[j].length		= p->index->len;
		mp[j].theta			= -p->index->theta;
		mp[j].single_E		= p->index->dataterm;
		mp[j].multiple_E	= 0;
		mp[j].sort_idx		= 0;
		mp[j].e[0]			= 0;
		mp[j].e[1]			= 0;
		mp[j].e[2]			= 0;
		mp[j].e[3]			= 0;
		j++;
	}

	p = C->link_d;
	for (int i = 0;i<n_d;i++)
	{
		p= p->next;
		enda = p->index->enda;
		endb = p->index->endb;

		mp[j].num			= j;
		mp[j].state			= STATE_EXIST;
		mp[j].type			= NDTYPE_CANDY;
		mp[j].center.x		= p->index->x;
		mp[j].center.y		= p->index->y;
		mp[j].enda.x		= p->index->enda.x;
		mp[j].enda.y		= p->index->enda.y;
		mp[j].endb.x		= p->index->endb.x;
		mp[j].endb.y		= p->index->endb.y;
		mp[j].width			= p->index->width;
		mp[j].length		= p->index->len;
		mp[j].theta			= -p->index->theta;
		mp[j].single_E		= p->index->dataterm;
		mp[j].multiple_E	= 0;
		mp[j].sort_idx		= 0;
		mp[j].e[0]			= 0;
		mp[j].e[1]			= 0;
		mp[j].e[2]			= 0;
		mp[j].e[3]			= 0;
		j++;
	}

	*num = j;
	return mp;
}




// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(site p, site q, site r)
{
    if (q.x <= MAX(p.x, r.x) && q.x >= MIN(p.x, r.x) &&
        q.y <= MAX(p.y, r.y) && q.y >= MIN(p.y, r.y))
       return true;
 
    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(site p, site q, site r)
{
    // See 10th slides from following link for derivation of the formula
    // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
 
    if (val == 0) return 0;  // colinear
 
    return (val > 0)? 1: 2; // clock or counterclock wise
}

bool doIntersect(site p1, site q1, site p2, site q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
 
    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    return false; // Doesn't fall in any of the above cases
}

double theta_from_two_ends(site enda, site endb)
{
	double theta;
	if(enda.y < endb.y)
	{
		if (enda.x-endb.x != 0)
		 theta = atan(double((endb.y-enda.y)/double(enda.x-endb.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;
	}
	else
	{
		if (enda.x-endb.x != 0)
		 theta = atan(double((enda.y-endb.y)/double(endb.x-enda.x)));
		else
			theta = _PI/2;
		if (theta < 0)
			theta += _PI;

	}
	return theta;
}

/*******************************************************
if_connect_a_freeEnd
Checks in end is connected to a free end of object obj

Inputs: 													
	end   : End Connecting to Object
	obj   : Line Object where end is connecting to

*********************************************************/
int if_connect_a_freeEnd(site end, lineObj *obj)
{
	int n = 0;
	if(obj->enda_C_Num == 0)
	{
		if((end.x == obj->enda.x)&&(end.y == obj->enda.y))
			n = 1;
	}
	if (obj->endb_C_Num == 0)
	{
		if((end.x == obj->endb.x)&&(end.y == obj->endb.y))
			n = 1;
	}

	return n;
}

int if_connect_a_connectedEnd(site end, lineObj *obj)
{
	int n = 0;
	if(obj->enda_C_Num == 1)
	{
		if((end.x == obj->enda.x)&&(end.y == obj->enda.y))
			n = 1;
	}
	if (obj->endb_C_Num == 1)
	{
		if((end.x == obj->endb.x)&&(end.y == obj->endb.y))
			n = 1;
	}

	return n;
}


/*******************************************************
Echange_from_neighboors_born
Once a line object is born, it changes the state of its neighbors
		ie. change free segment to single segment

	Inputs: 													
		M     : Linked List With All Objects.
		obj   : New object to be added. 
 
*********************************************************/
double Echange_from_neighboors_born(Candy *M, lineObj *obj)
{
	double w_s = M->w_s;
	double w_f = M->w_f;
	double w_d = M->w_d;

	double Echange = 0;
	if (obj->enda_C_Num != 0)
	{
		for (int i = 0;i<obj->enda_C_Num;i++)
		{
			if(obj->enda_C[i]->type == 0)  //a previous freeseg now changes to single
			{
				Echange += w_s-w_f;
			}
			else if(obj->enda_C[i]->type == 1)
			{
				if (if_connect_a_freeEnd(obj->enda,obj->enda_C[i]) == 1)  //a previous singseg now changes to double
					Echange += w_d-w_s;
			}
			else if(obj->enda_C[i]->type == 2)
			{
				Echange +=0;
			}

		}
	}

	if (obj->endb_C_Num != 0)
	{
		for (int i = 0;i<obj->endb_C_Num;i++)
		{
			if(obj->endb_C[i]->type == 0)  //a previous freeseg now changes to single
			{
				Echange += w_s-w_f;
			}
			if(obj->endb_C[i]->type == 1)
			{
				if (if_connect_a_freeEnd(obj->endb,obj->endb_C[i]) == 1)  //a previous singseg now changes to double
					Echange += w_d-w_s;
			}
		}
	}
	return Echange;
}

double Echange_from_neighboors_death(Candy *M, lineObj *obj)
{
	double w_s = M->w_s;
	double w_f = M->w_f;
	double w_d = M->w_d;

	double Echange = 0;
	if (obj->enda_C_Num != 0)
	{
		for (int i = 0;i<obj->enda_C_Num;i++)
		{
			if(obj->enda_C[i]->type == 2)  
			{
				if (if_connect_a_connectedEnd(obj->enda,obj->enda_C[i]) == 1)  //a previous double now changes to single
					Echange += w_s-w_d;
			}
			if(obj->enda_C[i]->type == 1)
			{
				if (if_connect_a_connectedEnd(obj->enda,obj->enda_C[i]) == 1)  //a previous singseg now changes to free
					Echange += w_f-w_s;
			}
		}
	}
	if (obj->endb_C_Num != 0)
	{
		for (int i = 0;i<obj->endb_C_Num;i++)
		{
			if(obj->endb_C[i]->type == 2)  
			{
				if (if_connect_a_connectedEnd(obj->endb,obj->endb_C[i]) == 1)  //a previous double now changes to single
					Echange += w_s-w_d;
			}
			if(obj->endb_C[i]->type == 1)
			{
				if (if_connect_a_connectedEnd(obj->endb,obj->endb_C[i]) == 1)  //a previous singseg now changes to free
					Echange += w_f-w_s;
			}
		}
	}

	//Echange = 0;
	return Echange;
}

void count_all_ends_in_link(LinkedList link, int n, site *end, lineObj **obj, int *num)
{
	Node *p = link;
	int found,j;
	site endx;

	for (int i = 0;i<n;i++)
	{
		p= p->next;
		endx = p->index->enda;
		found =0;
		for(j=0; j<(*num); j++)
			if((endx.x == end[j].x)&&(endx.y == end[j].y)) found = 1;
		if(!found){
			end[j] = endx;
			obj[j] = p->index;
			(*num)++;
		}
		endx = p->index->endb;
		found =0;
		for(j=0; j<(*num); j++)
			if((endx.x == end[j].x)&&(endx.y == end[j].y)) found = 1;
		if(!found){
			end[j] = endx;
			obj[j] = p->index;
			(*num)++;
		}
	}

}
int count_all_ends(Candy *M, site *end, lineObj **obj)
{
	int num=0;

	count_all_ends_in_link(M->link_f,M->n_f,end,obj,&num);
	count_all_ends_in_link(M->link_s,M->n_s,end,obj,&num);
	count_all_ends_in_link(M->link_d,M->n_d,end,obj,&num);

	return num;

}

/* Counts the Number of "Free" ends found so far in list link*/
void count_all_ends_in_link(LinkedList link, int n, site *end, int *num)
{
	Node *p = link;
	int found,j;
	site endx;


	for (int i = 0;i<n;i++)
	{
		/*Count one one side*/
		p= p->next;
		endx = p->index->enda;
		found =0;
		for(j=0; j<(*num); j++)
		{
			if((endx.x == end[j].x)&&(endx.y == end[j].y))
				found = 1;
		}
		if(!found)
		{
			end[j] = endx;
			(*num)++;
		}
		/*Count one side b*/
		endx = p->index->endb;
		found =0;
		for(j=0; j<(*num); j++)
		{
			if((endx.x == end[j].x)&&(endx.y == end[j].y))
				found = 1;
		}
		if(!found)
		{
			end[j] = endx;
			(*num)++;
		}
	}

}

/* Counts the Number of "Free" ends found so far*/
int count_all_ends(Candy *M, site *end)
{
	int num=0;

	count_all_ends_in_link(M->link_f,M->n_f,end,&num);
	count_all_ends_in_link(M->link_s,M->n_s,end,&num);
	count_all_ends_in_link(M->link_d,M->n_d,end,&num);

	return num;

}
