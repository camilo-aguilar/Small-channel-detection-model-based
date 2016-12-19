#include "QualityCandy.h"
#include "tiff.h"
#include "allocate.h"
#include <iostream> //KDW
#include <time.h>
#include "neck.h"


using namespace std;

void DEBUG_IMAGE(char* name, struct TIFF_img inter_img_mrf)
{
	FILE *fp;
			/* open image file */
			if ( ( fp = fopen ( name, "wb" ) ) == NULL ) 
			{
				fprintf ( stderr, "cannot open file image.tif\n");
				exit ( 1 );
			}

			/* write image */
			if ( write_TIFF ( fp, &inter_img_mrf ) ) 
			{
			fprintf ( stderr, "error writing TIFF file \n" );
			exit ( 1 );
			}

			/* close image file */
			fclose ( fp );
}


void DEBUG_ARRAY(char* name, double** my_image, int height, int width)
{	
	
	int i = 0, j=0;
	struct TIFF_img inter_img_mrf;
	get_TIFF ( &inter_img_mrf, height, width, 'g' );
	int pixel;
	for(i = 0; i < height; i++)
	{
		for(j=0; j< width; j++)
		{
			pixel = (int32_t)(my_image[i][j]);
			inter_img_mrf.mono[i][j] = pixel;
			
		}
	}	
	DEBUG_IMAGE(name, inter_img_mrf);
}



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
	Candy *C;
	C = (Candy *)malloc(sizeof(Candy));

	C->beta = 1;					//Beta
	C->lambda = 10;					//Lambda

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

	C->p_f_b = 0.5; 				/* Free Segment Birth */
	C->p_f_d = 0.5; 				/* Free Segment Death */
	
	C->p_s_b = 0.5; 				/* Single Segment Birth */
	C->p_s_d = 0.5; 				/* Single Segment Death */
	
	C->p_d_b = 0.5; 				/* Double Segment Birth */
	C->p_d_d = 0.5; 				/* Double Segment Death */


	C->gamma_d = mpp.gamma_d; 		//
	
	C->w_eo = mpp.w_eo; 			/* Weight External Orientation */
	C->w_f = mpp.w_f;   			/* Weight free Seg */
	C->w_s = mpp.w_s;   			/* Weight Single Seg */
	C->w_d = mpp.w_d;   			/* Weight Double Seg */
	C->w_io = mpp.w_io; 			/* Weight Internal Bad Orientation */

	C->p_b_d = 0.6;  				/* Birth and Death Step */
	C->p_t = 0.2;     				/* Translation Step */
	C->p_c = 0.2;     				/* Connection Step */

	C->p_c_CtoF = 0.3; 				/* Connected to Free */
	C->p_c_FtoC = 0.7; 				/* Free to Connected */

	C->p_t_f = 0.5; 				/* Transition Free Seg*/
	C->p_t_s = 0.5; 				/* Transition Single Seg*/
	C->p_t_d = 0.3; 				/* Transition Double Seg*/

	C->p_f = 0.6; 					/* Free (B&D) */
	C->p_s = 0.2;  					/* Single */
	C->p_d = 0.2;  					/*Double */

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
	int new_x = int(start_end.x -3 *cos(theta));
	int new_y = int(start_end.y+3*sin(theta));

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

	struct TIFF_img inter_img_mrf;
	FILE *fp;

	char name_buff2[2000];
	get_TIFF ( &inter_img_mrf, height, width, 'g' );

	Candy *C;
	C = CandyInit(mpp);

	double T = 1.0;

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
		
		/* Write Image for debugging */
		if(0)
		{	
			DEBUG_IMAGE("original_seg.tiff", inter_img_mrf);

		}
	 }
	
	int test = 0;
	time_t t;
	clock_t start_time=clock();

	T = mpp.T0;
	printf("Start RJMCMC\n");
	for (int i = 0;i<mpp.iter_num;i++)
	{
		if(i > 0 && i % (mpp.iter_num/10) == 0) printf(" %2.2f Percent Complete \n", (float)i/(float)mpp.iter_num * 100.00);
		double r = random2();
		T = T/(mpp.de_coeff); //KDW

		if (r <C->p_b_d)   //birth and death step
		{
			double rr = random2();
			double rrr = random2();
			if (rr < C->p_f)  // free segment
			{
				if (rrr < C->p_f_b)
				{
				
					AddFreeSeg(C,input_img,lm,output_seg, img_mpp_l,img_seg_l,height,width,T,test,patch,patch_len, Matrix, prior_pen1, prior_pen2,dilated_img);
					test++;
				}
				else
				{
					KillFreeSeg(C,input_img,height,width,T);
				}
			}
			else if (rr < C->p_s+C->p_f)  //single segment
			{
				if (rrr < C->p_s_b)
				{
					//AddSingleSeg_allends(C,input_img,height,width,T);
					AddSingleSeg(C,input_img,lm,output_seg, img_mpp_l,img_seg_l,height,width,T,patch,patch_len, Matrix, prior_pen1, prior_pen2);
				}
				else
				{
					//KillSingleSeg_allends(C,input_img,height,width,T);
					KillSingleSeg(C,input_img,height,width,T);
				}
			}
			else
			{
				if (rrr < C->p_d_b)
				{
					//AdddoubleSeg_allends(C,input_img,height,width,T);
					AdddoubleSeg(C,input_img,lm,output_seg,img_mpp_l,img_seg_l,height,width,T,patch,patch_len, Matrix, prior_pen1, prior_pen2,0);
				}
				else
				{
					//KillDoubleSeg_allends(C,input_img,height,width,T);
					KillDoubleSeg(C,input_img,height,width,T);
				}
			}
		}
		else if (r <C->p_b_d +C->p_t)   //transition step
		{
			double rr = random2();
			double rrr = random2();
			if (rr < C->p_t_f)
			{
				if (rrr < 0.2)
					FreeSeg_length_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height,width,T,patch,patch_len,Matrix, prior_pen1, prior_pen2);
				else if (rrr< 0.4)
					FreeSeg_theta_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height,width,T,patch,patch_len,Matrix, prior_pen1, prior_pen2);
				else  if (rrr< 0.6) //KDW 
					FreeSeg_width_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height,width,T,patch,patch_len,Matrix, prior_pen1, prior_pen2);
				else  if (rrr< 0.8)  //KDW
					FreeSeg_freeEnd_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);
				else //KDW
					FreeSeg_center_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);

			}
			else if (rr < C->p_t_s+C->p_t_f)
			{
				 SingleSeg_freeEnd_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);

			}
			else  // no disconnect & move
			{
				 SingleDoubleSeg_Connection_move(C, input_img, output_seg, img_mpp_l,img_seg_l, height, width, T, patch,patch_len,Matrix, prior_pen1, prior_pen2);
				 
			}
		}
		else   //connetion and separation step
		{
			double rr = random2();
			if (rr < C->p_c_FtoC)
			{
				Conncet_freeEnds(C,input_img,output_seg, img_mpp_l,img_seg_l,patch,patch_len,height,width,Matrix, prior_pen1,prior_pen2,T);
			}
			else
			{
				Seperate_connectedEnds(C,input_img,output_seg, img_mpp_l,img_seg_l,patch,patch_len,height,width,Matrix, prior_pen1,prior_pen2,T);
			}

		}

		
		
		//DELETED PARTS
		if (i!= 0 && i % 1000000 == 0 && 0)
		{
			printf("connections %d, neighbors %d\n",C->connection_n,C->neighbor_n);
			int pixel,pixel2;

			writeC(C,i);
			DrawLine_back(output_img,C);
			DrawLine_Seg(tmp_img,C,height,width,mu,cov,patch,patch_len);	

			for ( int ii = 0; ii < height; ii++ )
			{
				for (int jj = 0; jj < width; jj++ ) 
				{
					dilated_img[ii][jj] = 0;
				}
			}

			dilation(output_seg,dilated_img,width,height);


			FILE *fp;
			struct TIFF_img I_out,I_out_seg;
			get_TIFF ( &I_out, height, width, 'g' );
			get_TIFF(&I_out_seg,height,width, 'g');

		  for ( int ii = 0; ii < height; ii++ )
			for (int jj = 0; jj < width; jj++ ) 
			{
				pixel = (int)output_img[ii][jj];
			 	pixel2 = (int)output_seg[ii][jj]*256/3;
			if(pixel>255) {
				 I_out.mono[ii][jj] = 255;
				}
			else {
				if(pixel<128) I_out.mono[ii][jj] = (unsigned char)input_img[ii][jj]; //KDW
				else I_out.mono[ii][jj] = pixel;
			}
			if(pixel2>255) {
				 I_out_seg.mono[ii][jj] = 255;
				}
			else {
			 if(pixel2<0) I_out_seg.mono[ii][jj] = 0;
			 else I_out_seg.mono[ii][jj] = pixel2;
			}
			  }
				char a[360],aa[360];
				sprintf(a,"test_out%d_all_seg_520_HP_after_mpp_and_seg_new_offset1_t_%f_%f.tiff",i,PEN_1,PEN_2);
				sprintf(aa,"seg_out%d_all_seg_520_HP_after_mpp_and_seg_new_offset%f_t_%f_%f.tiff",i,SEG_OFFSET,PEN_1,PEN_2);
				/* open image file */
				if ( ( fp = fopen ( a, "wb" ) ) == NULL ) {
				fprintf ( stderr, "cannot open file MBD_img.tif\n");
				exit ( 1 );
			  }

				/* write image */
				if ( write_TIFF ( fp, &I_out ) ) {
				fprintf ( stderr, "error writing TIFF file" );
				exit ( 1 );
			  }

				/* close image file */
				fclose ( fp );
				if ( ( fp = fopen ( aa, "wb" ) ) == NULL ) {
				fprintf ( stderr, "cannot open file MBD_img.tif\n");
				exit ( 1 );
			  }
			  /* write image */
			  if ( write_TIFF ( fp, &I_out_seg ) ) {
				fprintf ( stderr, "error writing TIFF file" );
				exit ( 1 );
			  }

			  /* close image file */
			  fclose ( fp );
			  free_TIFF ( &(I_out) );
			  free_TIFF ( &(I_out_seg) );

		  for ( int ii = 0; ii < height; ii++ )
			for (int jj = 0; jj < width; jj++ ) 
				{
					output_img[ii][jj] = 0;
					tmp_img[ii][jj] = 0;
			}
			//update the segmentation

			printf("iter = %d\n",i / 1000000);
		  printf("free=%d,single=%d,double=%d\n",C->n_f,C->n_s,C->n_d);
		}

		
		
		
		if(i / 1000000 == 20)
		{
			clock_t end_time=clock();
			cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<endl;

		}

	}
	
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
     
   // return L;                           
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


int DistSquare(site a, site b)
{
	return((a.x - b.x)*(a.x - b.x)+(a.y - b.y)*(a.y - b.y));
}


double CheckConnection(site a, site b, site aa, site bb, double searchRatio, int neighboorR,int *if_eo,
					 int *if_connect, int *which_end_is_neighboor)
{
	double len_a = sqrt(double(DistSquare(a,b)));
	double len_aa = sqrt(double(DistSquare(aa,bb)));
	double dist;

	int searchR = int(MIN(len_a,len_aa)*searchRatio);
	*which_end_is_neighboor = 0;
	//if (searchR <= 5)
	//KDW	searchR = 5;
	if (searchR == 0) //KDW
		searchR = 1; //KDW

	searchR = neighboorR;

	if ((dist = DistSquare(a,aa)) < searchR*searchR )//&& DistSquare(a,bb) > searchR*searchR )
	{
		*if_eo = 1;
		if ((a.x == aa.x) && (a.y == aa.y))
			*if_connect = 1;

		if (DistSquare(a,bb) > searchR*searchR)
			if (DistSquare(b,aa) > searchR*searchR)
				*which_end_is_neighboor = 1;
			else
				*which_end_is_neighboor = 3; //KDW 3 5
		else
			*which_end_is_neighboor = 3; //KDW 3 5
	}
	else if ((dist = DistSquare(a,bb)) < searchR*searchR)// && DistSquare(a,aa) > searchR*searchR )
	{
		*if_eo = 1;
		if((a.x == bb.x)&&(a.y == bb.y))
			*if_connect = 1; //KDW 1 2

		if (DistSquare(a,aa) > searchR*searchR)
			if (DistSquare(b,bb) > searchR*searchR)
				*which_end_is_neighboor = 1;
			else
				*which_end_is_neighboor = 3; //KDW 3 5
		else
			*which_end_is_neighboor = 3; //KDW 3 5
	}
	else if ((dist = DistSquare(b,aa)) < searchR*searchR)// && DistSquare(b,bb) > searchR*searchR )
	{
		*if_eo = 1;
		if((b.x == aa.x) && (b.y == aa.y))
			(*if_connect) = 2; //KDW 2 3

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
				Echange += w_s-w_f;
			if(obj->enda_C[i]->type == 1)
			{
				if (if_connect_a_freeEnd(obj->enda,obj->enda_C[i]) == 1)  //a previous singseg now changes to double
					Echange += w_d-w_s;
			}
		}
	}
	if (obj->endb_C_Num != 0)
	{
		for (int i = 0;i<obj->endb_C_Num;i++)
		{
			if(obj->endb_C[i]->type == 0)  //a previous freeseg now changes to single
				Echange += w_s-w_f;
			if(obj->endb_C[i]->type == 1)
			{
				if (if_connect_a_freeEnd(obj->endb,obj->endb_C[i]) == 1)  //a previous singseg now changes to double
					Echange += w_d-w_s;
			}
		}
	}
	//Echange = 0;
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
void count_all_ends_in_link(LinkedList link, int n, site *end, int *num)
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
			(*num)++;
		}
		endx = p->index->endb;
		found =0;
		for(j=0; j<(*num); j++)
			if((endx.x == end[j].x)&&(endx.y == end[j].y)) found = 1;
		if(!found){
			end[j] = endx;
			(*num)++;
		}
	}

}
int count_all_ends(Candy *M, site *end)
{
	int num=0;

	count_all_ends_in_link(M->link_f,M->n_f,end,&num);
	count_all_ends_in_link(M->link_s,M->n_s,end,&num);
	count_all_ends_in_link(M->link_d,M->n_d,end,&num);

	return num;

}
