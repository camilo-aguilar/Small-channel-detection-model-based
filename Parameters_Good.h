#ifndef _PARAMETERS_
#define _PARAMETERS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

 
/*Input Parameters for Candy Model*/
typedef struct _input_params 
{
  int iterations;
  double gamma_d;//34;
  double w_eo;//1
  double w_f;//1
  double w_s;//1
  double w_d;//1
  double w_io;//1
} INPUT_PARAMS;


/* Function to Parse Input Parameters*/
INPUT_PARAMS parse_input_qc(int argc,char** argv);


#define INCLUDE_OPENCV 1



#define _PI					3.141592654
#define PRE_SEG_EM_MPM		
/****************************************************
 *			GEOMETRIC  PARAMETERS 					* 
 ****************************************************/
#define L_MAX							 50 
#define L_MIN							 10 
#define W_MIN							 3 
#define W_MAX							 5 
#define THETA_MIN						(_PI/2.0 - _PI/10)
#define THETA_MAX						(_PI/2.0 + _PI/10)


/****************************************************
 *				IMAGE BASED	PARAMETERS   			* 
 ****************************************************/
//#define ERROR_TH  	 	5
#define ERROR_TH  	 					 8

/****************************************************
 *	      CANDY MPP REGULARIZATION Parameters		* 
 ****************************************************/
#define GAMMA_D 		50.0

#define W_F 			10
#define W_S 			-10
#define W_D 			-20

#define W_EO 			10.0
#define W_IO 			1

#define SYM_TH			 0.5

//#define EN_NEW_SEG_AT_NON_FREE_END

/****************************************************
 *				QUALITY CANDY RJMCMC PARAMETERS		* 
 ****************************************************/
#define ITERATIONS 				 	 	  10000000

#define INITIAL_T0 				 	  	  0.1


#define DECREASE_COEFFICIENT 	 	  	  0.999
#define BETA_MPP                 	      14    

/* DEATH = 1 - BIRTH */
#define P_BIRTH_F 				 	  	  0.5
#define P_BIRTH_S					 	  0.5
#define P_BIRTH_D					 	  0.5

/* All these steps must add up to 1 */
#define P_BIRTH_DEATH_STEP       	  	  0.4
#define P_TRANSLATION_STEP		  	      0.2
#define P_CONNECTION_STEP		 	      0.4

/* All these steps must add up to 1 */
#define P_PICK_F 	         		      0.6
#define P_PICK_S            		      0.2
#define P_PICK_D            		      0.2
 
/* All these steps must add up to 1 */
#define CONNETECTED_TO_FREE  	 	      0.1
#define FREE_TO_CONNECTED 		 	      0.9

/* All these steps must add up to 1 */
#define TRANSITION_FREE_SEGMENT  	      0.3
#define TRANSITION_SINGLE_SEGMENT         0.3
#define TRANSITION_DOUBLE_SEGMENT	      0.4



/****************************************************
 *				PARAMETERS FOR EM/MPM				* 
 ****************************************************/
#define BETTA_FOR_MPP 		    	 27.0
#define GAUSSIAN_TAU    		     17.0




/****************************************************
 ****************************************************
 **************************************************** 
 *****************	NECK/DENT	*********************
 *****************	   MPP		********************* 
 ****************************************************
 ****************************************************
 ****************************************************/


/****************************************************
 *	      Neck/Dent MPP REGULARIZATION Parameters	* 
 ****************************************************/



/****************************************************
 *				NECK/DENT RJMCMC PARAMETERS			* 
 ****************************************************/
/* Multiple Birth and Death */
#define P_BIRTH                  0.5
#define P_DEATH                  0.5
#define P_DILATATION             0.1
#define P_TRANSLATION     		 0.1
#define P_ROTATION 				 0.8
#define P_SWITCHING 			 0

#define LAMBDA_RJMCMC 			 30000
#define LAMBDA_L 				 0
#define LAMBDA_INT  			 0.005
#define CHANNEL_TYPES 			 1



#endif