#ifndef _PARAMETERS_
#define _PARAMETERS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "em.h"
 

#define INCLUDE_OPENCV 1
#define IMAGE_IVERTED  1 
#define NUM_WINDOWS   100
#define DERIVATIVE_LIKELY 0

#if INCLUDE_OPENCV
	#define INCLUDE_OPENCV_FREE_SEG    1
	#define INCLUDE_OPENCV_SINGLE_SEG 1
	#define INCLUDE_OPENCV_DOUBLE_SEG  1
  #define INCLUDE_OPENCV_TRANSITION  1
#endif

#define _PI					3.141592654
#define PRE_SEG_EM_MPM		
/****************************************************
 *			GEOMETRIC  PARAMETERS 					* 
 ****************************************************/
#define L_MAX							          30 
#define L_MIN							          5
#define W_MIN							          3 
#define W_MAX							          5 
#define THETA_MIN						       (_PI/2.0 - _PI/10)
#define THETA_MAX		        			 (_PI/2.0 + _PI/10)

#define TAU_MAX                     _PI/4.0
#define DELTA_IO_MIN                0//_PI/8 (0 for no crossings at all)
/****************************************************
 *				IMAGE BASED	PARAMETERS   			* 
 ****************************************************/
//#define ERROR_TH  	 	5
#define ERROR_TH  	 					       7
#define DERIVATIVE_THRESHOLD			   500


/****************************************************
 *	      CANDY MPP REGULARIZATION Parameters		* 
 ****************************************************/
#define GAMMA_D 		                  50.0
  
#define W_F 			                    10
#define W_S 			                    -10
#define W_D 			                    -20

#define W_EO 			                    10.0
#define W_IO 			                    1

/* Symmetry Coefficient For Data Energy */
#define SYM_TH			                  0.5

/* Delta Dilatation/Rotation/Translation Kernels */
#define DELTA_LEN                      2 
#define DELTA_THETA                    _PI/20
#define DELTA_WIDTH                    2 
#define END_MOVE_RANGE                 2
//#define EN_NEW_SEG_AT_NON_FREE_END

/* NOT-SO IMPORTANT PARAMETERS */
#define MAX_CONNECTION_NUM              30
#define NEIGHBOORHOOD                   2//5

/****************************************************
 *				QUALITY CANDY RJMCMC PARAMETERS		* 
 ****************************************************/
#define ITERATIONS 				 	 	          5000000
#define INITIAL_T0 				 	  	        0.8

#define DECREASE_COEFFICIENT 	 	        0.999999
#define BETA_MPP                 	      14    

/* DEATH = 1 - BIRTH */
#define P_BIRTH_F 				 	  	        0.5
#define P_BIRTH_S					 	            0.5
#define P_BIRTH_D					 	            0.5

/* All these steps must add up to 1 */
#define P_BIRTH_DEATH_STEP       	  	  0.8
#define P_TRANSLATION_STEP		  	      0.2//0.2
#define P_CONNECTION_STEP		 	          0//0.4

/* All these steps must add up to 1 */
#define P_PICK_F 	         		          0.6
#define P_PICK_S            		        0.2
#define P_PICK_D            		        0.2
 
/* All these steps must add up to 1 */
#define CONNETECTED_TO_FREE  	 	        0.1
#define FREE_TO_CONNECTED 		 	        0.9

/* All these steps must add up to 1 */
#define TRANSITION_FREE_SEGMENT  	      1
#define TRANSITION_SINGLE_SEGMENT       0//0.3
#define TRANSITION_DOUBLE_SEGMENT	      0//0.4



/****************************************************
 *				PARAMETERS FOR EM/MPM				* 
 ****************************************************/
#define BETTA_FOR_MPP 		    	        27.0
#define GAUSSIAN_TAU    		            17.0 




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


/*Input Parameters for Candy Model*/
typedef struct _input_params 
{ 

  /* Delta Dilatation/Rotation/Translation Kernels */
  //#define DELTA_LEN                      2 
  //#define DELTA_THETA                    _PI/20
  //#define DELTA_WIDTH                    2 
  //#define END_MOVE_RANGE                 2
  //#define EN_NEW_SEG_AT_NON_FREE_END

  /* NOT-SO IMPORTANT PARAMETERS */
  //#define MAX_CONNECTION_NUM              30
  //#define NEIGHBOORHOOD                   2//5


  //#define BETA_MPP                        14    

/* DEATH = 1 - BIRTH */
  //#define P_BIRTH_F                       0.5
  //#define P_BIRTH_S                       0.5
  //#define P_BIRTH_D                       0.5

/* All these steps must add up to 1 */
  //#define P_BIRTH_DEATH_STEP              0.8
  //#define P_TRANSLATION_STEP              0.2//0.2
  //#define P_CONNECTION_STEP               0//0.4

/* All these steps must add up to 1 */
  //#define P_PICK_F                        0.6
  //#define P_PICK_S                        0.2
  //#define P_PICK_D                        0.2
 
/* All these steps must add up to 1 */
  //#define CONNETECTED_TO_FREE             0.1
  //#define FREE_TO_CONNECTED               0.9

/* All these steps must add up to 1 */
  //#define TRANSITION_FREE_SEGMENT         1
  //#define TRANSITION_SINGLE_SEGMENT       0//0.3
  //#define TRANSITION_DOUBLE_SEGMENT       0//0.4

} INPUT_PARAMS;



typedef struct mpp_parameters
{

  double lengthmax;
  double lengthmin;
  double widthmin;
  double widthmax;
  double thetamin;
  double thetamax;

  double tau_max;
  double tau_min;


/****************************************************
 *        IMAGE BASED PARAMETERS        * 
 ****************************************************/

  double error_th;
  double derivative_th;



  /* Symmetry Coefficient For Data Energy */
  double symmetry_error;
/****************************************************
 *        CANDY MPP REGULARIZATION Parameters   * 
 ****************************************************/
  double gamma_d;
  double w_f;
  double w_s;
  double w_d;
  double w_io;
  double w_eo;


  /****************************************************
 *        QUALITY CANDY RJMCMC PARAMETERS   * 
 ****************************************************/
  int iter_num;
  double T0;

  double de_coeff;

  int fixed_param_num;


  /*********************************************/
  /*    OLD PARAMETERS */
  /*********************************************/
  // fixed
  int     hard_repulsion;
  double    gaussian_tau;

  double    dent_l_w_ratio; 
  double    lambda_e;
  double    lambda_a;
  double    lambda_l;
  double    lambda_s;
  double    lambda_nc;
  double    lambda_dc;
  double    lambda_int;
  double    length_th;
  double    symmetry_th;
  double    p_birth;
  double    p_death;
  double    p_translation;
  double    p_dilation;
  double    p_rotation;
  double    p_switching;
    // unfixed
  double    alm;  // RJMCMC
  double    b_zero; // multiple birth and death
  double    betampp;
  double    vk;   // RJMCMC
  double    delta;  // multiple birth and death
  double    alpha;
  int       nd_type_num;
  double    amp_th;
    // etc
  double    mean[MAX_CLASSES];
  double    vari[MAX_CLASSES];
  double    **blur;
  int       blur_size;
  double    sigma;

  int     optimization_type;

} MPP_Parameters;

/* Function to Parse Input Parameters*/
MPP_Parameters parse_input_parameters(int argc,char** argv);
void print_help(void);
void _print_current_parameters(MPP_Parameters mpp);

#endif
