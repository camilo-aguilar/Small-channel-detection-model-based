/* 
 File "parameters.cpp"
 Added by Camilo Aguilar
 Created by Thomas Cool

 v1: Dec 06 2016

*/

#include "Parameters.h"

 
MPP_Parameters parse_input_parameters(int argc,char** argv)
{

	int   PARAM_OPTIMIZATION_TYPE = 2;	    			//3  0:RJMCMC, 1:Multiple Birth and Death	2:RJMCMC Quality Candy model
	float PARAM_LAMBDA_RJMCMC = 30000;	     			//5  Lambda(RJMCMC)  (intensity = lambda*vk), alm(MBND)0.5
	float PARAM__BETA_MPP = 14;		      				//7  Betampp(1)
	float PARAM_LAMBDA_l = 0;	            			//8	 Lambda_l		
	float PARAM_LAMBDA_int = 0.005;		      			//9  Lambda_int
	int   PARAM_CHANNEL_TYPES = 1;               		//10 0:simple channel	1:necking only	2: denting only 3:necking and denting
	float PARAM_AMPLITUDE_THRESHOLD = 21.0;	      		//12 Amp_th 	
	float PARAM_GAUSSIAN_TAU = 17;	

	MPP_Parameters inp; 
	

	inp.iter_num = 3000000; //	inp.iterations = 3000000;
  	inp.gamma_d = 1000; //inp.gamma_d/
  	inp.w_eo = 8;//1 inp.w_eo
  	inp.w_f = 6;//1 inp.w_f
  	inp.w_s = -10;//1 inp.w_s
  	inp.w_d = -7;//1 inp.w_d
  	inp.w_io = 1;//1 inp.w_io
	
	//Multiple Birth and Death
	inp.p_birth			= 0.3;
	inp.p_death			= 0.3;
	inp.p_translation	= 0.2;
	inp.p_dilation		= 0.2;
	inp.p_rotation		= 0.08;
	inp.p_switching		= 0.;

	inp.hard_repulsion	= 4;	
	inp.gaussian_tau	= 15;	
	//Denting and Necking Channels
	inp.symmetry_th		= 0;	// 8
	inp.lambda_s		= 0;	//Symmetry Potential
	inp.lambda_nc		= 0;
	

	//Denting
	inp.dent_l_w_ratio	= 1.3; // 1.25 (MBND) 1.3(RJMCMC)1.3
	

 
 
	inp.alpha			= 1;
	inp.sigma			= 0.9;
	inp.blur_size		= 7;
	
	//Very Important
	inp.lambda_e = 1;			//Object Potential
	inp.lambda_dc= 3;			//Discontinuity Potential 

	
	inp.lambda_a = 0.5;
	inp.length_th = inp.lengthmin;


	inp.T0 				= INITIAL_T0;
	inp.iter_num 		= ITERATIONS;
	inp.gamma_d  		= GAMMA_D;	
	inp.w_f 			= W_F;
	inp.w_s 			= W_S;
	inp.w_d 			= W_D;
	inp.w_io 			= W_IO;
	inp.w_eo 			= W_EO;
	inp.error_th 		= ERROR_TH;
	inp.widthmin		= W_MIN;	
	inp.widthmax		= W_MAX;	
	inp.lengthmin		= L_MIN;	
	inp.lengthmax		= L_MAX;	
	inp.thetamin 		=THETA_MIN;
	inp.thetamax 		=THETA_MAX;
	inp.de_coeff 		= DECREASE_COEFFICIENT;

	inp.optimization_type = PARAM_OPTIMIZATION_TYPE;	
	inp.alm = PARAM_LAMBDA_RJMCMC;
	inp.b_zero = inp.alm;
	inp.betampp = PARAM__BETA_MPP;
	inp.lambda_l = PARAM_LAMBDA_l;	// 0.12
	inp.lambda_int = PARAM_LAMBDA_int;
	inp.nd_type_num = PARAM_CHANNEL_TYPES;	// necking only = 1, 
	
	inp.amp_th = PARAM_AMPLITUDE_THRESHOLD;
	inp.gaussian_tau = PARAM_GAUSSIAN_TAU;


	inp.alm = LAMBDA_RJMCMC;
	inp.b_zero = inp.alm;
	inp.lambda_l = LAMBDA_L;	// 0.12
	inp.lambda_int = LAMBDA_INT;
	inp.nd_type_num = CHANNEL_TYPES;	// necking only = 1, 
	inp.gaussian_tau = GAUSSIAN_TAU;
	inp.betampp = BETA_MPP;

	
	
	inp = _parse_hard_parameters(inp);
	
	if(argc>2)
	{
	

		int cntr=2;
		while(cntr<argc)
		{
			if(strcmp(argv[cntr],"-n")==0)
			{

				inp.iter_num=atoi(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-gd")==0)
			{
				inp.gamma_d=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-wf")==0)
			{
				inp.w_f=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-ws")==0)
			{
				inp.w_s=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-wd")==0)
			{
				inp.w_d=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-weo")==0)
			{
				inp.w_eo=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-wio")==0)
			{
				inp.w_io=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr], "-t0")==0)
			{
				inp.T0 = atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr], "-de_coeff")==0)
			{
				inp.de_coeff = atof(argv[++cntr]);
			}

			else if(strcmp(argv[cntr], "-e_th")==0)
			{
				inp.error_th = atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr], "-h")==0)
			{
				print_help();
				exit(0);
			}
			else
			{
				fprintf(stderr, "\n\n\nERROR: %s is not a recognized argument\n\n\n",argv[cntr]);
				print_help();
                exit(100);

			}
			cntr++;
		}
	}
	
  return inp;
}

MPP_Parameters _parse_hard_parameters(MPP_Parameters input_p)
{
	FILE *f = fopen("hard_parameters.txt", "r");
	if(f == NULL)
	{
		printf("No hard parameters file, using default hard_parameters \n");
		return input_p;
	}

	char tmp_str[30];
	MPP_Parameters output_p = input_p;	
	int number_of_parameters = 0;
	int i = 0;

	fscanf(f, "%d", &number_of_parameters);
 	for(i=0; i< number_of_parameters; i++)
 	{
 		fscanf(f, "%s : %d",tmp_str , &output_p.optimization_type);
 		fscanf(f, "%s : %lf",tmp_str , &output_p.lengthmin);
 		fscanf(f, "%s : %lf",tmp_str , &output_p.lengthmax);
 		fscanf(f, "%s : %lf",tmp_str , &output_p.widthmin);
 		fscanf(f, "%s : %lf",tmp_str , &output_p.widthmax);
 		fscanf(f, "%s : %lf",tmp_str , &output_p.thetamin);
 		fscanf(f, "%s : %lf",tmp_str , &output_p.thetamax);

 		
 	}	

 	return output_p;
	
}

void print_help(void)
{
	printf("*******************************************************\n");
	printf("USAGE:\n");
	printf("ChannelMpp [name (w/o extension)] [OPTIONAL PARAMETERS]\n");
	printf("Example: \n");
	printf("ChannelMpp cracked -n 10000 -ws -20 -wd -30 ...\n");
	printf("PARAMETERS [and recommendations]:\n");
	printf("Iterations              :		-n [1 <->10000000]\n");
	printf("Data Potential Energy   :  	 	-gd [-50 <-> 0]\n");
	printf("Free Segment Weight     :  	 	-wf [-50 <-> 50]\n");
	printf("Single Segment Weight   :  	 	-ws [-50 <-> 50]\n");
	printf("Double Segment Weight   :  	 	-wd [-50 <-> 50]\n");
	printf("External Bad Orientation:		-weo [-50 <-> 50]\n");
	printf("Internal Bad Orientation:		-wio [-50 <-> 50]\n");
	
	printf("Initial Temperature 	: 		-t0 [~t_test(background,foreground)]\n");
	printf("Decrease Coefficient 	: 		-de_coeff [~t_test(background,foreground)]\n");
	printf("Error Threshold 		: 		-e_th [~t_test(background,foreground)]\n");
	printf("*******************************************************\n");


}




void _print_current_parameters(MPP_Parameters mpp)
{


	printf("Iterations:              		%d\n",mpp.iter_num);
	printf("Data Potential Energy:  	 	%f\n", mpp.gamma_d);
	printf("Free Segment Weight:    	 	%f\n",mpp.w_f);
	printf("Single Segment Weight:  	 	%f\n",mpp.w_s);
	printf("Double Segment Weight:  	 	%f\n",mpp.w_d);
	printf("External Bad Orientation:		%f\n",mpp.w_eo);
	printf("Internal Bad Orientation:		%f\n",mpp.w_io);

	printf("Initial Temperature: 	 		%f\n",mpp.T0);
	printf("Decrease Coefficient: 	 		%f\n",mpp.de_coeff);
	printf("Error_Threshold:				%f\n",mpp.error_th);

	printf("L_Min: 							%f\n", mpp.lengthmin);
	printf("L_Max: 							%f\n", mpp.lengthmax);
	printf("W_Min: 							%f\n", mpp.widthmin);
	printf("W_Max: 							%f\n", mpp.widthmax);
	printf("Theta_Min: 						%f\n", mpp.thetamin);
	printf("Theta_Max: 						%f\n", mpp.thetamax);
}


