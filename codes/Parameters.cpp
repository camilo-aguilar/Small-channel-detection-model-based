/* 
 File "parameters.cpp"
 Added by Camilo Aguilar
 Created by Thomas Cool

 v1: Dec 06 2016

*/

#include "parameters.h"
 

INPUT_PARAMS parse_input_qc(int argc,char** argv)
{
	INPUT_PARAMS inp; 
	

	inp.iterations = 3000000; //	inp.iterations = 3000000;
  	inp.gamma_d = 1000; //inp.gamma_d/
  	inp.w_eo = 8;//1 inp.w_eo
  	inp.w_f = 6;//1 inp.w_f
  	inp.w_s = -10;//1 inp.w_s
  	inp.w_d = -7;//1 inp.w_d
  	inp.w_io = 1;//1 inp.w_io
	



	if(argc>2)
	{
	

		int cntr=2;
		while(cntr<argc)
		{
			if(strcmp(argv[cntr],"-n")==0)
			{
				inp.iterations=atoi(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-gd")==0)
			{
				inp.gamma_d=atof(argv[++cntr]);
			}
			else if(strcmp(argv[cntr],"-weo")==0)
			{
				inp.w_eo=atof(argv[++cntr]);
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
			else if(strcmp(argv[cntr],"-wio")==0)
			{
				inp.w_io=atof(argv[++cntr]);
			}
			else
			{
				fprintf(stderr, "%s is not a recognized argument\n",argv[cntr]);
                exit(100);

			}
			cntr++;
		}
	}

	
  return inp;
}
