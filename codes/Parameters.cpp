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
	if(argc>2)
	{
	
		int remaining = argc -2;

		int cntr=2;
		while(cntr<argc)
		{
			if(strcmp(argv[cntr++],"-n")==0 && remaining--)
			{
				inp.iterations=atoi(argv[cntr++]);
			}
			else if(strcmp(argv[cntr++],"-gd")==0 && remaining--)
			{
				inp.gamma_d=atoi(argv[cntr++]);
			}
			else if(strcmp(argv[cntr++],"-weo")==0 && remaining--)
			{
				inp.w_eo=atoi(argv[cntr++]);
			}
			else if(strcmp(argv[cntr++],"-wf")==0 && remaining--)
			{
				inp.w_f=atoi(argv[cntr++]);
			}
			else if(strcmp(argv[cntr++],"-ws")==0 && remaining--)
			{
				inp.w_s=atoi(argv[cntr++]);
			}
			else if(strcmp(argv[cntr++],"-wd")==0 && remaining--)
			{
				inp.w_d=atoi(argv[cntr++]);
			}
			else if(strcmp(argv[cntr++],"-wio")==0 && remaining--)
			{
				inp.w_io=atoi(argv[cntr++]);
			}
		}
	}
	else
	{
		inp.iterations = 3000000;
  		inp.gamma_d= 34;//34;
  		inp.w_eo= 8;//1
  		inp.w_f= 6;//1
  		inp.w_s= -10;//1
  		inp.w_d= -7;//1
  		inp.w_io= 1;//1
	}
	
  return inp;
}
