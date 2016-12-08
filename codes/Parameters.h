#ifndef _PARAMETERS_
#define _PARAMETERS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

 
/*Input Parameters for Candy Model*/
typedef struct _input_params 
{
  int iterations;
  int gamma_d;//34;
  int w_eo;//1
  int w_f;//1
  int w_s;//1
  int w_d;//1
  int w_io;//1
} INPUT_PARAMS;


/* Function to Parse Input Parameters*/
INPUT_PARAMS parse_input_qc(int argc,char** argv);

#endif