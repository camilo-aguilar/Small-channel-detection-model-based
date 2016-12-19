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

#endif