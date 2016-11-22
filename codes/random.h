/* random.h */

/* Created by Joel Dumke on 9/4/06
	This is just a header file for random.c. It was created for 
	compatibility with VS .NET */

#ifndef _RANDOM_
#define	_RANDOM_

double random2();
int random3();
void readseed();
void writeseed();
double normal();
double dexprand();
void srandom2(unsigned int num);

#endif
