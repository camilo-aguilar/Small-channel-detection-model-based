#ifndef _DEBUG
	#undef MPP_PARAMETER_TEST
#else
	#undef MPP_PARAMETER_TEST
#endif
#ifdef MPP_PARAMETER_TEST
	#define SHOW_WINDOW	0
	#define LEVEL0		0
	#define LEVEL1		0
	#define LEVEL2		0
	#define LEVEL3		0
	#define SHOW_TEXT	0
	#define DEBUG_TEXT	0
#else
#if 0
	#define SHOW_WINDOW	0
	#define LEVEL0		0
	#define LEVEL1		0
	#define LEVEL2		0
	#define LEVEL3		0
	#define SHOW_TEXT	0
	#define DEBUG_TEXT	0
#else
	#define SHOW_WINDOW	1
	#define LEVEL0		1
	#define LEVEL1		1
	#define LEVEL2		1
	#define LEVEL3		0
	#define SHOW_TEXT	1
	#define DEBUG_TEXT	0
#endif
#endif

