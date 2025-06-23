#include "dkiss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>

/******************************************************************/
/*  variables for the kiss random number generator                */
/******************************************************************/
unsigned int z, c;
unsigned long long x, y,t; 
/******************************************************************/
/*  kickstart the kiss random number generator                    */
/******************************************************************/
void start_kiss_pID(int seed1, int seed2)
{
	x = 123456789123ULL,y = 987654321987ULL; /* Seed variables */
	z = seed1,c =seed2; /* Seed variables */
	x = 1490024343005336237ULL * x + 123456789;
	y ^= y << 21; y ^= y >> 17; y ^= y << 30; /* Do not set y=0! */
	
	t = 4294584393ULL * z + c; c = t >> 32; z = t; /* Avoid z=c=0! */
}

void start_kiss(int seed)
{
	x = 123456789123ULL,y = 987654321987ULL; /* Seed variables */
	z = seed,c =6543217; /* Seed variables */
	x = 1490024343005336237ULL * x + 123456789;
	y ^= y << 21; y ^= y >> 17; y ^= y << 30; /* Do not set y=0! */
	
	t = 4294584393ULL * z + c; c = t >> 32; z = t; /* Avoid z=c=0! */
}

unsigned int kiss(void)
{
x = 1490024343005336237ULL * x + 123456789;
y ^= y << 21; y ^= y >> 17; y ^= y << 30; /* Do not set y=0! */
t = 4294584393ULL * z + c; c = t >> 32; z = t; /* Avoid z=c=0! */
return (unsigned int)(x>>32) + (unsigned int)y + z; /* Return 32-bit result */
}
/******************************************************************/
double dkiss(void)
{
    return ((double)kiss()+0.5)/4294967296.0;
}
/******************************************************************/


