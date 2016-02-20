#include "include.h"

float Li2(float x)
/* calculates the dilogarithm function of x */
{
	float pisq6=pow((4.*atan(1.)),2.)/6.;
	float x_0 = -0.3;
	float x_1 = 0.25;
	float x_2 = 0.51;
	if (x == 1.) return pisq6;
	if (x <= x_0)
	{ 
		float temp = log(fabs(1.-x));
		return -Li2(-x/(1.-x)) - temp*temp/2. ; 
	}
	else if (x < x_1)
	{
		float z = - log(1.-x);
		float temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.*(1.-5.*z*z/294.*(1.-7.*z*z/360.*(1.-5.*z*z/242.*(1.-7601.*z*z/354900.*(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)))))))));
		return temp; 
	}
	else if (x < x_2) return - Li2(-x) + Li2(x*x)/2. ;
	else 
	{ 
		return pisq6 - Li2(1.-x) - log(fabs(x))*log(fabs(1.-x)) ; 
	}
}

/*--------------------------------------------------------------------*/

float max(float x, float y)
{
	if(x<y) return y; else return x;
}

/*--------------------------------------------------------------------*/

float min(float x, float y)
{
	if(x<y) return x; else return y;
}
