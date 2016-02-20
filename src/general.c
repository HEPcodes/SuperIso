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

float Li3(float x)
/* calculates the trilogarithm function of x */
{
	float pisq6=16.*pow(atan(1.),2.)/6.;
	float x_0 = -1.;
	float x_1 = -0.85;
	float x_2 = 0.25;
	float x_3 = 0.63;
	float x_4 =  1.;
	if (x == 1.) return zeta3;
	if (x == -1.) return - 0.75 * zeta3;
	if (x <= x_0)
	{ 
		float lnx = log(-x);
		return Li3(1./x) - pisq6*lnx - lnx*lnx*lnx/6.;
	}
	else if (x < x_1)
	{
		return Li3(x*x)/4. - Li3(-x);
	}
	else if (x < x_2)
	{
		float z = - log(1.-x);
		float temp = z*(1.-3.*z/8.*(1.-17.*z/81.*(1.-15.*z/136.
                    *(1.-28.*z/1875.*(1.+5.*z/8.*(1.-304.*z/7203.
                    *(1.+945.*z/2432.*(1.-44.*z/675.*(1.+7.*z/24.
                    *(1.-26104.*z/307461.*(1.+1925.*z/8023.
                    *(1.-53598548.*z/524808375.
                    *(1.+22232925.*z/107197096.
                     )))))))))))));
		return temp; 
	}
	else if (x < x_3)
	{
		return Li3(x*x)/4. - Li3(-x);
	}
	else if (x < x_4)
	{
		float ln1x = log(1.-x); 
		return -Li3(1.-x) - Li3(-x/(1.-x)) + zeta3 + pisq6*ln1x - log(x)*ln1x*ln1x/2. + ln1x*ln1x*ln1x/6.; 
	}
	else 
	{ 
		float lnx = log(x);
		return Li3(1./x) + 2.*pisq6*lnx - lnx*lnx*lnx/6.;
	}
}

/*--------------------------------------------------------------------*/

complex float CLi2(complex float x)
/* calculates the dilogarithm function of x, extended to complex numbers */
{
	float pisq6=pow((4.*atan(1.)),2.)/6.;

	float x_0 = -0.30;
	float x_1 = 0.25;
	float x_2 = 0.51;
	if (x == 1.) return pisq6;
	if (creal(x) >= x_2) 
	{ 
		return pisq6 - CLi2(1.-x) - clog(x)*clog(1.-x);
	}
	if ((fabs(cimag(x)) > 1.) || (creal(x)*creal(x) + cimag(x)*cimag(x) > 1.2))
	{
		return - CLi2(1./x) - 0.5 * clog(-x) * clog(-x) - pisq6;
	}
	if (creal(x) <= x_0)
	{ 
		complex float zz = clog(1.-x);
   		return -CLi2(-x/(1.-x)) - zz*zz/2. ;
	}
 	else if (creal(x) < x_1)
	{
		complex float z = - clog(1.-x);
		complex float temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.
		*(1.-5.*z*z/294.*(1.-7.*z*z/360.
		*(1.-5.*z*z/242.*(1.-7601.*z*z/354900.
		*(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)))))))));
   		return temp;
	}
   	else return - CLi2(-x) + CLi2(x*x)/2.;
}

/*--------------------------------------------------------------------*/

float Cl2(float x)
/* calculates the Cl2 function of x */
{
	return cimag(CLi2(cos(x)+I*sin(x)));
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
