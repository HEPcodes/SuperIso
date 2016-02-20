#include "include.h"


double Li2(double x)
/* calculates the dilogarithm function of x */
{
	double pisq6=pow((4.*atan(1.)),2.)/6.;
	double x_0 = -0.3;
	double x_1 = 0.25;
	double x_2 = 0.51;
	if (x == 1.) return pisq6;
	if (x <= x_0)
	{ 
		double temp = log(fabs(1.-x));
		return -Li2(-x/(1.-x)) - temp*temp/2. ; 
	}
	else if (x < x_1)
	{
		double z = - log(1.-x);
		double temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.*(1.-5.*z*z/294.*(1.-7.*z*z/360.*(1.-5.*z*z/242.*(1.-7601.*z*z/354900.*(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)))))))));
		return temp; 
	}
	else if (x < x_2) return - Li2(-x) + Li2(x*x)/2. ;
	else 
	{ 
		return pisq6 - Li2(1.-x) - log(fabs(x))*log(fabs(1.-x)) ; 
	}
}

/*--------------------------------------------------------------------*/

double Li3(double x)
/* calculates the trilogarithm function of x */
{
	double pisq6=16.*pow(atan(1.),2.)/6.;
	double x_0 = -1.;
	double x_1 = -0.85;
	double x_2 = 0.25;
	double x_3 = 0.63;
	double x_4 =  1.;
	if (x == 1.) return zeta3;
	if (x == -1.) return - 0.75 * zeta3;
	if (x <= x_0)
	{ 
		double lnx = log(-x);
		return Li3(1./x) - pisq6*lnx - lnx*lnx*lnx/6.;
	}
	else if (x < x_1)
	{
		return Li3(x*x)/4. - Li3(-x);
	}
	else if (x < x_2)
	{
		double z = - log(1.-x);
		double temp = z*(1.-3.*z/8.*(1.-17.*z/81.*(1.-15.*z/136.
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
		double ln1x = log(1.-x); 
		return -Li3(1.-x) - Li3(-x/(1.-x)) + zeta3 + pisq6*ln1x - log(x)*ln1x*ln1x/2. + ln1x*ln1x*ln1x/6.; 
	}
	else 
	{ 
		double lnx = log(x);
		return Li3(1./x) + 2.*pisq6*lnx - lnx*lnx*lnx/6.;
	}
}

/*--------------------------------------------------------------------*/

complex double CLi2(complex double x)
/* calculates the dilogarithm function of x, extended to complex numbers */
{
	double pisq6=pow((4.*atan(1.)),2.)/6.;

	double x_0 = -0.30;
	double x_1 = 0.25;
	double x_2 = 0.51;
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
		complex double zz = clog(1.-x);
   		return -CLi2(-x/(1.-x)) - zz*zz/2. ;
	}
 	else if (creal(x) < x_1)
	{
		complex double z = - clog(1.-x);
		complex double temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.
		*(1.-5.*z*z/294.*(1.-7.*z*z/360.
		*(1.-5.*z*z/242.*(1.-7601.*z*z/354900.
		*(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)))))))));
   		return temp;
	}
   	else return - CLi2(-x) + CLi2(x*x)/2.;
}

/*--------------------------------------------------------------------*/

double Cl2(double x)
/* calculates the Cl2 function of x */
{
	return cimag(CLi2(cos(x)+I*sin(x)));
}

/*--------------------------------------------------------------------*/

double max(double x, double y)
{
	if(x<y) return y; else return x;
}

/*--------------------------------------------------------------------*/

double min(double x, double y)
{
	if(x<y) return x; else return y;
}

/*-----------------------------------------------------------------------*/
/* Modified Bessel functions of second type and of order 0,1,2 */
/*-----------------------------------------------------------------------*/

double I0(double x)
{
	double y;

	double absx=fabs(x);

	if (absx < 3.75) 
	{
		y=x/3.75;
		y=y*y;
		return 1.+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
 	} 
	else
	{
		y=3.75/absx;
		return (exp(absx)/sqrt(absx))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
 	}
}

double I1(double x)
{
	double y,tmp;
	double I1=0.;

	double absx=fabs(x);

	if (absx < 3.75) 
	{
		y=x/3.75;
		y=y*y;
      		I1 = absx*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	}
	else
	{
		y=3.75/absx;
     		tmp=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
      		tmp=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*tmp))));
      		I1*=(exp(absx)/sqrt(absx));
	}
	if(x<0.) return -I1; else return I1;
}

double K0(double x)
{
	double y;

	if (x <= 2.) 
	{
		y=x*x/4.;
		return (-log(x/2.)*I0(x))+(-0.57721566+y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5))))));
	} 
	else 
	{
		y=2./x;
      		return (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))));
   	}
}

double K1(double x)
{
 	double y;

	if (x <= 2.) 
	{
		y=x*x/4.;
		return (log(x/2.)*I1(x))+(1./x)*(1.+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	}
	else 
	{
		y=2./x;
		return (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
   	}
}

double K2(double x)
{
   return K0(x)+2./x*K1(x);
}


double K0exp(double x, double z)
{
	double y;

	if (x <= 2.) 
	{
		y=x*x/4.;
		return ((-log(x/2.)*I0(x))+(-0.57721566+y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5)))))))*exp(z)*sqrt(z);
	} 
	else 
	{
		y=2./x;
      		return (exp(z-x)/sqrt(x/z))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))));
   	}
}

double K1exp(double x,double z)
{
 	double y;

	if (x <= 2.) 
	{
		y=x*x/4.;
		return ((log(x/2.)*I1(x))+(1./x)*(1.+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4))))))))*exp(z)*sqrt(z);
	}
	else 
	{
		y=2./x;
		return (exp(z-x)/sqrt(x/z))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
   	}
}

double K2exp(double x,double z)
{
   return K0exp(x,z)+2./x*K1exp(x,z);
}

/*--------------------------------------------------------------------*/

int test_integer(char name[])
/* tests if the string "name" is an integer, and return 1 if so, 0 otherwise */
{
	char *testint;
	strtol(name,&testint,10);
	return (*testint=='\0');
}
