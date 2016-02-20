#include "include.h"

/*---------------------------------------------------------------------*/

float f(float x)
{
	float f;
	f=1.-8.*x+8.*x*x*x-x*x*x*x-12.*x*x*log(x);
	return f;
}
/*---------------------------------------------------------------------*/
	
float f77(float x)
{
	float f77;
	f77=(10.+x-2.*x*x/3.+(x-4.)*log(x))*x/3.;
	return f77;
}
/*---------------------------------------------------------------------*/

float f78(float x)
{
	float f78;
	float pi=4.*atan(1.);
	f78=8.*(Li2(1.-x)-pi*pi/6.-x*log(x)+9.*x/4.-x*x/4.+x*x*x/12.)/9.;
	return f78;
}
/*---------------------------------------------------------------------*/

float f88(float x, float b)
{
	float f88;
	float pi=4.*atan(1.);
	f88=(4.*Li2(1.-x)-2.*pi*pi/3.+8.*log(1.-x)-x*(2.+x)*log(x)+
	+7.*x+3.*x*x-2.*x*x*x/3.-2.*(2.*x+x*x+4.*log(1.-x))*log(b))/27.;
	return f88;
}
/*---------------------------------------------------------------------*/

float f22(float delta, float z)
{
	float pi=4.*atan(1.);
	int ie;

	int n1=500.;
	int n2=500.;

	float x=0.;
	float t;

	float int1=0;
	float int2=0.;

	for(ie=1;ie<=n1-1;ie++)
	{
		x+=4.*z/n1;
		t=x/z;
		int1 +=(1-x)*(1.-max(x,1.-delta))*pow((-2./t*pow(atan(sqrt(t/(4.-t))),2.)+1./2.),2.);
	}
	x=4.*z;
	int1 +=1./2.*(1-x)*(1.-max(x,1.-delta))*pow(-pi*pi/8.+1./2.,2.);
	int1 *= (4.*z)/n1;

	for(ie=1;ie<=n2-1;ie++)
	{
		x+=(1.-4.*z)/n2;
		t=x/z;
		int2 +=(1-x)*(1-max(x,1.-delta))*(pow(2./t*(pow(log(sqrt(t)/2.+sqrt(t-4.)/2.),2.)-pi*pi/4.)+1./2.,2.)+pow(2.*pi/t*log(sqrt(t)/2.+sqrt(t-4.)/2.),2.));
	}
	x=4.*z;
	int2 +=1./2.*(1.-x)*(1.-max(x,1.-delta))*pow(-pi*pi/8.+1./2.,2.);
	int2 *= (1.-4.*z)/n2;

	float f22=16./27.*(int1+int2);
	return f22;
}
/*---------------------------------------------------------------------*/

float f27(float delta, float z)
{
	/* z < 1/4 */
	float pi=4.*atan(1.);
	int ie;

	int n1=500.;
	int n2=500.;

	float x=0.;
	float t;

	float int1=0;
	float int2=0.;

	for(ie=1;ie<=n1-1;ie++)
	{
		x+=4*z/n1;
		t=x/z;
		int1 +=(1.-max(x,1.-delta))*(-2.*pow(atan(sqrt(t/(4.-t))),2.)+t/2.);
	}
	x=4.*z;
	int1 +=1./2.*(1.-max(x,1.-delta))*pow(-pi*pi/2.+2.,2.);
	int1 *= (4.*z)/n1;

	for(ie=1;ie<=n2-1;ie++)
	{
		x+=(1.-4.*z)/n2;
		t=x/z;
		int2 +=(1.-max(x,1.-delta))*(2.*(pow(log(sqrt(t)/2.+sqrt(t-4.)/2.),2.)-pi*pi/4.)+t/2.);
	}
	x=4.*z;
	int2 +=1./2.*(1.-max(x,1.-delta))*pow(-pi*pi/2.+2.,2.);
	int2 *= (1.-4.*z)/n2;

	float f27=-8./9.*z*(int1+int2);
	return f27;
}

float BR_bsgamma(float C0[], float C1[], float Cem[], float mu, struct parameters* param)
/* computes the inclusive branching ratio of B -> Xs gamma */
/* C0 and C1: respectively LO and NLO contributions of the Wilson coefficients at scale mu=O(mb) */
{

	float pi=4.*atan(1.);
	float alpha=1./137.;

	float alpha_s_mb=alpha_s_running(param->mass_b,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);	

 	float mass_b_pole=param->mass_b*(1+alpha_s_mb/pi*(4./3.+alpha_s_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*((param->mass_u+param->mass_d+param->mass_s+param->mass_c)/param->mass_b))+alpha_s_mb/pi*(190.595-4.*(26.655-4.*0.6527)))));

	float alpha_s_mu=alpha_s_running(mu,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);	

	float z0=pow(param->mass_c/param->mass_b,2.); 	  
	float z1=pow(param->mass_c/mass_b_pole,2.);   
	
	float mbs=mass_b_pole/param->mass_s;	
	 
	float lambda2=0.12;      /* lambda2=(m_Bstar^2-m_B^2)/4  hadronic parameter */
	float delta=0.9;	
	
	float Vts=-0.0405;
	float Vtb=0.99915;
	float Vcb=0.041;
		
	float r7=-10./3.-8.*pi*pi/9.;
	float Re_r8=44./9.-8.*pi*pi/27.;
	float Re_r2=-4.092+12.78*(sqrt(z1)-.29);

	float kabar=3.382-4.14*(sqrt(z1)-.29);


	float gamma77=32./3.;
	float gamma27=416./81.;
	float gamma87=-32./9.;
	
	float S_delta=exp(-2.*alpha_s_mu*log(delta)*(log(delta)+7./2.)/3./pi);

	float k77 = S_delta*(1+alpha_s_mu*(r7+gamma77*log(mass_b_pole/mu)-16./3.)/2./pi
+6.*(pow((1.-z0),4.)/f(z0)-1.)*lambda2/param->mass_b/param->mass_b)+alpha_s_mu*f77(delta)/pi+S_delta*alpha_s_mu*kabar/pi/2.;

	float k27=S_delta*(alpha_s_mu*(Re_r2+gamma27*log(mass_b_pole/mu))/2./pi-lambda2/9./param->mass_c/param->mass_c)+alpha_s_mu*(f27(delta,z1))/pi;

	float k78=S_delta*alpha_s_mu*(Re_r8+gamma87*log(mass_b_pole/mu))/2./pi+alpha_s_mu*f78(delta)/pi;
		
	float k22=alpha_s_mu*f22(delta,z1)/pi;
	float k88=alpha_s_mu*f88(delta,mbs)/pi;
	float k28=alpha_s_mu*(-f27(delta,z1))/3./pi;
	
	float kem_SL=2.*alpha_s_mu*log(param->mass_W/mu)/pi;
		
	float Rth=6.*0.105*alpha*pow(Vts*Vtb/Vcb,2.)/pi/f(z0);
	
	float K_NLO= k22*C0[1]*C0[1]+ k77*C0[7]*C0[7]+ k88*C0[8]*C0[8]+ k27*C0[1]*C0[7]+ k28*C0[1]*C0[8]+k78*C0[7]*C0[8] + S_delta*alpha_s_mu*C1[7]*C0[7]/2./pi+ S_delta*alpha*(2.*Cem[7]*C0[7]-C0[7]*C0[7]*kem_SL)/alpha_s_mu;

	float BRinc=Rth*K_NLO;
	return BRinc;
}
/*---------------------------------------------------------------------*/

float BRbsgamma_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the inclusive branching ratio of b -> s gamma */
{
	float C0[9],C1[9],Cem[9];
	struct parameters param;
	
	float pi=4.*atan(1.);
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	float alpha_s_mb=alpha_s_running(param.mass_b,param.mass_top_pole,param.mass_b,param.alpha_s_MZ,param.mass_Z);	

 	float mass_b_pole=param.mass_b*(1+alpha_s_mb/pi*(4./3.+alpha_s_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*((param.mass_u+param.mass_d+param.mass_s+param.mass_c)/param.mass_b))+alpha_s_mb/pi*(190.595-4.*(26.655-4.*0.6527)))));

	float mu_b=mass_b_pole;

	C_calculator(C0,C1, mu_b, &param);
	Cem_calculator(Cem, mu_b, &param);
	
	return BR_bsgamma(C0,C1,Cem,mu_b,&param);
}

