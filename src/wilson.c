#include "include.h"


float E(float x)
{
	float E;
	E = x*(-18.+11.*x+x*x)/12./pow(x-1.,3.)+x*x*(15.-16.*x+4.*x*x)/6./pow(x-1.,4.)*log(x) -2./3.*log(x);
	return E;
}

/*----------------------------------------------------------------------*/

float F4H(float x)
{
	float F4H;
	F4H = x*(16.-29.*x+7.*x*x)/36./pow(x-1.,3.)+x*(3.*x-2.)/6./pow(x-1.,4.)*log(x);
	return F4H;
}

/*----------------------------------------------------------------------*/

float F4ch(float x)
{
	float F4ch;
	F4ch = x*(11.-7.*x+2.*x*x)/18./pow(x-1.,3.)-x/3./pow(x-1.,4.)*log(fabs(x));
	return F4ch;
}

/*----------------------------------------------------------------------*/

float F7_1(float x)
{
	float F7_1;
	F7_1 = x*(7.-5.*x-8.*x*x)/24./pow(x-1.,3.)+x*x*(3.*x-2.)/4./pow(x-1.,4.)*log(x);
	return F7_1;
}

/*----------------------------------------------------------------------*/

float F7_2(float x)
{
	float F7_2;
	F7_2 = x*(3.-5.*x)/12./pow(x-1.,2.)+x*(3.*x-2.)/6./pow(x-1.,3.)*log(x);
	return F7_2;
}

/*----------------------------------------------------------------------*/

float F7_3(float x)
{
	float F7_3;
	F7_3 = (5.-7.*x)/6./pow(x-1.,2.)+x*(3.*x-2.)/3./pow(x-1.,3.)*log(x);
	return F7_3;
}

/*----------------------------------------------------------------------*/

float F8_1(float x)
{
	float F8_1;
	F8_1 = x*(2.+5.*x-x*x)/8./pow(x-1.,3.)-3.*x*x/4./pow(x-1.,4.)*log(x);
	return F8_1;
}

/*----------------------------------------------------------------------*/

float F8_2(float x)
{
	float F8_2;
	F8_2 = x*(3.-x)/4./pow(x-1.,2.)-x/2./pow(x-1.,3.)*log(x);
	return F8_2;
}

/*----------------------------------------------------------------------*/

float F8_3(float x)
{
	float F8_3;
	F8_3 = (1.+x)/2./pow(x-1.,2.)-x/pow(x-1.,3.)*log(x);
	return F8_3;
}

/*----------------------------------------------------------------------*/

float G7(float x)
{
	float G7;
	G7 = (-16.*pow(x,4.)-122.*pow(x,3.)+80.*x*x-8.*x)/9./pow(x-1.,4.)*Li2(1.-1./x)+(6.*pow(x,4.)+46.*pow(x,3.)-28.*x*x)/3./pow(x-1.,5.)*pow(log(x),2.)+(-102.*pow(x,5.)-588.*pow(x,4.)-2262.*pow(x,3.)+3244.*x*x-1364.*x+208.)/81./pow(x-1.,5.)*log(x)+(1646.*pow(x,4.)+12205.*pow(x,3.)-10740.*x*x+2509.*x-436.)/486./pow(x-1.,4.);
	return G7;
}

/*----------------------------------------------------------------------*/

float G8(float x)
{
	float G8;
	G8 = (-4.*pow(x,4.)+40.*pow(x,3.)+41.*x*x+x)/6./pow(x-1.,4.)*Li2(1.-1./x)+(-17.*pow(x,3.)-31.*x*x)/2./pow(x-1.,5.)*pow(log(x),2.)+(-210.*pow(x,5.)+1086.*pow(x,4.)+4893.*pow(x,3.)+2857.*x*x-1994.*x+280.)/216./pow(x-1.,5.)*log(x)+(737.*pow(x,4.)-14102.*pow(x,3.)-28209.*x*x+610.*x-508.)/1296./pow(x-1.,4.);
	return G8;
}

/*----------------------------------------------------------------------*/

float H2(float x, float y)
{
	float H2;

	if(fabs(x-y)<1.e-5)
	{
		H2=1./pow(x-1.,2.)*(1.-x+log(x));
	}
	else
	{
		H2=x/(1.-x)/(x-y)*log(x)+y/(1.-y)/(y-x)*log(y);
	}
	return H2;
}

/*----------------------------------------------------------------------*/

float G7H(float x, float tb)
{
	float G7H;
	
	float G71=4.*(-3.+7.*x-2.*x*x)*Li2(1.-1./x)+(8.-14.*x-3.*x*x)*log(x)*log(x)/(x-1.)
	+2.*(-3.-x+12.*x*x-2*x*x*x)*log(x)/(x-1.)+3.*(7.-13.*x+2*x*x);
	float G72=x*(18.-37.*x+8.*x*x)*Li2(1.-1./x)+x*(-14.+23.*x+3.*x*x)*log(x)*log(x)/(x-1.)+
	(-50.+251.*x-174.*x*x-192.*x*x*x+21.*x*x*x*x)*log(x)/9./(x-1.)+(797.-5436.*x+7569.*x*x-1202*x*x*x)/108.;
	G7H=-4.*x*G71/9./pow((x-1.),3.)+2.*x*G72/pow((x-1.),4.)/9./tb/tb;

	return G7H;
}

/*----------------------------------------------------------------------*/
	
float Delta7H(float x, float tb)
{
	float Delta7H;
	float Delta71=(x-1.)*(21.-47.*x+8.*x*x)+2.*(-8.+14.*x+3.*x*x)*log(x);
	float Delta72=(x-1.)*(-31.-18.*x+135.*x*x-14.*x*x*x)/6.+x*(14.-23.*x-3.*x*x)*log(x);
	Delta7H=2.*x/9./pow((x-1.),4.)*(-Delta71+Delta72/(x-1.)/tb/tb);
	
	return Delta7H;
}

/*----------------------------------------------------------------------*/

float EH(float x, float tb)
{
	float EH;

	EH=x*((x-1.)*(16.-29.*x+7.*x*x)+6.*(3.*x-2.)*log(x))/36./pow((x-1.),4.)/tb/tb;
	return EH;
}

/*----------------------------------------------------------------------*/

float G8H(float x, float tb)
{
	float G8H;

	float G81=0.5*(-36.+25.*x-17*x*x)*Li2(1.-1./x)+(19.+17.*x)*log(x)*log(x)/(x-1.)
	+0.25*(-3.-187.*x+12.*x*x-14.*x*x*x)*log(x)/(x-1.)+3.*(143.-44.*x+29.*x*x)/8.;
	float G82=x*(30.-17.*x+13.*x*x)*Li2(1.-1./x)-x*(31.+17.*x)*log(x)*log(x)/(x-1.)+
	(-226.+817.*x+1353.*x*x+318.*x*x*x+42.*x*x*x*x)*log(x)/36./(x-1.)+(1130.-18153.*x+7650.*x*x-4451.*x*x*x)/216.;	
	G8H=-x*G81/3./pow((x-1.),3.)+ x*G82/pow((x-1.),4.)/6./tb/tb;

	return G8H;
}

/*----------------------------------------------------------------------*/

float Delta8H(float x, float tb)
{
	float Delta8H;

	float Delta81=(x-1.)*(81.-16.*x+7.*x*x)-2.*(19.+17.*x)*log(x);
	float Delta82=(x-1.)*(-38.-261.*x+18.*x*x-7.*x*x*x)/6.+x*(31.+17.*x)*log(x);
	Delta8H=(x/6./pow((x-1.),4.))*(-Delta81+Delta82/(x-1.)/tb/tb);	

	return Delta8H;
}

/*----------------------------------------------------------------------*/

void C_calculator(float C0[], float C1[],float mu, struct parameters* param)
/* calculates the LO (C0) and NLO (C1) contributions to the Wilson coefficients at scale mu, using the parameters of the structure param */
{
	int ie;
	
	float pi=4.*atan(1.);

	for(ie=1;ie<=8;ie++) C0[ie]=C1[ie]=0.;

	float alpha_s_mb=alpha_s_running(param->mass_b,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);	

	float mass_b_pole=param->mass_b*(1+alpha_s_mb/pi*(4./3.+alpha_s_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*(1.61/4.62))+alpha_s_mb/pi*(190.595-4.*(26.655-4.*0.6527)))));
/* pole mass of b, from PDG */
	
	float alpha_s_MW=alpha_s_running(param->mass_W,param->mass_top_pole,mass_b_pole,param->alpha_s_MZ,param->mass_Z); 
	
	float alpha_s_mtop=alpha_s_running(param->mass_top_pole,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z); 
	
	float mtmt=param->mass_top_pole*(1.-4./3.*alpha_s_mtop/pi);
	
	float mass_top=running_mass(mtmt,mtmt,param->mass_W,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);
/* mt(MW) */

	float sw=sin(atan(param->gp/param->g2));

/*----------------------------------------------------------------------*/
	/* EPSILON_B */
/*----------------------------------------------------------------------*/
	
	float alpha_s_MSOFT=alpha_s_running(param->MSOFT_Q,param->mass_top_pole,mass_b_pole,param->alpha_s_MZ,param->mass_Z);
		
	float epsilonb=-2.0/3.0*alpha_s_MSOFT/pi*(param->mu_Q-param->A_b/param->tan_beta)/param->mass_gluino*H2(param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) -param->yut*param->yut/16.0/pi/pi*(param->A_t-param->mu_Q/param->tan_beta)*(param->charg_Umix[1][2]*param->charg_Vmix[1][2]/param->mass_cha1 *H2(param->mass_t1*param->mass_t1/param->mass_cha1/param->mass_cha1,param->mass_t2*param->mass_t2/param->mass_cha1/param->mass_cha1) +param->charg_Umix[2][2]*param->charg_Vmix[2][2]/param->mass_cha2*H2(param->mass_t1*param->mass_t1/param->mass_cha2/param->mass_cha2,param->mass_t2*param->mass_t2/param->mass_cha2/param->mass_cha2));

	
	if(fabs(param->mu_Q*param->mu_Q/param->M2_Q/param->M2_Q-1)<0.01)
	{
		epsilonb+=1./param->inv_alpha_em/sw/sw/4.0/pi*(-param->mu_Q*param->M2_Q)*(param->stop_mix[1][1]*param->stop_mix[1][1]*param->mass_t1*param->mass_t1/(param->M2_Q*param->M2_Q-param->mass_t1*param->mass_t1)/(param->mu_Q*param->mu_Q-param->mass_t1*param->mass_t1)*(log(param->mass_t1*param->mass_t1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1-1.0)) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->mass_t2*param->mass_t2/(param->M2_Q*param->M2_Q-param->mass_t2*param->mass_t2)/(param->mu_Q*param->mu_Q-param->mass_t2*param->mass_t2)*(log(param->mass_t2*param->mass_t2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2-1.0)) +param->sbot_mix[1][1]*param->sbot_mix[1][1]*param->mass_b1*param->mass_b1/(param->M2_Q*param->M2_Q-param->mass_b1*param->mass_b1)/(param->mu_Q*param->mu_Q-param->mass_b1*param->mass_b1)*(log(param->mass_b1*param->mass_b1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1-1.0))/2. +param->sbot_mix[1][2]*param->sbot_mix[1][2]*param->mass_b2*param->mass_b2/(param->M2_Q*param->M2_Q-param->mass_b2*param->mass_b2)/(param->mu_Q*param->mu_Q-param->mass_b2*param->mass_b2)*(log(param->mass_b2*param->mass_b2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2-1.0))/2.);
	}
	else
	{
		epsilonb+=1./param->inv_alpha_em/sw/sw/4.0/pi*(param->mu_Q*param->M2_Q)*((param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2) +param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1/2. +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2/2.);	
	} 

/*----------------------------------------------------------------------*/
	/* EPSILON_B' */
/*----------------------------------------------------------------------*/
	float epsilonbp=-2.0/3.0*alpha_s_MSOFT/pi*(param->mu_Q-param->A_b/param->tan_beta)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino));
		
	for(ie=1;ie<=4;ie++)
		epsilonbp+=param->yut*param->yut/16.0/pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->A_t-param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]));

	if(fabs(param->mu_Q*param->mu_Q/param->M2_Q/param->M2_Q-1.)<0.01)
	{
		epsilonbp+=1./param->inv_alpha_em/sw/sw/4.0/pi*(-param->mu_Q*param->M2_Q)* (param->stop_mix[1][1]*param->stop_mix[1][1]*param->mass_t1*param->mass_t1/(param->M2_Q*param->M2_Q-param->mass_t1*param->mass_t1)/(param->mu_Q*param->mu_Q-param->mass_t1*param->mass_t1)* (log(param->mass_t1*param->mass_t1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1-1.0))/2. +param->stop_mix[1][2]*param->stop_mix[1][2]*param->mass_t2*param->mass_t2/(param->M2_Q*param->M2_Q-param->mass_t2*param->mass_t2)/(param->mu_Q*param->mu_Q-param->mass_t2*param->mass_t2)* (log(param->mass_t2*param->mass_t2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2-1.0))/2. +param->sbot_mix[1][1]*param->sbot_mix[1][1]*param->mass_b1*param->mass_b1/(param->M2_Q*param->M2_Q-param->mass_b1*param->mass_b1)/(param->mu_Q*param->mu_Q-param->mass_b1*param->mass_b1)* (log(param->mass_b1*param->mass_b1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1-1.0)) +param->sbot_mix[1][2]*param->sbot_mix[1][2]*param->mass_b2*param->mass_b2/(param->M2_Q*param->M2_Q-param->mass_b2*param->mass_b2)/(param->mu_Q*param->mu_Q-param->mass_b2*param->mass_b2)* (log(param->mass_b2*param->mass_b2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2-1.0)));	
	}
	else
	{
		epsilonbp+=1./param->inv_alpha_em/sw/sw/4.0/pi*(param->mu_Q*param->M2_Q)*( (param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2)/2. +(param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1 +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2));	
	}

/*----------------------------------------------------------------------*/
	/* EPSILON_T' */
/*----------------------------------------------------------------------*/

	float epsilontp=-2.0/3.0*alpha_s_MSOFT/pi*(param->mu_Q+param->A_t/param->tan_beta)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino));
						
	for(ie=1;ie<=4;ie++)
		epsilontp+=param->yub*param->yub/16.0/pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]));

/*----------------------------------------------------------------------*/
/* STANDARD MODEL */ 
/*----------------------------------------------------------------------*/

	float xt= pow(mass_top/param->mass_W,2.);
	float yt=pow(mass_top/param->mass_H,2.);

	float C2SM_0 = 1.;
	float C7SM_0 = F7_1(xt);
	float C8SM_0 = F8_1(xt);

	float C7SMeps_0 = (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(xt);
	float C8SMeps_0 = (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(xt);
	
	float C7SM_1 = G7(xt); 
	float C8SM_1 = G8(xt); 

/*----------------------------------------------------------------------*/
/* CHARGED HIGGS */
/*----------------------------------------------------------------------*/
		
	float C7H_0=1./3./pow(param->tan_beta,2.)*F7_1(yt) + F7_2(yt);
	float C8H_0=1./3./pow(param->tan_beta,2.)*F8_1(yt) + F8_2(yt);

	float C7Heps_0=-(epsilontp+epsilonb)*param->tan_beta/(1.+epsilonb*param->tan_beta)*F7_2(yt);
	float C8Heps_0=-(epsilontp+epsilonb)*param->tan_beta/(1.+epsilonb*param->tan_beta)*F8_2(yt);

	float C4H_1=1./pow(param->tan_beta,2.)*F4H(yt);
	
	float C7H_1= G7H(yt,param->tan_beta)+Delta7H(yt,param->tan_beta)*log(yt/xt)-4./9.*EH(yt,param->tan_beta);
	float C8H_1= G8H(yt,param->tan_beta)+Delta8H(yt,param->tan_beta)*log(yt/xt)-1./6.*EH(yt,param->tan_beta);

/*----------------------------------------------------------------------*/
/* CHARGINOS */	
/*----------------------------------------------------------------------*/	
	
	float mass_top_SUSYscale=running_mass(mass_top,mass_top,param->MSOFT_Q,mass_top,param->mass_b,param->alpha_s_MZ,param->mass_Z);
		    
 	float r11= param->stop_mix[1][1]*param->charg_Vmix[1][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][1]*param->charg_Vmix[1][2];
	float r12= param->stop_mix[1][1]*param->charg_Vmix[2][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][1]*param->charg_Vmix[2][2];
	float r21= param->stop_mix[1][2]*param->charg_Vmix[1][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][2]*param->charg_Vmix[1][2];
	float r22= param->stop_mix[1][2]*param->charg_Vmix[2][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][2]*param->charg_Vmix[2][2];
	
	float rt11= param->charg_Vmix[1][1];
	float rt12= param->charg_Vmix[2][1];

	float rp11= param->stop_mix[1][1]*param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rp12= param->stop_mix[1][1]*param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rp21= param->stop_mix[1][2]*param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rp22= param->stop_mix[1][2]*param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	
	float rpt11= param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rpt12= param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
		
	float C7charg_0_SUSYscale = -( 2./3.*pow(r11,2.)*pow(param->mass_W/param->mass_t1,2.)*F7_1(pow(param->mass_t1/param->mass_cha1,2.))
			        +2./3.*pow(r12,2.)*pow(param->mass_W/param->mass_t1,2.)*F7_1(pow(param->mass_t1/param->mass_cha2,2.))
			        +2./3.*pow(r21,2.)*pow(param->mass_W/param->mass_t2,2.)*F7_1(pow(param->mass_t2/param->mass_cha1,2.))
			        +2./3.*pow(r22,2.)*pow(param->mass_W/param->mass_t2,2.)*F7_1(pow(param->mass_t2/param->mass_cha2,2.))
				+r11*rp11*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_t1/param->mass_cha1,2.))
				+r12*rp12*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_t1/param->mass_cha2,2.))
				+r21*rp21*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_t2/param->mass_cha1,2.))
				+r22*rp22*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_t2/param->mass_cha2,2.)))
				+( 2./3.*pow(rt11,2.)*pow(param->mass_W/param->mass_upl,2.)*F7_1(pow(param->mass_upl/param->mass_cha1,2.))
			        +2./3.*pow(rt12,2.)*pow(param->mass_W/param->mass_upl,2.)*F7_1(pow(param->mass_upl/param->mass_cha2,2.))
				+rt11*rpt11*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_upl/param->mass_cha1,2.))
				+rt12*rpt12*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_upl/param->mass_cha2,2.)));	
	
	float C8charg_0_SUSYscale = -( 2./3.*pow(r11,2.)*pow(param->mass_W/param->mass_t1,2.)*F8_1(pow(param->mass_t1/param->mass_cha1,2.))
			        +2./3.*pow(r12,2.)*pow(param->mass_W/param->mass_t1,2.)*F8_1(pow(param->mass_t1/param->mass_cha2,2.))
			        +2./3.*pow(r21,2.)*pow(param->mass_W/param->mass_t2,2.)*F8_1(pow(param->mass_t2/param->mass_cha1,2.))
				+2./3.*pow(r22,2.)*pow(param->mass_W/param->mass_t2,2.)*F8_1(pow(param->mass_t2/param->mass_cha2,2.))
				+r11*rp11*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_t1/param->mass_cha1,2.))
				+r12*rp12*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_t1/param->mass_cha2,2.))
				+r21*rp21*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_t2/param->mass_cha1,2.))
				+r22*rp22*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_t2/param->mass_cha2,2.)))
				+( 2./3.*pow(rt11,2.)*pow(param->mass_W/param->mass_upl,2.)*F8_1(pow(param->mass_upl/param->mass_cha1,2.))
			        +2./3.*pow(rt12,2.)*pow(param->mass_W/param->mass_upl,2.)*F8_1(pow(param->mass_upl/param->mass_cha2,2.))
				+rt11*rpt11*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_upl/param->mass_cha1,2.))
				+rt12*rpt12*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_upl/param->mass_cha2,2.)));	
	
 	float eta_MSOFTQ =alpha_s_MSOFT/alpha_s_MW; 
	
	float beta0=-7.;
	
	float C7charg_0 = pow(eta_MSOFTQ,-16./3./beta0)*C7charg_0_SUSYscale +8./3.*(pow(eta_MSOFTQ,-14./3./beta0)-pow(eta_MSOFTQ,-16./3./beta0))*C8charg_0_SUSYscale;

	float C8charg_0 = pow(eta_MSOFTQ,-14./3./beta0)*C8charg_0_SUSYscale;

	float C7charg_1=0.;
	float C8charg_1=0.;

	float t21=param->stop_mix[1][2]*param->charg_Vmix[1][1]+mass_top/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[1][1]*param->charg_Vmix[1][2];
	float t22=param->stop_mix[1][2]*param->charg_Vmix[2][1]+mass_top/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[1][1]*param->charg_Vmix[2][2];
	
	float C4charg_1 = pow(param->mass_W/param->mass_t2,2.)*(t21*t21*F4ch(param->mass_t2/param->mass_cha1)+t22*t22*F4ch(param->mass_t2/param->mass_cha2));
	
/*----------------------------------------------------------------------*/
/* C TOTAL */	
/*----------------------------------------------------------------------*/

#ifdef SMONLY
	C7SMeps_0=C7H_0=C7Heps_0=C7charg_0=C7H_1=C7charg_1=0.;
	C8SMeps_0=C8H_0=C8Heps_0=C8charg_0=C8H_1=C8charg_1=0.;
	C4H_1=C4charg_1=0.;
#endif

        float C2_0 = C2SM_0;
        float C7_0 = C7SM_0 + C7SMeps_0 +C7H_0+C7Heps_0+C7charg_0;
        float C8_0 = C8SM_0 + C8SMeps_0 +C8H_0+C8Heps_0+C8charg_0;
 		
        float C7_1 = C7SM_1 +C7H_1+C7charg_1;
        float C8_1 = C8SM_1 +C8H_1+C8charg_1;

/*----------------------------------------------------------------------*/
/* MW -> mu */
/*----------------------------------------------------------------------*/		
	float alpha_s_mu=alpha_s_running(mu,param->mass_top_pole,mass_b_pole,param->alpha_s_MZ,param->mass_Z);	
	
	float eta_mu=alpha_s_MW/alpha_s_mu;
		

 	float C1_0mu= 1./2.*pow(eta_mu,6./23.) -1./2.*pow(eta_mu,-12./23.);
	float C2_0mu= 1./2.*pow(eta_mu,6./23.) +1./2.*pow(eta_mu,-12./23.);
	float C3_0mu= -1./14.*pow(eta_mu,6./23.) +1./6.*pow(eta_mu,-12./23.) +0.0510*pow(eta_mu,0.4086) -0.1403*pow(eta_mu,-0.4230) -0.0113*pow(eta_mu,-0.8994) +0.0054*pow(eta_mu,0.1456);
	float C4_0mu= -1./14.*pow(eta_mu,6./23.) -1./6.*pow(eta_mu,-12./23.) +0.0984*pow(eta_mu,0.4086) +0.1214*pow(eta_mu,-0.4230) +0.0156*pow(eta_mu,-0.8994) +0.0026*pow(eta_mu,0.1456);
	float C5_0mu= -0.0397*pow(eta_mu,0.4086) +0.0117*pow(eta_mu,-0.4230) -0.0025*pow(eta_mu,-0.8994) +0.0304*pow(eta_mu,0.1456);
	float C6_0mu= 0.0335*pow(eta_mu,0.4086) +0.0239*pow(eta_mu,-0.4230) -0.0462*pow(eta_mu,-0.8994) -0.0112*pow(eta_mu,0.1456);
 
	float C7_0mu= pow(eta_mu,16./23.)*C7_0 + 8./3.*(pow(eta_mu,14./23.)-pow(eta_mu,16./23.))*C8_0 + C2_0 * (2.2996*pow(eta_mu,14./23.) -1.0880*pow(eta_mu,16./23.) -3./7.*pow(eta_mu,6./23.) -1./14.*pow(eta_mu,-12./23.) -0.6494*pow(eta_mu,0.4086) -0.0380*pow(eta_mu,-0.4230) -0.0185*pow(eta_mu,-0.8994) -0.0057*pow(eta_mu,0.1456));

	float C8_0mu= pow(eta_mu,14./23.)*C8_0 + C2_0 * (0.8623*pow(eta_mu,14./23.) -0.9135*pow(eta_mu,0.4086) +0.0873*pow(eta_mu,-0.4230) -0.0571*pow(eta_mu,-0.8994) +0.0209*pow(eta_mu,0.1456));
		
	/* NLO */
	float C1_1mu= (0.8136+1.0197*eta_mu)*pow(eta_mu,6./23.)+(0.7142+2.9524*eta_mu)*pow(eta_mu,-12./23.);
	
	float C2_1mu= (0.8136+1.0197*eta_mu)*pow(eta_mu,6./23.)-(0.7142+2.9524*eta_mu)*pow(eta_mu,-12./23.);
	
	float C3_1mu= (-0.0766-0.1457*eta_mu)*pow(eta_mu,6./23.)+(-0.1455-0.9841*eta_mu)*pow(eta_mu,-12./23.)
	        +(0.1494*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.8848+0.2303*eta_mu)*pow(eta_mu,0.4086)
		+(-0.3726*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.4137+1.4672*eta_mu)*pow(eta_mu,-0.4230)
		+(0.0738*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.0114+0.0971*eta_mu)*pow(eta_mu,-0.8994)
		+(-0.0173*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.1722-0.0213*eta_mu)*pow(eta_mu,0.1456);
	
	float C4_1mu= (-0.2353-0.1457*eta_mu)*pow(eta_mu,6./23.)+(-0.0397+0.9841*eta_mu)*pow(eta_mu,-12./23.)
	        +(0.2885*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.4920+0.4447*eta_mu)*pow(eta_mu,0.4086)
		+(0.3224*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.2758-1.2696*eta_mu)*pow(eta_mu,-0.4230)
		+(-0.1025*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.0019-0.1349*eta_mu)*pow(eta_mu,-0.8994)
		+(-0.0084*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.1449-0.0104*eta_mu)*pow(eta_mu,0.1456);
	
	float C5_1mu= 0.0397*pow(eta_mu,6./23.)+0.0926*pow(eta_mu,-12./23.)
	        +(-0.1163*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.7342-0.1792*eta_mu)*pow(eta_mu,0.4086)
		+(0.0310*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.1262-0.1221*eta_mu)*pow(eta_mu,-0.4230)
		+(0.0162*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.1209+0.0213*eta_mu)*pow(eta_mu,-0.8994)
		+(-0.0975*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.1085-0.1197*eta_mu)*pow(eta_mu,0.1456);
	
	float C6_1mu= -0.1191*pow(eta_mu,6./23.)-0.2778*pow(eta_mu,-12./23.)
	        +(0.0982*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.5544+0.1513*eta_mu)*pow(eta_mu,0.4086)
		+(0.0634*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.1915-0.2497*eta_mu)*pow(eta_mu,-0.4230)
		+(0.3026*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.2744+0.3983*eta_mu)*pow(eta_mu,-0.8994)
		+(0.0358*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.3568+0.0440*eta_mu)*pow(eta_mu,0.1456);

	float C7_1mu= pow(eta_mu,39./23.)*C7_1+8./3.*(pow(eta_mu,37./23.)-pow(eta_mu,39./23.))*C8_1    
		+(297664./14283.*pow(eta_mu,16./23.)-7164416./357075.*pow(eta_mu,14./23.)+256868./14283.*pow(eta_mu,37./23.)-6698884./357075.*pow(eta_mu,39./23.))*C8_0
		+37208./4761.*(pow(eta_mu,39./23.)-pow(eta_mu,16./23.))*C7_0
		+(4661194./816831.*eta_mu*(E(xt)+(C4H_1+C4charg_1))-17.3023+14.8088*eta_mu)*pow(eta_mu,14./23.)
		+(-8516./2217.*eta_mu*(E(xt)+(C4H_1+C4charg_1))+8.5027-10.8090*eta_mu)*pow(eta_mu,16./23.)
		+ (4.5508-0.8740*eta_mu)*pow(eta_mu,6./23.)
		+ (0.7519+0.4218*eta_mu)*pow(eta_mu,-12./23.)
		+ (-1.9043*eta_mu*(E(xt)+(C4H_1+C4charg_1))+2.0040-2.9347*eta_mu)*pow(eta_mu,0.4086)
		+ (-0.1008*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.7476+0.3971*eta_mu)*pow(eta_mu,-0.4230)
		+ (0.1216*eta_mu*(E(xt)+(C4H_1+C4charg_1))-0.5385+0.1600*eta_mu)*pow(eta_mu,-0.8994)
		+ (0.0183*eta_mu*(E(xt)+(C4H_1+C4charg_1))+0.0914+0.0225*eta_mu)*pow(eta_mu,0.1456);
	
	float C8_1mu=0.; 

	/* inverting the conventions of C1 and C2 */
	
	C0[1] = C2_0mu;
	C0[2] = C1_0mu;
	C0[3] = C3_0mu;
	C0[4] = C4_0mu;
	C0[5] = C5_0mu;
	C0[6] = C6_0mu;
	C0[7] = C7_0mu;
	C0[8] = C8_0mu;
	C1[1] = C2_1mu;
	C1[2] = C1_1mu;
	C1[3] = C3_1mu;
	C1[4] = C4_1mu;
	C1[5] = C5_1mu;
	C1[6] = C6_1mu;
	C1[7] = C7_1mu;
	C1[8] = C8_1mu;

	return;	
}

/*--------------------------------------------------------------------*/

void Cem_calculator(float Cem[], float mu, struct parameters* param)
/* calculates of the "electromagnetic Wilson coefficients" needed for the calculation of the branching ratio, at scale mu, using the parameters of the structure param */
{	
	int ie;
	
	float pi=4.*atan(1.);

	for(ie=1;ie<=8;ie++) Cem[ie]=0.;
		
	float alpha_s_mb=alpha_s_running(param->mass_b,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);	

	float mass_b_pole=param->mass_b*(1+alpha_s_mb/pi*(4./3. +alpha_s_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*(1.61/4.62))+alpha_s_mb/pi*(190.595-4.*(26.655-4.*0.6527)))));

	float alpha_s_MW=alpha_s_running(param->mass_W,param->mass_top_pole,mass_b_pole,param->alpha_s_MZ,param->mass_Z); 
	
	float alpha_s_mtop=alpha_s_running(param->mass_top_pole,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z); 
	
	float mtmt=param->mass_top_pole*(1.-4./3.*alpha_s_mtop/pi);
	
	float mass_top=running_mass(mtmt,mtmt,param->mass_W,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);

	float sw=sin(atan(param->gp/param->g2));

/*----------------------------------------------------------------------*/
	/* EPSILON_B */
/*----------------------------------------------------------------------*/

	float alpha_s_MSOFT=alpha_s_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);

	float epsilonb=-2.0/3.0*alpha_s_MSOFT/pi*(param->mu_Q-param->A_b/param->tan_beta)/param->mass_gluino*H2(param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) -param->yut*param->yut/16.0/pi/pi*(param->A_t-param->mu_Q/param->tan_beta)*(param->charg_Umix[1][2]*param->charg_Vmix[1][2]/param->mass_cha1 *H2(param->mass_t1*param->mass_t1/param->mass_cha1/param->mass_cha1,param->mass_t2*param->mass_t2/param->mass_cha1/param->mass_cha1) +param->charg_Umix[2][2]*param->charg_Vmix[2][2]/param->mass_cha2*H2(param->mass_t1*param->mass_t1/param->mass_cha2/param->mass_cha2,param->mass_t2*param->mass_t2/param->mass_cha2/param->mass_cha2));

	
	if(fabs(param->mu_Q*param->mu_Q/param->M2_Q/param->M2_Q-1)<0.01)
	{
		epsilonb+=1./param->inv_alpha_em/sw/sw/4.0/pi*(-param->mu_Q*param->M2_Q)*(param->stop_mix[1][1]*param->stop_mix[1][1]*param->mass_t1*param->mass_t1/(param->M2_Q*param->M2_Q-param->mass_t1*param->mass_t1)/(param->mu_Q*param->mu_Q-param->mass_t1*param->mass_t1)*(log(param->mass_t1*param->mass_t1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1-1.0)) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->mass_t2*param->mass_t2/(param->M2_Q*param->M2_Q-param->mass_t2*param->mass_t2)/(param->mu_Q*param->mu_Q-param->mass_t2*param->mass_t2)*(log(param->mass_t2*param->mass_t2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2-1.0)) +param->sbot_mix[1][1]*param->sbot_mix[1][1]*param->mass_b1*param->mass_b1/(param->M2_Q*param->M2_Q-param->mass_b1*param->mass_b1)/(param->mu_Q*param->mu_Q-param->mass_b1*param->mass_b1)*(log(param->mass_b1*param->mass_b1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1-1.0))/2. +param->sbot_mix[1][2]*param->sbot_mix[1][2]*param->mass_b2*param->mass_b2/(param->M2_Q*param->M2_Q-param->mass_b2*param->mass_b2)/(param->mu_Q*param->mu_Q-param->mass_b2*param->mass_b2)*(log(param->mass_b2*param->mass_b2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2-1.0))/2.);
	}
	else
	{
		epsilonb+=1./param->inv_alpha_em/sw/sw/4.0/pi*(param->mu_Q*param->M2_Q)*((param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2) +param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1/2. +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2/2.);	
	} 


/*----------------------------------------------------------------------*/
	/* EPSILON_B' */
/*----------------------------------------------------------------------*/

	float epsilonbp=-2.0/3.0*alpha_s_MSOFT/pi*(param->mu_Q-param->A_b/param->tan_beta)/param->mass_gluino*(
	    param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino));
		
	for(ie=1;ie<=4;ie++)
		epsilonbp+=param->yut*param->yut/16.0/pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->A_t-param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]));

	if(fabs(param->mu_Q*param->mu_Q/param->M2_Q/param->M2_Q-1.)<0.01)
	{
		epsilonbp+=1./param->inv_alpha_em/sw/sw/4.0/pi*(-param->mu_Q*param->M2_Q)* (param->stop_mix[1][1]*param->stop_mix[1][1]*param->mass_t1*param->mass_t1/(param->M2_Q*param->M2_Q-param->mass_t1*param->mass_t1)/(param->mu_Q*param->mu_Q-param->mass_t1*param->mass_t1)* (log(param->mass_t1*param->mass_t1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1-1.0))/2. +param->stop_mix[1][2]*param->stop_mix[1][2]*param->mass_t2*param->mass_t2/(param->M2_Q*param->M2_Q-param->mass_t2*param->mass_t2)/(param->mu_Q*param->mu_Q-param->mass_t2*param->mass_t2)* (log(param->mass_t2*param->mass_t2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2-1.0))/2. +param->sbot_mix[1][1]*param->sbot_mix[1][1]*param->mass_b1*param->mass_b1/(param->M2_Q*param->M2_Q-param->mass_b1*param->mass_b1)/(param->mu_Q*param->mu_Q-param->mass_b1*param->mass_b1)* (log(param->mass_b1*param->mass_b1/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1-1.0)) +param->sbot_mix[1][2]*param->sbot_mix[1][2]*param->mass_b2*param->mass_b2/(param->M2_Q*param->M2_Q-param->mass_b2*param->mass_b2)/(param->mu_Q*param->mu_Q-param->mass_b2*param->mass_b2)* (log(param->mass_b2*param->mass_b2/param->M2_Q/param->M2_Q)+(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2-1.0)));	
	}
	else
	{
		epsilonbp+=1./param->inv_alpha_em/sw/sw/4.0/pi*(param->mu_Q*param->M2_Q)*( (param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2)/2. +(param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1 +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2));	
	}

/*----------------------------------------------------------------------*/
	/* EPSILON_T' */
/*----------------------------------------------------------------------*/

	float epsilontp=-2.0/3.0*alpha_s_MSOFT/pi*(param->mu_Q+param->A_t/param->tan_beta)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino));
						
	for(ie=1;ie<=4;ie++)
		epsilontp+=param->yub*param->yub/16.0/pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]));
	 
/*----------------------------------------------------------------------*/
/* STANDARD MODEL */ 
/*----------------------------------------------------------------------*/

	float xt= pow(mass_top/param->mass_W,2.);
	float yt=pow(mass_top/param->mass_H,2.);

	float C7SM_0 = F7_1(xt);
	float C8SM_0 = F8_1(xt);

	float C7SMeps_0 = (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(xt);
	float C8SMeps_0 = (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(xt);

/*----------------------------------------------------------------------*/
/* CHARGED HIGGS */
/*----------------------------------------------------------------------*/
		
	float C7H_0=1./3./pow(param->tan_beta,2.)*F7_1(yt) + F7_2(yt);
	float C8H_0=1./3./pow(param->tan_beta,2.)*F8_1(yt) + F8_2(yt);

	float C7Heps_0=-(epsilontp+epsilonb)*param->tan_beta/(1.+epsilonb*param->tan_beta)*F7_2(yt);
	float C8Heps_0=-(epsilontp+epsilonb)*param->tan_beta/(1.+epsilonb*param->tan_beta)*F8_2(yt);

/*----------------------------------------------------------------------*/
/* CHARGINOS */	
/*----------------------------------------------------------------------*/	
	
	float mass_top_SUSYscale=running_mass(mass_top,mass_top,param->MSOFT_Q,mass_top,param->mass_b,param->alpha_s_MZ,param->mass_Z);			
	    
 	float r11= param->stop_mix[1][1]*param->charg_Vmix[1][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][1]*param->charg_Vmix[1][2];
	float r12= param->stop_mix[1][1]*param->charg_Vmix[2][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][1]*param->charg_Vmix[2][2];
	float r21= param->stop_mix[1][2]*param->charg_Vmix[1][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][2]*param->charg_Vmix[1][2];
	float r22= param->stop_mix[1][2]*param->charg_Vmix[2][1] + mass_top_SUSYscale/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][2]*param->charg_Vmix[2][2];
	
	float rt11= param->charg_Vmix[1][1];
	float rt12= param->charg_Vmix[2][1];

	float rp11= param->stop_mix[1][1]*param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rp12= param->stop_mix[1][1]*param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rp21= param->stop_mix[1][2]*param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rp22= param->stop_mix[1][2]*param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	
	float rpt11= param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
	float rpt12= param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta))/(1.+epsilonb*param->tan_beta);
		
	float C7charg_0_SUSYscale = -( 2./3.*pow(r11,2.)*pow(param->mass_W/param->mass_t1,2.)*F7_1(pow(param->mass_t1/param->mass_cha1,2.))
			        +2./3.*pow(r12,2.)*pow(param->mass_W/param->mass_t1,2.)*F7_1(pow(param->mass_t1/param->mass_cha2,2.))
			        +2./3.*pow(r21,2.)*pow(param->mass_W/param->mass_t2,2.)*F7_1(pow(param->mass_t2/param->mass_cha1,2.))
			        +2./3.*pow(r22,2.)*pow(param->mass_W/param->mass_t2,2.)*F7_1(pow(param->mass_t2/param->mass_cha2,2.))
				+r11*rp11*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_t1/param->mass_cha1,2.))
				+r12*rp12*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_t1/param->mass_cha2,2.))
				+r21*rp21*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_t2/param->mass_cha1,2.))
				+r22*rp22*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_t2/param->mass_cha2,2.)))
				+( 2./3.*pow(rt11,2.)*pow(param->mass_W/param->mass_upl,2.)*F7_1(pow(param->mass_upl/param->mass_cha1,2.))
			        +2./3.*pow(rt12,2.)*pow(param->mass_W/param->mass_upl,2.)*F7_1(pow(param->mass_upl/param->mass_cha2,2.))
				+rt11*rpt11*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_upl/param->mass_cha1,2.))
				+rt12*rpt12*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_upl/param->mass_cha2,2.)));	
	
	float C8charg_0_SUSYscale = -( 2./3.*pow(r11,2.)*pow(param->mass_W/param->mass_t1,2.)*F8_1(pow(param->mass_t1/param->mass_cha1,2.))
			        +2./3.*pow(r12,2.)*pow(param->mass_W/param->mass_t1,2.)*F8_1(pow(param->mass_t1/param->mass_cha2,2.))
			        +2./3.*pow(r21,2.)*pow(param->mass_W/param->mass_t2,2.)*F8_1(pow(param->mass_t2/param->mass_cha1,2.))
			        +2./3.*pow(r22,2.)*pow(param->mass_W/param->mass_t2,2.)*F8_1(pow(param->mass_t2/param->mass_cha2,2.))
				+r11*rp11*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_t1/param->mass_cha1,2.))
				+r12*rp12*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_t1/param->mass_cha2,2.))
				+r21*rp21*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_t2/param->mass_cha1,2.))
				+r22*rp22*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_t2/param->mass_cha2,2.)))
				+( 2./3.*pow(rt11,2.)*pow(param->mass_W/param->mass_upl,2.)*F8_1(pow(param->mass_upl/param->mass_cha1,2.))
			        +2./3.*pow(rt12,2.)*pow(param->mass_W/param->mass_upl,2.)*F8_1(pow(param->mass_upl/param->mass_cha2,2.))
				+rt11*rpt11*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_upl/param->mass_cha1,2.))
				+rt12*rpt12*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_upl/param->mass_cha2,2.)));	
	

 	float eta_MSOFTQ =alpha_s_MSOFT/alpha_s_MW; 
	
	float beta0=-7.;
	
	float C7charg_0 = pow(eta_MSOFTQ,-16./3./beta0)*C7charg_0_SUSYscale +8./3.*(pow(eta_MSOFTQ,-14./3./beta0)-pow(eta_MSOFTQ,-16./3./beta0))*C8charg_0_SUSYscale;

	float C8charg_0 = pow(eta_MSOFTQ,-14./3./beta0)*C8charg_0_SUSYscale;
	

/*----------------------------------------------------------------------*/
/* C TOTAL */	
/*----------------------------------------------------------------------*/

#ifdef SMONLY
	C7SMeps_0=C7H_0=C7Heps_0=C7charg_0=0.;
	C8SMeps_0=C8H_0=C8Heps_0=C8charg_0=0.;
#endif


        float C7_0 = C7SM_0 + C7SMeps_0 +C7H_0+C7Heps_0+C7charg_0;
        float C8_0 = C8SM_0 + C8SMeps_0 +C8H_0+C8Heps_0+C8charg_0;
 		
/*----------------------------------------------------------------------*/
/* MW -> mu */
/*----------------------------------------------------------------------*/		
	float alpha_s_mu=alpha_s_running(mu,mass_top,param->mass_b,param->alpha_s_MZ,param->mass_Z);	
	
	float eta_mu=alpha_s_MW/alpha_s_mu;
		
	Cem[2]= -190./8073.*pow(eta_mu,-35./23.)-359./3105. *pow(eta_mu,-17./23.)+4276./121095. *pow(eta_mu,-12./23.) +350531./1009125. *pow(eta_mu,-9./23.)+2./4347. *pow(eta_mu,-7./23.)-5956./15525.*pow(eta_mu,6./23.) +38380./169533. *pow(eta_mu,14./23.)-748/8625. *pow(eta_mu,16./23.);
	
	Cem[8]= -32./575.*pow(eta_mu,-9./23.)+32./1449.*pow(eta_mu,-7./23.)+640./1449.*pow(eta_mu,14./23.)-704./1725.*pow(eta_mu,16./23.);
			
	Cem[7]= (32./75.*pow(eta_mu,-9/23.)-40/69.*pow(eta_mu,-7/23.)+88./575.*pow(eta_mu,16/23.))*C7_0+Cem[8]*C8_0+Cem[2];

	return;	
}
