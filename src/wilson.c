#include "include.h"

/*--------------------------------------------------------------------*/

double A0t(double x)
{
	return (-46.+159.*x-153.*x*x+22.*x*x*x)/36./pow(1.-x,3.)+x*x*(-3.*x+2.)/2./pow(1.-x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double A1t(double x, double l)
{
	return (32.*pow(x,4.)+244.*pow(x,3.)-160.*x*x+16.*x)/9./pow(1.-x,4.)*Li2(1.-1./x)
	+(-774.*x*x*x*x-2826.*x*x*x+1994.*x*x-130.*x+8.)/81./pow(1.-x,5.)*log(x)
	+(-94.*x*x*x*x-18665.*x*x*x+20682.*x*x-9113.*x+2006.)/243./pow(1.-x,4.)
	+((-12.*x*x*x*x-92.*x*x*x+56.*x*x)/3./pow(1.-x,5.)*log(x)+(-68.*x*x*x*x-202.*x*x*x-804.*x*x+794.*x-152.)/27./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double F0t(double x)
{
	return  (5.*x*x*x-9.*x*x+30.*x-8.)/12./pow(1.-x,3.)+3.*x*x/2./pow(1.-x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F1t(double x,double l)
{
	return (4.*pow(x,4.)-40.*pow(x,3.)-41.*x*x-x)/3./pow(1.-x,4.)*Li2(1.-1./x)
	+(-144.*x*x*x*x+3177.*x*x*x+3661.*x*x+250.*x-32.)/108./pow(1.-x,5.)*log(x)
	+(-247.*x*x*x*x+11890.*x*x*x+31779.*x*x-2966.*x+1016.)/648./pow(1.-x,4.)
	+((17.*x*x*x+31.*x*x)/pow(1.-x,5.)*log(x)+(-35.*x*x*x*x+170.*x*x*x+447.*x*x+338.*x-56.)/18./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double E0t(double x)
{
	return (-9.*x*x+16.*x-4.)/6./pow(1.-x,4.)*log(x)+(-7.*x*x*x-21.*x*x+42.*x+4.)/36./pow(1.-x,3.);
}

/*----------------------------------------------------------------------*/

double G1t(double x, double l)
{
	return (10.*pow(x,4.)-100.*pow(x,3.)+30.*x*x+160.*x-40.)/27./pow(1.-x,4.)*Li2(1.-1./x)
	+(30.*pow(x,3.)-42.*x*x-332.*x+68.)/81./pow(1.-x,4.)*log(x)
	+(-6.*pow(x,3.)-293.*x*x+161.*x+42.)/81./pow(1.-x,3.)
	+((90.*x*x-160.*x+40.)/27./pow(1.-x,4.)*log(x)+(35.*pow(x,3.)+105.*x*x-210.*x-20.)/81./pow(1.-x,3.))*l;
}

/*----------------------------------------------------------------------*/

double E1t(double x, double l)
{
	return (515.*pow(x,4.)-614.*pow(x,3.)-81.*x*x-190.*x+40.)/54./pow(1.-x,4.)*Li2(1.-1./x)
	+(-1030.*pow(x,4.)+435.*pow(x,3.)+1373.*x*x+1950.*x-424.)/108./pow(1.-x,5.)*log(x)
	+(-29467.*pow(x,4.)+45604.*pow(x,3.)-30237.*x*x+66532.*x-10960.)/1944./pow(1.-x,4.)
	+((-1125.*pow(x,3.)+1685.*x*x+380.*x-76.)/54./pow(1.-x,5.)*log(x)
	+(133.*pow(x,4.)-2758.*pow(x,3.)-2061.*x*x+11522.*x-1652.)/324./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double T(double x)
{
	return -(16.*x+8.)*sqrt(4.*x-1.)*Cl2(2.*asin(0.5/sqrt(x)))+(16.*x+20./3.)*log(x)+32.*x+112./9.;
}

/*----------------------------------------------------------------------*/

double Ech(double x)
{
	return x*(11.-7.*x+2.*x*x)/18./pow(x-1.,3.)-x/3./pow(x-1.,4.)*log(fabs(x));
}

/*----------------------------------------------------------------------*/

double F7_1(double x)
{
	return x*(7.-5.*x-8.*x*x)/24./pow(x-1.,3.)+x*x*(3.*x-2.)/4./pow(x-1.,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F7_2(double x)
{
	return x*(3.-5.*x)/12./pow(x-1.,2.)+x*(3.*x-2.)/6./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double F7_3(double x)
{
	return (5.-7.*x)/6./pow(x-1.,2.)+x*(3.*x-2.)/3./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_1(double x)
{
	return x*(2.+5.*x-x*x)/8./pow(x-1.,3.)-3.*x*x/4./pow(x-1.,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_2(double x)
{
	return x*(3.-x)/4./pow(x-1.,2.)-x/2./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_3(double x)
{
	return (1.+x)/2./pow(x-1.,2.)-x/pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double H2(double x, double y)
{
	return D2(x,y);
}

/*----------------------------------------------------------------------*/

double B(double m1, double m2, double Q)
{
	double x=pow(m2/m1,2.);

	if(fabs(x-1.)<1.e-5) return -0.5*log(m2*m2/Q/Q);

	return 0.5*(0.5+1./(1.-x)+log(x)/pow(1.-x,2.)-log(m2*m2/Q/Q));
}

/*----------------------------------------------------------------------*/

double G7H(double x, double lu, double ld)
{
	return lu*(ld*4.*x*(4.*(-3.+7.*x-2.*x*x)*Li2(1.-1./x)+(8.-14.*x-3.*x*x)*log(x)*log(x)/(x-1.)
	+2.*(-3.-x+12.*x*x-2*x*x*x)*log(x)/(x-1.)+3.*(7.-13.*x+2*x*x))/9./pow((x-1.),3.)
	+lu*2.*x*(x*(18.-37.*x+8.*x*x)*Li2(1.-1./x)+x*(-14.+23.*x+3.*x*x)*log(x)*log(x)/(x-1.)+
	(-50.+251.*x-174.*x*x-192.*x*x*x+21.*x*x*x*x)*log(x)/9./(x-1.)+(797.-5436.*x+7569.*x*x-1202.*x*x*x)/108.)/pow((x-1.),4.)/9.);
}

/*----------------------------------------------------------------------*/
	
double Delta7H(double x, double lu, double ld)
{
	return 2.*x/9./pow((x-1.),4.)*lu*
	(ld*((x-1.)*(21.-47.*x+8.*x*x)+2.*(-8.+14.*x+3.*x*x)*log(x))
	+lu*((-31.-18.*x+135.*x*x-14.*x*x*x)/6.+x*(14.-23.*x-3.*x*x)*log(x)/(x-1.)));
}

/*----------------------------------------------------------------------*/

double EH(double x, double lu)
{
	return lu*lu*x*((x-1.)*(16.-29.*x+7.*x*x)+6.*(3.*x-2.)*log(x))/36./pow((x-1.),4.);
}

/*----------------------------------------------------------------------*/

double G8H(double x, double lu, double ld)
{
	return lu*(ld*x*(0.5*(-36.+25.*x-17*x*x)*Li2(1.-1./x)+(19.+17.*x)*log(x)*log(x)/(x-1.)
	+0.25*(-3.-187.*x+12.*x*x-14.*x*x*x)*log(x)/(x-1.)+3.*(143.-44.*x+29.*x*x)/8.)/3./pow((x-1.),3.)
	+ lu*x*(x*(30.-17.*x+13.*x*x)*Li2(1.-1./x)-x*(31.+17.*x)*log(x)*log(x)/(x-1.)+
	(-226.+817.*x+1353.*x*x+318.*x*x*x+42.*x*x*x*x)*log(x)/36./(x-1.)+(1130.-18153.*x+7650.*x*x-4451.*x*x*x)/216.)/pow((x-1.),4.)/6.);
}

/*----------------------------------------------------------------------*/

double Delta8H(double x, double lu, double ld)
{
	return (x/6./pow((x-1.),4.))*lu*(ld*((x-1.)*(81.-16.*x+7.*x*x)-2.*(19.+17.*x)*log(x))
	+lu*((-38.-261.*x+18.*x*x-7.*x*x*x)/6.+x*(31.+17.*x)*log(x)/(x-1.)));	
}

/*----------------------------------------------------------------------*/

double C7t2mt(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.4) return 12.06+12.93*z+3.013*z*log(z)+96.71*z*z+52.73*z*z*log(z)+147.9*pow(z,3.)
	+187.7*pow(z,3.)*log(z)-144.9*pow(z,4.)+236.1*pow(z,4.)*log(z);
	
	else return 11.74+0.3642*w+0.1155*w*w-0.003145*pow(w,3.)-0.03263*pow(w,4.)-0.03528*pow(w,5.)
	-0.03076*pow(w,6.)-0.02504*pow(w,7.)-0.01985*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C7c2MW(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.4) return 1.525-0.1165*z+0.01975*z*log(z)+0.06283*z*z+0.005349*z*z*log(z)
	+0.01005*pow(z*log(z),2.)-0.04202*pow(z,3.)+0.01535*pow(z,3.)*log(z)-0.00329*z*pow(z*log(z),2.)
	+0.002372*pow(z,4.)-0.0007910*pow(z,4.)*log(z);
		
	else return 1.432+0.06709*w+0.01257*w*w+0.004710*pow(w,3.)+0.002373*pow(w,4.)
	+0.001406*pow(w,5.)+0.0009216*pow(w,6.)+0.00064730*pow(w,7.)+0.0004779*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C8t2mt(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.35) return -0.8954-7.043*z-98.34*z*z-46.21*z*z*log(z)-127.1*pow(z,3.)
	-181.6*pow(z,3.)*log(z)+535.8*pow(z,4.)-76.76*pow(z,4.)*log(z);
	
	else return -0.6141-0.8975*w-0.03492*w*w+0.06791*pow(w,3.)+0.07966*pow(w,4.)
	+0.07226*pow(w,5.)+0.06132*pow(w,6.)+0.05096*pow(w,7.)+0.04216*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C8c2MW(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.35) return -1.870+0.1010*z-0.1218*z*log(z)+0.1045*z*z-0.03748*z*z*log(z)
	+0.01151*pow(z*log(z),2.)-0.01023*pow(z,3.)+0.004342*pow(z,3.)*log(z)+0.0003031*z*pow(z*log(z),2.)
	-0.001537*pow(z,4.)+0.0007532*pow(z,4.)*log(z);
	
	else return -1.676-0.1179*w-0.02926*w*w-0.01297*pow(w,3.)-0.007296*pow(w,4.)
	-0.004672*pow(w,5.)-0.003248*pow(w,6.)-0.002389*pow(w,7.)-0.001831*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double epsilon_0(struct parameters* param)
{
	if(param->SM==1) return 0;

#ifdef SM_ChargedHiggs
	return 0;
#endif

 	double sw=sin(atan(param->gp/param->g2));
	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);

	return 2./3.*alphas_MSOFT/pi*((param->A_b/param->tan_beta-param->mu_Q)/param->mass_gluino*
H2(param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino)	
	-0.5*(B(param->mass_gluino,param->mass_b1, param->MSOFT_Q)+B(param->mass_gluino,param->mass_b2, param->MSOFT_Q))/param->tan_beta )
	 +1./param->inv_alpha_em/sw/sw/4./pi*(param->mu_Q*param->M2_Q)*( +param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1/2. +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2/2.);
}

/*----------------------------------------------------------------------*/

double epsilon_2(struct parameters* param)
{
	if(param->SM==1) return 0;

#ifdef SM_ChargedHiggs
	return 0;
#endif

 	double sw=sin(atan(param->gp/param->g2));

	return param->yut[3]*param->yut[3]/16./pi/pi*(param->mu_Q/param->tan_beta-param->A_t)*(param->charg_Umix[1][2]*param->charg_Vmix[1][2]/param->mass_cha1 *H2(param->mass_t1*param->mass_t1/param->mass_cha1/param->mass_cha1,param->mass_t2*param->mass_t2/param->mass_cha1/param->mass_cha1) +param->charg_Umix[2][2]*param->charg_Vmix[2][2]/param->mass_cha2*H2(param->mass_t1*param->mass_t1/param->mass_cha2/param->mass_cha2,param->mass_t2*param->mass_t2/param->mass_cha2/param->mass_cha2))	+1./param->inv_alpha_em/sw/sw/4./pi*(param->mu_Q*param->M2_Q)*(param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2);
}

/*----------------------------------------------------------------------*/

double epsilon_b(struct parameters* param)
{
	return epsilon_0(param)+epsilon_2(param);
}

/*----------------------------------------------------------------------*/

double epsilon_bp(struct parameters* param)
{
	if(param->SM==1) return 0;

#ifdef SM_ChargedHiggs
	return 0;
#endif
 	double sw=sin(atan(param->gp/param->g2));
	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	int ie;
	int nb_neut;
	if(param->mass_neut[5]==0.) nb_neut=4; else nb_neut=5;

	double epsilonbp=2./3.*alphas_MSOFT/pi*(param->A_b/param->tan_beta-param->mu_Q)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino));	
	for(ie=1;ie<=nb_neut;ie++)	epsilonbp+=param->yut[3]*param->yut[3]/16./pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->A_t-param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]));
	epsilonbp+=1./param->inv_alpha_em/sw/sw/4./pi*(param->mu_Q*param->M2_Q)*( (param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2)/2. +(param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1 +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2));

	return epsilonbp;
}

/*----------------------------------------------------------------------*/

double epsilon_0p(struct parameters* param)
{
	if(param->SM==1) return 0;

#ifdef SM_ChargedHiggs
	return 0;
#endif

	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	int ie;
	int nb_neut;
	if(param->mass_neut[5]==0.) nb_neut=4; else nb_neut=5;
	
	double epsilon0p=-2./3.*alphas_MSOFT/pi*(param->mu_Q+param->A_t/param->tan_beta)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino));			
	for(ie=1;ie<=nb_neut;ie++)	epsilon0p+=param->yub[3]*param->yub[3]/16./pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]));

	return epsilon0p;
}

/*----------------------------------------------------------------------*/

double epsilon_1p(struct parameters* param)
{
	if(param->SM==1) return 0;

#ifdef SM_ChargedHiggs
	return 0;
#endif

	return 1./16./pi/pi*(param->yub[3]*param->yub[3]*param->A_b/param->mu_Q*H2(pow(param->MqL3_Q/param->mu_Q,2.),pow(param->MbR_Q/param->mu_Q,2.))
	-param->g2*param->g2*param->M2_Q/param->mu_Q*H2(pow(param->MqL3_Q/param->mu_Q,2.),pow(param->M2_Q/param->mu_Q,2.)));
}

/*----------------------------------------------------------------------*/

void CW_calculator(double C0w[], double C1w[], double C2w[], double mu_W, struct parameters* param)
/* calculates the LO (C0w), NLO (C1w) and NNLO (C2w) contributions to the Wilson coefficients at scale mu_W, using the parameters of the structure param */
{
	int ie;
	for(ie=1;ie<=8;ie++) C0w[ie]=C1w[ie]=C2w[ie]=0.;
	
	double mass_top_muW=running_mass(param->mtmt,param->mtmt,mu_W,param->mass_top_pole,param->mass_b,param);

	double epsilonbp,epsilon0p,epsilon0,epsilon2,epsilon1p,epsilonb;
	if(param->THDM_model==0)
	{	
		epsilonbp=epsilon_bp(param);
		epsilon0p=epsilon_0p(param);
		epsilon0=epsilon_0(param);
		epsilon2=epsilon_2(param);
		epsilon1p=epsilon_1p(param);
		epsilonb=epsilon0+epsilon2;	
	}

/*----------------------------------------------------------------------*/
/* STANDARD MODEL */ 
/*----------------------------------------------------------------------*/

	double L=log(mu_W*mu_W/param->mass_W/param->mass_W);

	double xt= pow(mass_top_muW/param->mass_W,2.);
	double yt= pow(mass_top_muW/param->mass_H,2.);

	double C2SM_0 = 1.;
	double C7SM_0 = -0.5*A0t(xt)-23./36.;
	double C8SM_0 = -0.5*F0t(xt)-1./3.;
	
	double C7SMeps_0,C8SMeps_0;
	if(param->THDM_model==0)
	{
		C7SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(xt);
		C8SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(xt);
	}

	double C1SM_1 = 15.+6.*L;
	double C4SM_1 = E0t(xt)-7./9.+2./3.*L;
	double C7SM_1 = -0.5*A1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+713./243.+4./81.*L-4./9.*C4SM_1;
	double C8SM_1 = -0.5*F1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+91./324.-4./27.*L-C4SM_1/6.;

	double C1SM_2 = -T(xt)+7987./72.+17.*pi*pi/3.+475./6.*L+17.*L*L;
	double C2SM_2 = 127./18.+4./3.*pi*pi+46./3.*L+4.*L*L;
	double C3SM_2 = G1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))-680./243.-20./81.*pi*pi-68./81.*L-20./27.*L*L;
	double C4SM_2 = E1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+950./243.+10./81.*pi*pi+124./27.*L+10./27.*L*L;
	double C5SM_2 = -G1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))/10.+2./15.*E0t(xt)+68./243.+2./81.*pi*pi+14./81.*L+2./27.*L*L;
	double C6SM_2 = -3./16.*G1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+E0t(xt)/4.+85./162.+5./108.*pi*pi+35./108.*L+5./36.*L*L;

	double xtW=pow(running_mass(param->mtmt,param->mtmt,param->mass_W,param->mass_top_pole,param->mass_b,param)/param->mass_W,2.);
	double xtt=pow(param->mtmt/param->mass_W,2.);
	
	double C7SM_2 = (C7t2mt(xtt)+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-592.*pow(xt,5.)-22.*pow(xt,4.)+12814.*pow(xt,3.)-6376.*xt*xt+512.*xt)/27./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-26838.*pow(xt,5.)+25938.*pow(xt,4.)+627367.*pow(xt,3.)-331956.*xt*xt+16989.*xt-460.)/729./pow(xt-1.,6.)*log(xt)
	+(34400.*pow(xt,5.)+276644.*pow(xt,4.)-2668324.*pow(xt,3.)+1694437.*xt*xt-323354.*xt+53077.)/2187./pow(xt-1.,5.)
	+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-63.*pow(xt,5.)+532.*pow(xt,4.)+2089.*pow(xt,3.)-1118.*xt*xt)/9./pow(xt-1.,6.)*log(xt)
	+(1186.*pow(xt,5.)-2705.*pow(xt,4.)-24791.*pow(xt,3.)-16099.*xt*xt+19229.*xt-2740.)/162./pow(xt-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(mu_W*mu_W/param->mass_W/param->mass_W)+814./729.*pow(log(mu_W*mu_W/param->mass_W/param->mass_W),2.));

	double C8SM_2 = (C8t2mt(xtt)+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-148.*pow(xt,5.)+1052.*pow(xt,4.)-4811.*pow(xt,3.)-3520.*xt*xt-61.*xt)/18./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-15984.*pow(xt,5.)+152379.*pow(xt,4.)-1358060.*pow(xt,3.)-1201653.*xt*xt-74190.*xt+9188.)/1944./pow(xt-1.,6.)*log(xt)
	+(109669.*pow(xt,5.)-1112675.*pow(xt,4.)+6239377.*pow(xt,3.)+8967623.*xt*xt+768722.*xt-42796.)/11664./pow(xt-1.,5.)
	+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-139.*pow(xt,4.)-2938.*pow(xt,3.)-2683.*xt*xt)/12./pow(xt-1.,6.)*log(xt)
	+(1295.*pow(xt,5.)-7009.*pow(xt,4.)+29495.*pow(xt,3.)+64513.*xt*xt+17458.*xt-2072.)/216./pow(xt-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(mu_W*mu_W/param->mass_W/param->mass_W)+397./486.*pow(log(mu_W*mu_W/param->mass_W/param->mass_W),2.));
	
/*----------------------------------------------------------------------*/
/* CHARGED HIGGS */
/*----------------------------------------------------------------------*/
	
	double C7Heps_0,C8Heps_0,C7Heps2_0,C8Heps2_0,mass_b_muW;
	double lu,ld;
	if(param->THDM_model==0)
	{
		C7Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(yt);
		C8Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(yt);

		C7Heps2_0;
		C8Heps2_0;
	
		mass_b_muW=running_mass(param->mass_b, param->mass_b,mu_W,param->mass_top_pole,param->mass_b,param);

		C7Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.		*((cos(param->alpha)+sin(param->alpha)*param->tan_beta)*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*param->tan_beta)*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan(param->tan_beta))-sin(atan(param->tan_beta))*param->tan_beta)*(sin(atan(param->tan_beta))-epsilonb*cos(atan(param->tan_beta)))*pow(mass_b_muW/param->mass_A0,2.));
	
		C8Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F8_2(yt);
		C8Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.		*((cos(param->alpha)+sin(param->alpha)*param->tan_beta)*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*param->tan_beta)*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan(param->tan_beta))-sin(atan(param->tan_beta))*param->tan_beta)*(sin(atan(param->tan_beta))-epsilonb*cos(atan(param->tan_beta)))*pow(mass_b_muW/param->mass_A0,2.));

		lu=1./param->tan_beta;
		ld=-param->tan_beta;
	}
	else
	{
		lu=param->lambda_u[3][3];
		ld=param->lambda_d[3][3];
	}

	double C7H_0=1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
	double C8H_0=1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

	double C4H_1=EH(yt,lu);
	
 	double C7H_1= G7H(yt,lu,ld)+Delta7H(yt,lu,ld)*log(pow(mu_W/param->mass_H,2.))-4./9.*C4H_1;
	double C8H_1= G8H(yt,lu,ld)+Delta8H(yt,lu,ld)*log(pow(mu_W/param->mass_H,2.))-1./6.*C4H_1;

/*----------------------------------------------------------------------*/
/* CHARGINOS */	
/*----------------------------------------------------------------------*/	
	double r11,r12,r21,r22,C4charg_1,rt11,rt12,rp11,rp12,rp21,rp22,rpt11,rpt12;
	double C7charg_0,C8charg_0,C7_chargeps_0,C8_chargeps_0,C7charg_1,C8charg_1;
	if(param->THDM_model==0)
	{
		r11= param->stop_mix[1][1]*param->charg_Vmix[1][1] + mass_top_muW/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][1]*param->charg_Vmix[1][2];
		r12= param->stop_mix[1][1]*param->charg_Vmix[2][1] + mass_top_muW/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][1]*param->charg_Vmix[2][2];
		r21= param->stop_mix[1][2]*param->charg_Vmix[1][1] + mass_top_muW/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][2]*param->charg_Vmix[1][2];
		r22= param->stop_mix[1][2]*param->charg_Vmix[2][1] + mass_top_muW/sqrt(2.)/param->mass_W/sin(atan(param->tan_beta))*param->stop_mix[2][2]*param->charg_Vmix[2][2];

		C4charg_1 =  pow(param->mass_W/param->mass_t2,2.)*(r21*r21*Ech(param->mass_t2/param->mass_cha1)+r22*r22*Ech(param->mass_t2/param->mass_cha2));
	
		rt11= param->charg_Vmix[1][1];
		rt12= param->charg_Vmix[2][1];

		rp11= param->stop_mix[1][1]*param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta));
		rp12= param->stop_mix[1][1]*param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta));
		rp21= param->stop_mix[1][2]*param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta));
		rp22= param->stop_mix[1][2]*param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta));
	
		rpt11= param->charg_Umix[1][2]/sqrt(2.)/cos(atan(param->tan_beta));
		rpt12= param->charg_Umix[2][2]/sqrt(2.)/cos(atan(param->tan_beta));	
		
		C7charg_0 = -( 2./3.*pow(r11,2.)*pow(param->mass_W/param->mass_t1,2.)*F7_1(pow(param->mass_t1/param->mass_cha1,2.))
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
	
		C8charg_0 = -( 2./3.*pow(r11,2.)*pow(param->mass_W/param->mass_t1,2.)*F8_1(pow(param->mass_t1/param->mass_cha1,2.))
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
 	
		C7_chargeps_0= 	-epsilonb/(1.+epsilonb*param->tan_beta)*param->tan_beta*(-r11*rp11*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_t1/param->mass_cha1,2.))
		-r12*rp12*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_t1/param->mass_cha2,2.))
		-r21*rp21*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_t2/param->mass_cha1,2.))
		-r22*rp22*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_t2/param->mass_cha2,2.))
		+(rt11*rpt11*param->mass_W/param->mass_cha1*F7_3(pow(param->mass_upl/param->mass_cha1,2.))
		+rt12*rpt12*param->mass_W/param->mass_cha2*F7_3(pow(param->mass_upl/param->mass_cha2,2.))));
	
		C8_chargeps_0= 	-epsilonb/(1.+epsilonb*param->tan_beta)*param->tan_beta*(-r11*rp11*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_t1/param->mass_cha1,2.))
		-r12*rp12*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_t1/param->mass_cha2,2.))
		-r21*rp21*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_t2/param->mass_cha1,2.))
		-r22*rp22*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_t2/param->mass_cha2,2.))
		+(rt11*rpt11*param->mass_W/param->mass_cha1*F8_3(pow(param->mass_upl/param->mass_cha1,2.))
		+rt12*rpt12*param->mass_W/param->mass_cha2*F8_3(pow(param->mass_upl/param->mass_cha2,2.))));
	
		C7charg_1=0.;
		C8charg_1=0.;
	}

/*----------------------------------------------------------------------*/
/* C TOTAL */	
/*----------------------------------------------------------------------*/

	if(param->SM==1) 
	{
		C7SMeps_0=C7H_0=C7Heps_0=C7Heps2_0=C7charg_0=C7_chargeps_0=C7H_1=C7charg_1=0.;
		C8SMeps_0=C8H_0=C8Heps_0=C8Heps2_0=C8charg_0=C8_chargeps_0=C8H_1=C8charg_1=0.;
		C4H_1=C4charg_1=0.;
	}

#ifdef SM_ChargedHiggs
	C7SMeps_0=C7Heps_0=C7Heps2_0=C7charg_0=C7_chargeps_0=C7charg_1=0.;
	C8SMeps_0=C8Heps_0=C8Heps2_0=C8charg_0=C8_chargeps_0=C8charg_1=0.;
	C4charg_1=0.;
#endif
	
	if(param->THDM_model==0)
	{
		C0w[2]= C2SM_0;
		C0w[7]= C7SM_0+C7SMeps_0+C7H_0+C7Heps_0+C7Heps2_0+C7charg_0+C7_chargeps_0;
		C0w[8]= C8SM_0+C8SMeps_0+C8H_0+C8Heps_0+C8Heps2_0+C8charg_0+C8_chargeps_0;
		C1w[1]= C1SM_1;
		C1w[4]= C4SM_1+C4H_1+C4charg_1;	
		C1w[7]= C7SM_1+C7H_1+C7charg_1;
		C1w[8]= C8SM_1+C8H_1+C8charg_1;
		C2w[1]= C1SM_2;
		C2w[2]= C2SM_2;
		C2w[3]= C3SM_2;
		C2w[4]= C4SM_2;
		C2w[5]= C5SM_2;
		C2w[6]= C6SM_2;
		C2w[7]= C7SM_2;
		C2w[8]= C8SM_2;
		return;
	}

	else
	{
		C0w[2]= C2SM_0;
		C0w[7]= C7SM_0+C7H_0;
		C0w[8]= C8SM_0+C8H_0;
		C1w[1]= C1SM_1;
		C1w[4]= C4SM_1+C4H_1;	
		C1w[7]= C7SM_1+C7H_1;
		C1w[8]= C8SM_1+C8H_1;
		C2w[1]= C1SM_2;
		C2w[2]= C2SM_2;
		C2w[3]= C3SM_2;
		C2w[4]= C4SM_2;
		C2w[5]= C5SM_2;
		C2w[6]= C6SM_2;
		C2w[7]= C7SM_2;
		C2w[8]= C8SM_2;
		return;
	}


}

/*-----------------------------------------------------------------------------*/

void C_calculator_base1(double C0w[], double C1w[], double C2w[], double mu_W, double C0b[], double C1b[], double C2b[], double mu, struct parameters* param)
/* calculates the LO (C0b), NLO (C1b) and NNLO (C2b) contributions to the Wilson coefficients at scale mu, using the LO (C0w), NLO (C1w) and NNLO (C2w) contributions to the Wilson coefficients at scale mu_W and the parameters of the structure param, in the standard operator basis */
{
	int ie;

	for(ie=1;ie<=8;ie++) C0b[ie]=C1b[ie]=C2b[ie]=0.;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	
	double eta_mu=alphas_muW/alphas_mu;
		
 	C0b[1]= C0w[2]* (pow(eta_mu,6./23.) - pow(eta_mu,-12./23.));
	
 	C0b[2]= C0w[2]* (0.6667*pow(eta_mu,6./23.) + 0.3333*pow(eta_mu,-12./23.));
	
 	C0b[3]= C0w[2]* (0.0317*pow(eta_mu,6./23.) - 0.0370*pow(eta_mu,-12./23.)
	-0.0659*pow(eta_mu,0.4086)+0.0595*pow(eta_mu,-0.4230)-0.0218*pow(eta_mu,-0.8994)+0.0335*pow(eta_mu,0.1456));
	
	C0b[4]= C0w[2]* (0.0476*pow(eta_mu,6./23.) + 0.1111*pow(eta_mu,-12./23.)
	+0.0237*pow(eta_mu,0.4086)-0.0173*pow(eta_mu,-0.4230)-0.1336*pow(eta_mu,-0.8994)-0.0316*pow(eta_mu,0.1456));
	
	C0b[5]= C0w[2]* (-0.0079*pow(eta_mu,6./23.) + 0.0093*pow(eta_mu,-12./23.)
	+0.0094*pow(eta_mu,0.4086)-0.0100*pow(eta_mu,-0.4230)+0.0010*pow(eta_mu,-0.8994)-0.0017*pow(eta_mu,0.1456));

	C0b[6]= C0w[2]* (-0.0119*pow(eta_mu,6./23.) - 0.0278*pow(eta_mu,-12./23.)
	+0.0108*pow(eta_mu,0.4086)+0.0163*pow(eta_mu,-0.4230)+0.0103*pow(eta_mu,-0.8994)+0.0023*pow(eta_mu,0.1456));
 
	C0b[7]= C0w[2]*(2.2996*pow(eta_mu,14./23.) - 1.0880*pow(eta_mu,16./23.)
	-0.4286*pow(eta_mu,6./23.) - 0.0714*pow(eta_mu,-12./23.)
	-0.6494*pow(eta_mu,0.4086)-0.0380*pow(eta_mu,-0.4230)
	-0.0185*pow(eta_mu,-0.8994)-0.0057*pow(eta_mu,0.1456))
	+C0w[7]*pow(eta_mu,16./23.)
	+C0w[8]*(2.6667*pow(eta_mu,14./23.) - 2.6667*pow(eta_mu,16./23.));
	
	C0b[8]= C0w[2]*(0.8623*pow(eta_mu,14./23.) 
	-0.9135*pow(eta_mu,0.4086)+0.0873*pow(eta_mu,-0.4230)
	-0.0571*pow(eta_mu,-0.8994)+0.0209*pow(eta_mu,0.1456))
	+C0w[8]*pow(eta_mu,14./23.);
			
	
	C1b[1]=eta_mu*(
	C1w[1]*(0.3333*pow(eta_mu,6./23.) +0.6667*pow(eta_mu,-12./23.))
	+C0w[2]*(-2.9606*pow(eta_mu,6./23.) - 4.0951*pow(eta_mu,-12./23.)
	+5.9606*pow(eta_mu,6./23.-1.) +1.0951*pow(eta_mu,-12./23.-1.)));	
	
	C1b[2]=eta_mu*(
	C1w[1]*(0.2222*pow(eta_mu,6./23.) -0.2222*pow(eta_mu,-12./23.))
	+C0w[2]*(-1.9737*pow(eta_mu,6./23.) +1.3650*pow(eta_mu,-12./23.)
	+1.9737*pow(eta_mu,6./23.-1.) -1.3650*pow(eta_mu,-12./23.-1.)));
	
	C1b[3]= eta_mu*(
	C1w[1]*(0.0106*pow(eta_mu,6./23.) +0.0247*pow(eta_mu,-12./23.)
	-0.0129*pow(eta_mu,0.4086)-0.0497*pow(eta_mu,-0.4230)
	+0.0092*pow(eta_mu,-0.8994)+0.0182*pow(eta_mu,0.1456))
	+C1w[4]*(-0.1933*pow(eta_mu,0.4086)+0.1579*pow(eta_mu,-0.4230)
	+0.1428*pow(eta_mu,-0.8994)-0.1074*pow(eta_mu,0.1456))
	+C0w[2]*(-0.0940*pow(eta_mu,6./23.) -0.1517*pow(eta_mu,-12./23.)
	-0.2327*pow(eta_mu,0.4086)+0.2288*pow(eta_mu,-0.4230)
	+0.1455*pow(eta_mu,-0.8994)-0.4760*pow(eta_mu,0.1456)
	-0.5409*pow(eta_mu,6./23.-1.)+1.6332*pow(eta_mu,-12./23.-1.)
	+1.6406*pow(eta_mu,0.4086-1.)-1.6702*pow(eta_mu,-0.4230-1.)
	-0.2576*pow(eta_mu,-0.8994-1.)-0.2250*pow(eta_mu,0.1456-1.)));
	
	C1b[4]= eta_mu*(
	C1w[1]*(0.0159*pow(eta_mu,6./23.) -0.0741*pow(eta_mu,-12./23.)
	+0.0046*pow(eta_mu,0.4086)+0.0144*pow(eta_mu,-0.4230)
	+0.0562*pow(eta_mu,-0.8994)-0.0171*pow(eta_mu,0.1456))
	+C1w[4]*(0.0695*pow(eta_mu,0.4086)-0.0459*pow(eta_mu,-0.4230)
	+0.8752*pow(eta_mu,-0.8994)+0.1012*pow(eta_mu,0.1456))
	+C0w[2]*(-0.1410*pow(eta_mu,6./23.) +0.4550*pow(eta_mu,-12./23.)
	+0.0836*pow(eta_mu,0.4086)-0.0664*pow(eta_mu,-0.4230)
	+0.8919*pow(eta_mu,-0.8994)+0.4485*pow(eta_mu,0.1456)
	+2.2203*pow(eta_mu,6./23.-1.)+2.0265*pow(eta_mu,-12./23.-1.)
	-4.1830*pow(eta_mu,0.4086-1.)-0.7135*pow(eta_mu,-0.4230-1.)
	-1.8215*pow(eta_mu,-0.8994-1.)+0.7996*pow(eta_mu,0.1456-1.)));
	
	C1b[5]= eta_mu*(
	C1w[1]*(-0.0026*pow(eta_mu,6./23.) -0.0062*pow(eta_mu,-12./23.)
	+0.0018*pow(eta_mu,0.4086)+0.0083*pow(eta_mu,-0.4230)
	-0.0004*pow(eta_mu,-0.8994)-0.0009*pow(eta_mu,0.1456))
	+C1w[4]*(0.0274*pow(eta_mu,0.4086)-0.0264*pow(eta_mu,-0.4230)
	-0.0064*pow(eta_mu,-0.8994)+0.0055*pow(eta_mu,0.1456))
	+C0w[2]*(0.0235*pow(eta_mu,6./23.) +0.0379*pow(eta_mu,-12./23.)
	+0.0330*pow(eta_mu,0.4086)-0.0383*pow(eta_mu,-0.4230)
	-0.0066*pow(eta_mu,-0.8994)+0.0242*pow(eta_mu,0.1456)
	+0.0400*pow(eta_mu,6./23.-1.)-0.1861*pow(eta_mu,-12./23.-1.)
	-0.1669*pow(eta_mu,0.4086-1.)+0.1887*pow(eta_mu,-0.4230-1.)
	+0.0201*pow(eta_mu,-0.8994-1.)+0.0304*pow(eta_mu,0.1456-1.)));
	
	C1b[6]= eta_mu*(
	C1w[1]*(-0.0040*pow(eta_mu,6./23.) +0.0185*pow(eta_mu,-12./23.)
	+0.0021*pow(eta_mu,0.4086)-0.0136*pow(eta_mu,-0.4230)
	-0.0043*pow(eta_mu,-0.8994)+0.0012*pow(eta_mu,0.1456))
	+C1w[4]*(0.0317*pow(eta_mu,0.4086)+0.0432*pow(eta_mu,-0.4230)
	-0.0675*pow(eta_mu,-0.8994)-0.0074*pow(eta_mu,0.1456))
	+C0w[2]*(0.0352*pow(eta_mu,6./23.) -0.1138*pow(eta_mu,-12./23.)
	+0.0382*pow(eta_mu,0.4086)+0.0625*pow(eta_mu,-0.4230)
	-0.0688*pow(eta_mu,-0.8994)-0.0327*pow(eta_mu,0.1456)
	-0.2614*pow(eta_mu,6./23.-1.)-0.1918*pow(eta_mu,-12./23.-1.)
	+0.4197*pow(eta_mu,0.4086-1.)+0.0295*pow(eta_mu,-0.4230-1.)
	+0.1474*pow(eta_mu,-0.8994-1.)-0.0640*pow(eta_mu,0.1456-1.)));
	
	C1b[7]=eta_mu*(
	C1w[1]*(0.5784*pow(eta_mu,14./23.) - 0.3921*pow(eta_mu,16./23.)
	-0.1429*pow(eta_mu,6./23.) +0.0476*pow(eta_mu,-12./23.)
	-0.1275*pow(eta_mu,0.4086)+0.0317*pow(eta_mu,-0.4230)
	+0.0078*pow(eta_mu,-0.8994)-0.0031*pow(eta_mu,0.1456))
	+C1w[4]*(5.7064*pow(eta_mu,14./23.) - 3.8412*pow(eta_mu,16./23.)
	-1.9043*pow(eta_mu,0.4086)-0.1008*pow(eta_mu,-0.4230)
	+0.1216*pow(eta_mu,-0.8994)+0.0183*pow(eta_mu,0.1456))
	+C1w[7]*pow(eta_mu,16./23.)
	+C1w[8]*(2.6667*pow(eta_mu,14./23.) - 2.6667*pow(eta_mu,16./23.))
	+C0w[2]*(9.9372*pow(eta_mu,14./23.) - 7.4878*pow(eta_mu,16./23.)
	+1.2688*pow(eta_mu,6./23.) -0.2925*pow(eta_mu,-12./23.)
	-2.2923*pow(eta_mu,0.4086)-0.1461*pow(eta_mu,-0.4230)
	+0.1239*pow(eta_mu,-0.8994)+0.0812*pow(eta_mu,0.1456)
	-17.3023*pow(eta_mu,14./23.-1.) +8.5027*pow(eta_mu,16./23.-1.)
	+4.5508*pow(eta_mu,6./23.-1.)+0.7519*pow(eta_mu,-12./23.-1.)
	+2.0040*pow(eta_mu,0.4086-1.)+0.7476*pow(eta_mu,-0.4230-1.)
	-0.5385*pow(eta_mu,-0.8994-1.)+0.0914*pow(eta_mu,0.1456-1.))
	+C0w[7]*(7.8152*pow(eta_mu,16./23.)-7.8152*pow(eta_mu,16./23.-1))
	+C0w[8]*(17.9842*pow(eta_mu,14./23.) - 18.7604*pow(eta_mu,16./23.)
	-20.0642*pow(eta_mu,14./23.-1.) +20.8404*pow(eta_mu,16./23.-1.)));
	
	C1b[8]=eta_mu*(
	C1w[1]*(0.2169*pow(eta_mu,14./23.)
	-0.1793*pow(eta_mu,0.4086)-0.0730*pow(eta_mu,-0.4230)
	+0.0240*pow(eta_mu,-0.8994)+0.0113*pow(eta_mu,0.1456))
	+C1w[4]*(2.1399*pow(eta_mu,14./23.)
	-2.6788*pow(eta_mu,0.4086)+0.2318*pow(eta_mu,-0.4230)
	+0.3741*pow(eta_mu,-0.8994)-0.0670*pow(eta_mu,0.1456))
	+C1w[8]*pow(eta_mu,14./23.)
	+C0w[2]*(3.7264*pow(eta_mu,14./23.)
	-3.2247*pow(eta_mu,0.4086)+0.3359*pow(eta_mu,-0.4230)
	+0.3812*pow(eta_mu,-0.8994)-0.2968*pow(eta_mu,0.1456)
	-5.8157*pow(eta_mu,14./23.-1.)
	+1.4062*pow(eta_mu,6./23.-1.)-3.9895*pow(eta_mu,-12./23.-1.)
	+3.2850*pow(eta_mu,0.4086-1.)+3.6851*pow(eta_mu,-0.4230-1.)
	-0.1424*pow(eta_mu,-0.8994-1.)+0.6492*pow(eta_mu,0.1456-1.))
	+C0w[8]*(6.7441*pow(eta_mu,14./23.) - 6.7441*pow(eta_mu,14./23.-1.)));
	
			
	C2b[7]=eta_mu*eta_mu*(
		
	C2w[1]*(0.5784*pow(eta_mu,14./23.) - 0.3921*pow(eta_mu,16./23.)
	-0.1429*pow(eta_mu,6./23.) +0.0476*pow(eta_mu,-12./23.)
	-0.1275*pow(eta_mu,0.4086)+0.0317*pow(eta_mu,-0.4230)
	+0.0078*pow(eta_mu,-0.8994)-0.0031*pow(eta_mu,0.1456))
	+C2w[2]*(2.2996*pow(eta_mu,14./23.) - 1.0880*pow(eta_mu,16./23.)
	-0.4286*pow(eta_mu,6./23.) - 0.0714*pow(eta_mu,-12./23.)
	-0.6494*pow(eta_mu,0.4086)-0.0380*pow(eta_mu,-0.4230)
	-0.0185*pow(eta_mu,-0.8994)-0.0057*pow(eta_mu,0.1456))
	+C2w[3]*(8.0780*pow(eta_mu,14./23.) - 5.2777*pow(eta_mu,16./23.)
	-2.8536*pow(eta_mu,0.4086)+0.1281*pow(eta_mu,-0.4230)
	+0.1495*pow(eta_mu,-0.8994)-0.2244*pow(eta_mu,0.1456))
	+C2w[4]*(5.7064*pow(eta_mu,14./23.) - 3.8412*pow(eta_mu,16./23.)
	-1.9043*pow(eta_mu,0.4086)-0.1008*pow(eta_mu,-0.4230)
	+0.1216*pow(eta_mu,-0.8994)+0.0183*pow(eta_mu,0.1456))
	+C2w[5]*(202.9010*pow(eta_mu,14./23.) - 149.4668*pow(eta_mu,16./23.)
	-55.2813*pow(eta_mu,0.4086)+2.6494*pow(eta_mu,-0.4230)
	+0.7191*pow(eta_mu,-0.8994)-1.5213*pow(eta_mu,0.1456))
	+C2w[6]*(86.4618*pow(eta_mu,14./23.) - 59.6604*pow(eta_mu,16./23.)
	-25.4430*pow(eta_mu,0.4086)-1.2894*pow(eta_mu,-0.4230)
	+0.0228*pow(eta_mu,-0.8994)-0.0917*pow(eta_mu,0.1456))
	+C2w[7]*pow(eta_mu,16./23.)
	+C2w[8]*(2.6667*pow(eta_mu,14./23.) - 2.6667*pow(eta_mu,16./23.))
	
	+C1w[1]*(0.0021*pow(eta_mu,14./23.)-1.4498*pow(eta_mu,16./23.)
	+0.8515*pow(eta_mu,6./23.) +0.0521*pow(eta_mu,-12./23.)
	+0.6707*pow(eta_mu,0.4086)+0.1220*pow(eta_mu,-0.4230)
	-0.0578*pow(eta_mu,-0.8994)+0.0355*pow(eta_mu,0.1456)
	-4.3519*pow(eta_mu,14./23.-1.) +3.0646*pow(eta_mu,16./23.-1.)
	+1.5169*pow(eta_mu,6./23.-1.)-0.5013*pow(eta_mu,-12./23.-1.)
	+0.3934*pow(eta_mu,0.4086-1.)-0.6245*pow(eta_mu,-0.4230-1.)
	+0.2268*pow(eta_mu,-0.8994-1.)+0.0496*pow(eta_mu,0.1456-1.))
	+C1w[4]*(-8.6840*pow(eta_mu,14./23.) +8.5586*pow(eta_mu,16./23.)
	+0.7579*pow(eta_mu,0.4086)+0.4446*pow(eta_mu,-0.4230)
	+0.3093*pow(eta_mu,-0.8994)+0.4318*pow(eta_mu,0.1456)
	-42.9356*pow(eta_mu,14./23.-1.) +30.0198*pow(eta_mu,16./23.-1.)
	+5.8768*pow(eta_mu,0.4086-1.)+1.9845*pow(eta_mu,-0.4230-1.)
	+3.5291*pow(eta_mu,-0.8994-1.)-0.2929*pow(eta_mu,0.1456-1.))
	+C1w[7]*(7.8152*pow(eta_mu,16./23.)-7.8152*pow(eta_mu,16./23.-1))
	+C1w[8]*(17.9842*pow(eta_mu,14./23.) - 18.7604*pow(eta_mu,16./23.)
	-20.0642*pow(eta_mu,14./23.-1.) +20.8404*pow(eta_mu,16./23.-1.))
	
	+C0w[2]*(-212.4136*pow(eta_mu,14./23.)+167.6577*pow(eta_mu,16./23.)
	+5.7465*pow(eta_mu,6./23.) -3.7262*pow(eta_mu,-12./23.)
	+28.8574*pow(eta_mu,0.4086)-2.1262*pow(eta_mu,-0.4230)
	+2.2903*pow(eta_mu,-0.8994)+0.1462*pow(eta_mu,0.1456)
	-74.7681*pow(eta_mu,14./23.-1.)+58.5182*pow(eta_mu,16./23.-1.)
	-13.4731*pow(eta_mu,6./23.-1.) +3.0791*pow(eta_mu,-12./23.-1.)
	+7.0744*pow(eta_mu,0.4086-1.)+2.8757*pow(eta_mu,-0.4230-1.)
	+3.5962*pow(eta_mu,-0.8994-1.)-1.2982*pow(eta_mu,0.1456-1.)
	+31.4443*pow(eta_mu,14./23.-2.)-18.1165*pow(eta_mu,16./23.-2.)
	+23.2117*pow(eta_mu,6./23.-2.) +13.2771*pow(eta_mu,-12./23.-2.)
	-19.8699*pow(eta_mu,0.4086-2.)+4.0279*pow(eta_mu,-0.4230-2.)
	-8.6259*pow(eta_mu,-0.8994-2.)+2.6149*pow(eta_mu,0.1456-2.))
	+C0w[7]*(44.4252*pow(eta_mu,16./23.)-61.0768*pow(eta_mu,16./23.-1)+16.6516*pow(eta_mu,16./23.-2.))
	+C0w[8]*(15.4051*pow(eta_mu,14./23.)-18.7662*pow(eta_mu,16./23.)
	-135.3141*pow(eta_mu,14./23.-1.)+146.6159*pow(eta_mu,16./23.-1.)
	+36.4636*pow(eta_mu,14./23.-2.)-44.4043*pow(eta_mu,16./23.-2.)));
	
	return;	
}

/*-----------------------------------------------------------------------------*/

void C_calculator_base2(double C0w[], double C1w[], double mu_W, double C0b[], double C1b[], double mu, struct parameters* param)
/* calculates the LO (C0b) and NLO (C1b) contributions to the Wilson coefficients at scale mu, using the LO (C0w) and NLO (C1w) contributions to the Wilson coefficients at scale mu_W and the parameters of the structure param, in the traditional operator basis */
{
	int ie;
	for(ie=1;ie<=8;ie++) C0b[ie]=C1b[ie]=0.;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	

	double eta_mu=alphas_muW/alphas_mu;
	
	C1w[7]-=-4./9.*C1w[4];
	C1w[8]-=-C1w[4]/6.;
			

 	C0b[1]= (1./2.*pow(eta_mu,6./23.) -1./2.*pow(eta_mu,-12./23.))*C0w[2];
	C0b[2]= (1./2.*pow(eta_mu,6./23.) +1./2.*pow(eta_mu,-12./23.))*C0w[2];
	C0b[3]= (-1./14.*pow(eta_mu,6./23.) +1./6.*pow(eta_mu,-12./23.) +0.0509*pow(eta_mu,0.4086) -0.1403*pow(eta_mu,-0.4230) -0.01126*pow(eta_mu,-0.8994) +0.0054*pow(eta_mu,0.1456))*C0w[2];
	C0b[4]= (-1./14.*pow(eta_mu,6./23.) -1./6.*pow(eta_mu,-12./23.) +0.0984*pow(eta_mu,0.4086) +0.1214*pow(eta_mu,-0.4230) +0.0156*pow(eta_mu,-0.8994) +0.0026*pow(eta_mu,0.1456))*C0w[2];
	C0b[5]= (-0.0397*pow(eta_mu,0.4086) +0.0117*pow(eta_mu,-0.4230) -0.0025*pow(eta_mu,-0.8994) +0.0304*pow(eta_mu,0.1456))*C0w[2];
	C0b[6]= (0.0335*pow(eta_mu,0.4086) +0.0239*pow(eta_mu,-0.4230) -0.0462*pow(eta_mu,-0.8994) -0.0112*pow(eta_mu,0.1456))*C0w[2];
 
	C0b[7]= pow(eta_mu,16./23.)*C0w[7] + 8./3.*(pow(eta_mu,14./23.)-pow(eta_mu,16./23.))*C0w[8] + C0w[2] * (2.2996*pow(eta_mu,14./23.) -1.0880*pow(eta_mu,16./23.) -3./7.*pow(eta_mu,6./23.) -1./14.*pow(eta_mu,-12./23.) -0.6494*pow(eta_mu,0.4086) -0.0380*pow(eta_mu,-0.4230) -0.0185*pow(eta_mu,-0.8994) -0.0057*pow(eta_mu,0.1456));

	C0b[8]= pow(eta_mu,14./23.)*C0w[7] + C0w[2] * (0.8623*pow(eta_mu,14./23.) -0.9135*pow(eta_mu,0.4086) +0.0873*pow(eta_mu,-0.4230) -0.0571*pow(eta_mu,-0.8994) +0.0209*pow(eta_mu,0.1456));
		

	C1b[1]= (C0w[2] *0.8136+1.0197*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(C0w[2] *0.7142+2.9524*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.);
	
	C1b[2]= (C0w[2] *0.8136+1.0197*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)-(C0w[2] *0.7142+2.9524*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.);
	
	C1b[3]= (-0.0766*C0w[2]-0.1457*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(-0.1455*C0w[2]-0.9841*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
	        +(0.1494*eta_mu*C1w[4]-0.8848*C0w[2]+0.2303*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(-0.3726*eta_mu*C1w[4]+0.4137*C0w[2]+1.4672*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.0738*eta_mu*C1w[4]-0.0114*C0w[2]+0.0971*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0173*eta_mu*C1w[4]+0.1722*C0w[2]-0.0213*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[4]= (-0.2353*C0w[2]-0.1457*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(-0.0397*C0w[2]+0.9841*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
	        +(0.2885*eta_mu*C1w[4]+0.4920*C0w[2]+0.4447*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.3224*eta_mu*C1w[4]-0.2758*C0w[2]-1.2696*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(-0.1025*eta_mu*C1w[4]+0.0019*C0w[2]-0.1349*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0084*eta_mu*C1w[4]-0.1449*C0w[2]-0.0104*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[5]= 0.0397*C0w[2]*pow(eta_mu,6./23.)+0.0926*C0w[2]*pow(eta_mu,-12./23.)
	        +(-0.1163*eta_mu*C1w[4]+0.7342*C0w[2]-0.1792*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.0310*eta_mu*C1w[4]-0.1262*C0w[2]-0.1221*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.0162*eta_mu*C1w[4]-0.1209*C0w[2]+0.0213*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0975*eta_mu*C1w[4]-0.1085*C0w[2]-0.1197*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[6]= -0.1191*C0w[2]*pow(eta_mu,6./23.)-0.2778*C0w[2]*pow(eta_mu,-12./23.)
	        +(0.0982*eta_mu*C1w[4]-0.5544*C0w[2]+0.1513*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.0634*eta_mu*C1w[4]+0.1915*C0w[2]-0.2497*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.3026*eta_mu*C1w[4]-0.2744*C0w[2]+0.3983*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(0.0358*eta_mu*C1w[4]+0.3568*C0w[2]+0.0440*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);

	C1b[7]= pow(eta_mu,39./23.)*C1w[7] *+8./3.*(pow(eta_mu,37./23.)-pow(eta_mu,39./23.))*C1w[8] *    
		+(297664./14283.*pow(eta_mu,16./23.)-7164416./357075.*pow(eta_mu,14./23.)+256868./14283.*pow(eta_mu,37./23.)-6698884./357075.*pow(eta_mu,39./23.))*C0w[8]
		+37208./4761.*(pow(eta_mu,39./23.)-pow(eta_mu,16./23.))*C0w[7]
		+(4661194./816831.*eta_mu*C1w[4]-17.3023*C0w[2]+14.8088*eta_mu*C1w[1]/15.)*pow(eta_mu,14./23.)
		+(-8516./2217.*eta_mu*C1w[4]+8.5027*C0w[2]-10.8090*eta_mu*C1w[1]/15.)*pow(eta_mu,16./23.)
		+ (4.5508*C0w[2]-0.8740*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)
		+ (0.7519*C0w[2]+0.4218*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
		+ (-1.9043*eta_mu*C1w[4]+2.0040*C0w[2]-2.9347*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+ (-0.1008*eta_mu*C1w[4]+0.7476*C0w[2]+0.3971*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+ (0.1216*eta_mu*C1w[4]-0.5385*C0w[2]+0.1600*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+ (0.0183*eta_mu*C1w[4]+0.0914*C0w[2]+0.0225*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
		
	return;	
}
