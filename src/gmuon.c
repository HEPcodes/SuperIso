#include "include.h"


double F1N(double x)
{
	if(x==1.) return 1.;
	return 2./pow(1.-x,4.)*(1.-6.*x+3.*x*x+2.*pow(x,3.)-6.*x*x*log(x));
}

/*--------------------------------------------------------------------*/

double F2N(double x)
{
	if(x==1.) return 1.;
	return 3./pow(1.-x,3.)*(1.-x*x+2.*x*log(x));
}

/*--------------------------------------------------------------------*/

double F1C(double x)
{
	if(x==1.) return 1.;
	return 2./pow(1.-x,4.)*(2.+3.*x-6.*x*x+pow(x,3.)+6.*x*log(x));
}

/*--------------------------------------------------------------------*/

double F2C(double x)
{
	if(x==1.) return 1.;
	return -3./2./pow(1.-x,3.)*(3.-4.*x+x*x+2.*log(x));
}

/*--------------------------------------------------------------------*/

double fPS(double x)
{
	if(x<0.25)
	{
		double y=sqrt(1.-4.*x);
		return 2.*x/y*(Li2(1.-(1.-y)/2./x)-Li2(1.-(1.+y)/2./x));
	}
	else if(x<2.5) return -1.365013496e-1-1.858455623*log(1.+x)-5.996763746e-1*log(1.+x)*log(1.+x)+4.390843985e-1*sqrt(x)*log(1.+x)-1.444359743e-1*x*log(1.+x)+3.852425143*sqrt(x);
	else if(x<100.)  return 4.304425955e-1+6.766323794e-2*log(1.+x)-1.584446296e-1*log(1.+x)*log(1.+x)-2.787080541e-1*sqrt(x)*log(1.+x)+1.557845370e-3*x*log(1.+x)+2.139180566*sqrt(x);
	else if(x<10000.)  return 2.025445594+9.960866255e-1*log(x)+1.122896720e-4*sqrt(x);
	else return 2.000835136+9.9992369e-1*log(x)+2.327105016e-7*sqrt(x);
}

/*--------------------------------------------------------------------*/

double fS(double x)
{
	return (2.*x-1.)*fPS(x)-2.*x*(2.+log(x));
}

/*--------------------------------------------------------------------*/

double fft(double x)
{
	return x/2.*(2.+log(x)-fPS(x));
}

/*--------------------------------------------------------------------*/

double muonf(double z)
{
	if(z<1.e-8) return 0.;
	
	double int1=0.;

	int ie;
	
	int ne=5000;
	double dx=0.5/ne;
	
	double x=0.5;
	int1+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	for(ie=1;ie<=ne-2;ie++)
	{	
		x-=dx;
		int1+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int1+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	int1*=dx;
	
	double int2=0.;

	x=dx;
	dx=dx/ne;
	int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x-=dx;
		int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	x=1.e-10;
	int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	int2*=dx;
		
	return int1+int2;
}

/*--------------------------------------------------------------------*/

double muong(double z)
{
	if(z<1.e-8) return 0.;
	
	double int1=0.;

	int ie;
	
	int ne=5000;
	double dx=0.5/ne;
	
	double x=0.5;
	int1+=1./(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	for(ie=1;ie<=ne-2;ie++)
	{	
		x-=dx;
		int1+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int1+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	
	int1*=dx;
	
	double int2=0.;

	x=dx;
	dx=dx/ne;
	int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x-=dx;
		int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	
	x=1.e-10;
	int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	int2*=dx;
		
	return int1+int2;
}

/*--------------------------------------------------------------------*/

double muon_gm2(struct parameters* param)
/* computes the muon anomalous magnetic moment a_mu */
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif	
	int ie,me,ke,je;

	if(param->THDM_model>0) 
	{
		double gmuon_2l=0.;
		double mass_f[9];
		double T3_f[9];
		double charge_f[9];
		double ncol_f[9];
		double sba=sin(atan(param->tan_beta)-param->alpha);
		double cba=cos(atan(param->tan_beta)-param->alpha);
		double v=246;
		double rho_f[9];
		double kappa_f[9];
		
		mass_f[0]=param->mass_u;
		
		for(ie=0;ie<=2;ie++) 
		{
			T3_f[ie]=0.5;
			charge_f[ie]=2./3.;
			ncol_f[ie]=3.;
		}
		
		mass_f[3]=param->mass_d;
		mass_f[4]=param->mass_s;
		
		for(ie=3;ie<=5;ie++) 
		{
			T3_f[ie]=-0.5;
			charge_f[ie]=-1./3.;
			ncol_f[ie]=3.;
		}
		
		mass_f[6]=param->mass_e;
		mass_f[7]=param->mass_mu;
		mass_f[8]=param->mass_tau_pole;
		
		for(ie=6;ie<=8;ie++) 
		{
			T3_f[ie]=-0.5;
			charge_f[ie]=-1.;
			ncol_f[ie]=1.;
		}
		
		/* h0 case */
		
		mass_f[1]=running_mass(param->mass_c,param->mass_c,param->mass_h0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[2]=running_mass(param->mtmt,param->mtmt,param->mass_h0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[5]=running_mass(param->mass_b,param->mass_b,param->mass_h0,param->mass_top_pole,param->mass_b_pole,param);
		
		for(ie=0;ie<=8;ie++) kappa_f[ie]=sqrt(2.)*mass_f[ie]/v;
		
		rho_f[0]=kappa_f[0]*param->lambda_u[1][1];
		rho_f[1]=kappa_f[1]*param->lambda_u[2][2];
		rho_f[2]=kappa_f[2]*param->lambda_u[3][3];
		rho_f[3]=kappa_f[3]*param->lambda_d[1][1];
		rho_f[4]=kappa_f[4]*param->lambda_d[2][2];
		rho_f[5]=kappa_f[5]*param->lambda_d[3][3];
		rho_f[6]=kappa_f[6]*param->lambda_l[1][1];
		rho_f[7]=kappa_f[7]*param->lambda_l[2][2];
		rho_f[8]=kappa_f[8]*param->lambda_l[3][3];
		
		for(ie=0;ie<=8;ie++) gmuon_2l+=
		ncol_f[ie]*charge_f[ie]*charge_f[ie]*mass_f[7]*mass_f[ie]/param->mass_h0/param->mass_h0
		*(-(kappa_f[ie]*sba+rho_f[ie]*cba)*(kappa_f[7]*sba+rho_f[7]*cba)/2.*muonf(pow(mass_f[ie]/param->mass_h0,2.)));
		
		/* H0 case */
		
		mass_f[1]=running_mass(param->mass_c,param->mass_c,param->mass_H0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[2]=running_mass(param->mtmt,param->mtmt,param->mass_H0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[5]=running_mass(param->mass_b,param->mass_b,param->mass_H0,param->mass_top_pole,param->mass_b_pole,param);
		
		for(ie=0;ie<=8;ie++) kappa_f[ie]=sqrt(2.)*mass_f[ie]/v;
		
		rho_f[0]=kappa_f[0]*param->lambda_u[1][1];
		rho_f[1]=kappa_f[1]*param->lambda_u[2][2];
		rho_f[2]=kappa_f[2]*param->lambda_u[3][3];
		rho_f[3]=kappa_f[3]*param->lambda_d[1][1];
		rho_f[4]=kappa_f[4]*param->lambda_d[2][2];
		rho_f[5]=kappa_f[5]*param->lambda_d[3][3];
		rho_f[6]=kappa_f[6]*param->lambda_l[1][1];
		rho_f[7]=kappa_f[7]*param->lambda_l[2][2];
		rho_f[8]=kappa_f[8]*param->lambda_l[3][3];
		
		for(ie=0;ie<=8;ie++) gmuon_2l+=
		ncol_f[ie]*charge_f[ie]*charge_f[ie]*mass_f[7]*mass_f[ie]/param->mass_H0/param->mass_H0
		*(-(kappa_f[ie]*cba-rho_f[ie]*sba)*(kappa_f[7]*cba-rho_f[7]*sba)/2.*muonf(pow(mass_f[ie]/param->mass_H0,2.)));
		
		
		/* A0 case */
		
		mass_f[1]=running_mass(param->mass_c,param->mass_c,param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[2]=running_mass(param->mtmt,param->mtmt,param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[5]=running_mass(param->mass_b,param->mass_b,param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);
		
		for(ie=0;ie<=8;ie++) kappa_f[ie]=sqrt(2.)*mass_f[ie]/v;
		
		rho_f[0]=kappa_f[0]*param->lambda_u[1][1];
		rho_f[1]=kappa_f[1]*param->lambda_u[2][2];
		rho_f[2]=kappa_f[2]*param->lambda_u[3][3];
		rho_f[3]=kappa_f[3]*param->lambda_d[1][1];
		rho_f[4]=kappa_f[4]*param->lambda_d[2][2];
		rho_f[5]=kappa_f[5]*param->lambda_d[3][3];
		rho_f[6]=kappa_f[6]*param->lambda_l[1][1];
		rho_f[7]=kappa_f[7]*param->lambda_l[2][2];
		rho_f[8]=kappa_f[8]*param->lambda_l[3][3];
		
		for(ie=0;ie<=8;ie++) gmuon_2l+=
		ncol_f[ie]*charge_f[ie]*charge_f[ie]*mass_f[7]*mass_f[ie]/param->mass_A0/param->mass_A0
		*(-2.*T3_f[ie]*rho_f[ie]*rho_f[7]/2.*muong(pow(mass_f[ie]/param->mass_A0,2.)));
		
		return gmuon_2l/param->inv_alpha_em/4./pi/pi/pi;
	
	}

	/* SUSY */

	double mass_smu[3],mass_charg[3];
	double n_L[6][3],n_R[6][3],xim[6][3],xk[3],c_L[3],c_R[3],X[3][3];	
	int nb_neut;
	if(param->mass_neut[5]==0.) nb_neut=4; else nb_neut=5;


	mass_charg[1]=param->mass_cha1;
	mass_charg[2]=param->mass_cha2;


	double M2smu11=pow(param->mass_mul,2.);
	double M2smu22=pow(param->mass_mur,2.);
	double M2smu12=(param->A_mu-param->mu_Q*param->tan_beta)*param->mass_mu;

	mass_smu[1]=sqrt((M2smu11+M2smu22-sqrt(pow(M2smu11-M2smu22,2.)+4.*M2smu12*M2smu12))/2.);
	mass_smu[2]=sqrt((M2smu11+M2smu22+sqrt(pow(M2smu11-M2smu22,2.)+4.*M2smu12*M2smu12))/2.);

	X[1][1]=cos(atan2(-2.*M2smu12,-M2smu11+M2smu22)/2.);
	X[1][2]=sin(atan2(-2.*M2smu12,-M2smu11+M2smu22)/2.);
	X[2][1]=-X[1][2];
	X[2][2]=X[1][1];
	
	double ymu=param->g2*param->mass_mu/sqrt(2.)/param->mass_W/cos(atan(param->tan_beta));
	
	for(ie=1;ie<=nb_neut;ie++) for(me=1;me<=2;me++)
	{
		n_L[ie][me]=1./sqrt(2.)*(param->gp*param->neut_mix[ie][1]+param->g2*param->neut_mix[ie][2])*X[me][1]
		-ymu*param->neut_mix[ie][3]*X[me][2];
		n_R[ie][me]=sqrt(2.)*param->gp*param->neut_mix[ie][1]*X[me][2]+ymu*param->neut_mix[ie][3]*X[me][1];
		xim[ie][me]=pow(param->mass_neut[ie]/mass_smu[me],2.);
	}
	
	for(ke=1;ke<=2;ke++)
	{
		c_L[ke]=-param->g2*param->charg_Vmix[ke][1];
		c_R[ke]=ymu*param->charg_Umix[ke][2];
		xk[ke]=pow(mass_charg[ke]/param->mass_numl,2.);
	}
		
	double a_neut=0.;
	for(ie=1;ie<=nb_neut;ie++) for(me=1;me<=2;me++)
	{
		a_neut+=
		-param->mass_mu/12./pow(mass_smu[me],2.)*(pow(n_L[ie][me],2.)+pow(n_R[ie][me],2.))*F1N(xim[ie][me])
		+param->mass_neut[ie]/3./pow(mass_smu[me],2.)*n_L[ie][me]*n_R[ie][me]*F2N(xim[ie][me]);
	}
	a_neut*=param->mass_mu/16./pi/pi;
	
		
	double a_charg=0.;
	for(ke=1;ke<=2;ke++)
	{
		a_charg+=
		param->mass_mu/12./pow(param->mass_numl,2.)*(pow(c_L[ke],2.)+pow(c_R[ke],2.))*F1C(xk[ke])
		+2.*mass_charg[ke]/3./pow(param->mass_numl,2.)*c_L[ke]*c_R[ke]*F2C(xk[ke]);
	}
	a_charg*=param->mass_mu/16./pi/pi;

	double a_SUSYQED=(a_neut+a_charg)*(1.-4./param->inv_alpha_em/pi*log(param->MSOFT_Q/param->mass_mu));
	

	double lmu[4],lcharg[4][3],lst[3][3],lsb[3][3],lstau[3][3];
	double mass_stop[3],mass_sbot[3],mass_stau[3],mass_H[3];
	
	mass_H[1]=param->mass_h0;
	mass_H[2]=param->mass_H0;

	mass_stop[1]=param->mass_t1;
	mass_stop[2]=param->mass_t2;

	mass_sbot[1]=param->mass_b1;
	mass_sbot[2]=param->mass_b2;

	mass_stau[1]=param->mass_tau1;
	mass_stau[2]=param->mass_tau2;

	double sa=sin(param->alpha);
	double ca=cos(param->alpha);
	double sb=sin(atan(param->tan_beta));
	double cb=cos(atan(param->tan_beta));
	
	double sw=sin(atan(param->gp/param->g2));

	lmu[1]=-sa/cb;
	lmu[2]=ca/cb;
	lmu[3]=param->tan_beta;
	
	for(ie=1;ie<=2;ie++)
	{
		lcharg[1][ie] = sqrt(2.)*param->mass_W/mass_charg[ie]*(param->charg_Umix[ie][1]*param->charg_Vmix[ie][2]*ca - param->charg_Umix[ie][2]*param->charg_Vmix[ie][1]*sa);

		lcharg[2][ie] = sqrt(2.)*param->mass_W/mass_charg[ie]*(param->charg_Umix[ie][1]*param->charg_Vmix[ie][2]*sa + param->charg_Umix[ie][2]*param->charg_Vmix[ie][1]*ca);
		
		lcharg[3][ie] = -sqrt(2.)*param->mass_W/mass_charg[ie]*(param->charg_Umix[ie][1]*param->charg_Vmix[ie][2]*cb + param->charg_Umix[ie][2]*param->charg_Vmix[ie][1]*sb);
		
		lst[1][ie] = 2.*param->mtmt/mass_stop[ie]/mass_stop[ie]/sb*(param->mu_Q*sa+param->A_t*ca)*param->stop_mix[ie][1]*param->stop_mix[ie][2];
		
		lst[2][ie] = 2.*param->mtmt/mass_stop[ie]/mass_stop[ie]/sb*(-param->mu_Q*ca+param->A_t*sa)*param->stop_mix[ie][1]*param->stop_mix[ie][2];

		lsb[1][ie] = 2.*param->mass_b/mass_sbot[ie]/mass_sbot[ie]/cb*(-param->mu_Q*ca-param->A_b*sa)*param->sbot_mix[ie][1]*param->sbot_mix[ie][2];
		
		lsb[2][ie] = 2.*param->mass_b/mass_sbot[ie]/mass_sbot[ie]/cb*(-param->mu_Q*sa+param->A_b*ca)*param->sbot_mix[ie][1]*param->sbot_mix[ie][2];
		
		lstau[1][ie] = 2.*param->mass_tau_pole/mass_stau[ie]/mass_stau[ie]/cb*(-param->mu_Q*ca-param->A_tau*sa)*param->stau_mix[ie][1]*param->stau_mix[ie][2];
		
		lstau[2][ie] = 2.*param->mass_tau_pole/mass_stau[ie]/mass_stau[ie]/cb*(-param->mu_Q*sa+param->A_tau*ca)*param->stau_mix[ie][1]*param->stau_mix[ie][2];
	}
	
	double a_charH=0.;
	double a_sfH=0.;
	
	for(ke=1;ke<=2;ke++)
	{
		for(je=1;je<=2;je++) a_charH+=lmu[je]*lcharg[je][ke] *fS(mass_charg[ke]*mass_charg[ke]/mass_H[je]/mass_H[je]);
		
		a_charH += lmu[3]*lcharg[3][ke]*fPS(mass_charg[ke]*mass_charg[ke]/param->mass_A0/param->mass_A0);	
	}
	double const1=param->mass_mu*param->mass_mu/8./pi/pi/param->inv_alpha_em/param->inv_alpha_em/param->mass_W/param->mass_W/sw/sw;
	a_charH*=const1;

	for(ie=1;ie<=2;ie++) for(ke=1;ke<=2;ke++)
	{
		a_sfH+=4./3.*lmu[ke]*lst[ke][ie]*fft(mass_stop[ie]*mass_stop[ie]/mass_H[ke]/mass_H[ke]);
		a_sfH+=1./3.*lmu[ke]*lsb[ke][ie]*fft(mass_sbot[ie]*mass_sbot[ie]/mass_H[ke]/mass_H[ke]);
		a_sfH+=lmu[ke]*lstau[ke][ie]*fft(mass_stau[ie]*mass_stau[ie]/mass_H[ke]/mass_H[ke]);
	}
	a_sfH*=const1;

	return a_SUSYQED+a_charH+a_sfH;
}

/*--------------------------------------------------------------------*/

double muon_gm2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating a_mu */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return muon_gm2(&param);
}

