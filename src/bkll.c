#include "include.h"

/*----------------------------------------------------------------------*/

double complex h_bkll(double q2, double mq, double mu)
{
	if(mq==0.) return 4./9.*(2./3.+I*pi+log(mu*mu/q2));
	
	double z=4.*mq*mq/q2;
	
	if(z>1.) return -4./9.*(log(mq*mq/mu/mu)-2./3.-z)
	-4./9.*(2.+z)*sqrt(z-1.)*atan(1./sqrt(z-1.));
	
	else return -4./9.*(log(mq*mq/mu/mu)-2./3.-z)
	-4./9.*(2.+z)*sqrt(1.-z)*(log((1.+sqrt(1.-z))/sqrt(z))-I*pi/2.);
}

/*----------------------------------------------------------------------*/

double phi_Kstar(double u, double a1, double a2)
{
	double x=2.*u-1.;
	double C1=3.*x;
	double C2=-1.5+15./2.*x*x;

	return 6.*u*(1.-u)*(1.+a1*C1+a2*C2);
}

/*----------------------------------------------------------------------*/

double complex B0_bkll(double s, double mq)
{
	double epsilon=1.e-10;
	return -2.*csqrt(4.*(mq*mq-I*epsilon)/s-1.)*catan(1./csqrt(4.*(mq*mq-I*epsilon)/s-1.));
}

/*----------------------------------------------------------------------*/

double complex L1_bkll(double complex x)
{
	return clog((x-1.)/x)*clog(1.-x)-pi*pi/6.+CLi2(x/(x-1.));
}

/*----------------------------------------------------------------------*/

double complex I1_bkll(double u, double mq, double q2, struct parameters* param)
{
	if(mq==0.) return 1.;
	
	double epsilon=1.e-10;

	double complex xp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd+u*q2));
	double complex xm=0.5-csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd+u*q2));
	double complex yp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/q2);
	double complex ym=0.5-csqrt(0.25-(mq*mq-I*epsilon)/q2);

	return 1.+2.*mq*mq/(1.-u)/(param->m_Bd*param->m_Bd-q2)*(L1_bkll(xp)+L1_bkll(xm)-L1_bkll(yp)-L1_bkll(ym));
}

/*----------------------------------------------------------------------*/

double complex tperp_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param)
{
	return 2.*param->m_Bd/(1.-u)/E_Kstar*I1_bkll(u,mq,q2,param)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B0_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mq)-B0_bkll(q2,mq));
}

/*----------------------------------------------------------------------*/

double complex tpar_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param)
{
	return 2.*param->m_Bd/(1.-u)/E_Kstar*I1_bkll(u,mq,q2,param)+((1.-u)*param->m_Bd*param->m_Bd+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B0_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mq)-B0_bkll(q2,mq));
}

/*----------------------------------------------------------------------*/

double complex F27u(double shat, double l)
{
	double z=4./shat;

	double complex A=
	208./243.*l+4.*shat/27./(1.-shat)*(Li2(shat)+log(shat)*log(1.-shat))
	+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*log(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
	-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(cabs(z-1.)))+9.*shat*shat*shat*log(shat)*log(shat)+18.*pi*I*shat*(1.-2.*shat)*log(shat))
	+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
	
	return -6.*A;
}

/*----------------------------------------------------------------------*/

double complex F19u(double shat, double l)
{
	double z=4./shat;
	
	double complex x1=0.5+0.5*I*csqrt(z-1.);
	double complex x2=0.5-0.5*I*csqrt(z-1.);
	double complex x3=0.5+0.5*I/csqrt(z-1.);
	double complex x4=0.5-0.5*I/csqrt(z-1.);

	double complex B=
	8./243./shat*(-2.*(4.-34.*shat-17.*pi*I*shat)*l+8.*shat*pow(2.*l,2.)-17.*shat*log(shat)*2.*l)
	+(2.+shat)*csqrt(z-1.)/729./shat*(48.*2.*l*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
	-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
	-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
	+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
	+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
	-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
	+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
	double complex C=-16./81.*(log(shat)-2.*l)+428./243.-64./27.*zeta3+16./81.*pi*I;
		
	return B+4.*C;
}

/*----------------------------------------------------------------------*/

double complex F29u(double shat, double l)
{
	double z=4./shat;

	double complex x1=0.5+0.5*I*csqrt(z-1.);
	double complex x2=0.5-0.5*I*csqrt(z-1.);
	double complex x3=0.5+0.5*I/csqrt(z-1.);
	double complex x4=0.5-0.5*I/csqrt(z-1.);

	double complex B=
	8./243./shat*(-2.*(4.-34.*shat-17.*pi*I*shat)*l+8.*shat*pow(2.*l,2.)-17.*shat*log(shat)*2.*l)
	+(2.+shat)*csqrt(z-1.)/729./shat*(48.*2.*l*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
	-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
	-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
	+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
	+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
	-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
	+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
	double complex C=-16./81.*(log(shat)-2.*l)+428./243.-64./27.*zeta3+16./81.*pi*I;	
	return -6.*B+3.*C;
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarmumu_dq2(double q2, double obs[][3], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double mc=mc_pole_1loop(param);
	double mbpole=mb_pole_1loop(param);
	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);

	int ie,je;
	
	double beta_l=sqrt(1.-4.*param->mass_mu*param->mass_mu/q2);

	double alpha_em=1./133.;

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(2.,param->mass_top_pole,param->mass_b_pole,param);

	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(param->m_Bd*param->m_Bd+param->m_Kstar*param->m_Kstar-q2)/2./param->m_Bd;

	int nf=5;
	double f_K_perp=param->f_K_perp;
	f_K_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_K_par=param->f_K_par;
	
	double complex ALperp=0.;
	double complex ARperp=0.;
	double complex ALpar=0.;
	double complex ARpar=0.;
	double complex AL0=0.;
	double complex AR0=0.;
	double complex At=0.;
	double complex AS=0.;

	double complex ALperp_bar=0.;
	double complex ARperp_bar=0.;
	double complex ALpar_bar=0.;
	double complex ARpar_bar=0.;
	double complex AL0_bar=0.;
	double complex AR0_bar=0.;
	double complex At_bar=0.;
	double complex AS_bar=0.;

	double complex ALperp_high=0.;
	double complex ARperp_high=0.;
	double complex ALpar_high=0.;
	double complex ARpar_high=0.;
	double complex AL0_high=0.;
	double complex AR0_high=0.;
	double complex At_high=0.;
	double complex AS_high=0.;

	double complex ALperp_bar_high=0.;
	double complex ARperp_bar_high=0.;
	double complex ALpar_bar_high=0.;
	double complex ARpar_bar_high=0.;
	double complex AL0_bar_high=0.;
	double complex AR0_bar_high=0.;
	double complex At_bar_high=0.;
	double complex AS_bar_high=0.;
	
	double V=0.;
	double A1=0.;
	double A2=0.;
		
	double xi_par=0.;
	double xi_perp=0.;	

	if(q2<14.)
	{
	 	double C7eff=Cmub[7];
		double C8eff=Cmub[8];
		double C9=Cmub[9];
		double C10=Cmub[10];
	
		double C7effp=Cpb[7];
		double C9p=Cpb[9];
		double C10p=Cpb[10];
	
		double C1bar=Cmub[1]/2.;
		double C2bar=Cmub[2]-Cmub[1]/6.;
		double C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
		double C4bar=Cmub[4]/2.+8.*Cmub[6];
		double C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
		double C6bar=Cmub[4]/2.+2.*Cmub[6];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi;

		double tau_plus=pow(param->m_Bd+param->m_Kstar,2.);
		double tau_minus=pow(param->m_Bd-param->m_Kstar,2.);
		double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
		
		double V0=0.36;
		double b1V=-4.8;
		double A00=0.29;
		double b1A0=-18.2*0.3;
		V=V0/(1.-q2/param->m_Bs/param->m_Bs)*(1.+b1V*(z_q2tau0-z_0tau0+0.5*(z_q2tau0*z_q2tau0-z_0tau0*z_0tau0)));
		double A0=A00/(1.-q2/param->m_Bs/param->m_Bs)*(1.+b1A0*(z_q2tau0-z_0tau0+0.5*(z_q2tau0*z_q2tau0-z_0tau0*z_0tau0)));
		A1=2.*param->m_Bd*E_Kstar/pow(param->m_Bd+param->m_Kstar,2.)*V;
		A2=((param->m_Bd+param->m_Kstar)*A1/2./E_Kstar-(param->m_Kstar/E_Kstar)*A0)*param->m_Bd/(param->m_Bd-param->m_Kstar);
		
		xi_par=(param->m_Bd+param->m_Kstar)/2./E_Kstar*A1-(param->m_Bd-param->m_Kstar)/param->m_Bd*A2;
		xi_perp=param->m_Bd/(param->m_Bd+param->m_Kstar)*V;
		
		double complex h_mc=h_bkll(q2,mc,mu_b);
		double complex h_mb=h_bkll(q2,mbpole,mu_b);
		double complex h_0=h_bkll(q2,0.,mu_b);

		double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
		+h_mc*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
		+h_mb*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
		+h_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

		double Yu=(h_mc-h_0)*(4./3.*Cmub[1]+Cmub[2]);

		double complex Cperpp0=C7eff+C7effp+q2/2./mb/param->m_Bd*Y;
		double complex Cperpm0=C7eff-C7effp+q2/2./mb/param->m_Bd*Y;
		double complex Cparp0=-(C7eff+C7effp)-param->m_Bd/2./mb*Y;
		double complex Cparm0=-(C7eff-C7effp)-param->m_Bd/2./mb*Y;

		double complex Cperpp0u=q2/2./mb/param->m_Bd*Yu;
		double complex Cperpm0u=q2/2./mb/param->m_Bd*Yu;
		double complex Cparp0u=-param->m_Bd/2./mb*Yu;
		double complex Cparm0u=-param->m_Bd/2./mb*Yu;

		double logb=log(mu_b/mb);
		double DeltaM=-6.*logb-4.*(1.-mu_f/mb);
		double L=-(mb*mb-q2)/q2*log(1.-q2/mb/mb);
	
		double shat=q2/mb/mb;

		double mchat=mc/mb;
		double z=mchat*mchat;	

		double complex Cperppf=(C7eff+C7effp)*(-2.*logb-L+DeltaM);
		double complex Cperpmf=(C7eff-C7effp)*(-2.*logb-L+DeltaM);
		double complex Cparpf=-(C7eff+C7effp)*(-2.*logb+2.*L+DeltaM);
		double complex Cparmf=-(C7eff-C7effp)*(-2.*logb+2.*L+DeltaM);

		double complex F27=F27_bsll(shat,z,logb);
		double complex F87=F87_bsll(shat,logb);
		double complex F29=F29_bkll(shat,z,logb);
		double complex F19=F19_bkll(shat,z,logb);
		double complex F89=F89_bsll(shat);
		double complex F27_u=F27u(shat,logb);
		double complex F29_u=F29u(shat,logb);
		double complex F19_u=F19u(shat,logb);

		double complex Cperpnf=(-C2bar*F27-C8eff*F87
-q2/2./mb/param->m_Bd*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		double complex Cparnf=(C2bar*F27+C8eff*F87
+param->m_Bd/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;

		double complex Cperpnfu=(-C2bar*(F27+F27_u)
		-q2/2./mb/param->m_Bd*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cparnfu=(C2bar*(F27+F27_u)
		+param->m_Bd/2./mb*(C2bar*(F29+F29_u)+2.*C1bar*(F19+F19_u+1./6.*(F29+F29_u))))/4.*3.;
		
		double complex Cperpp1=Cperppf+Cperpnf;
		double complex Cperpm1=Cperpmf+Cperpnf;
		double complex Cparp1=Cparpf+Cparnf;
		double complex Cparm1=Cparmf+Cparnf;
	
		double complex Cperpp=Cperpp0+alphas_mub*4./3./4./pi*Cperpp1;
		double complex Cperpm=Cperpm0+alphas_mub*4./3./4./pi*Cperpm1;
		double complex Cparp=Cparp0+alphas_mub*4./3./4./pi*Cparp1;
		double complex Cparm=Cparm0+alphas_mub*4./3./4./pi*Cparm1;
				
		double complex Cperppu=Cperpp0u+alphas_mub*4./3./4./pi*Cperpnfu;
		double complex Cperpmu=Cperpm0u+alphas_mub*4./3./4./pi*Cperpnfu;
		double complex Cparpu=Cparp0u+alphas_mub*4./3./4./pi*Cparnfu;
		double complex Cparmu=Cparm0u+alphas_mub*4./3./4./pi*Cparnfu;
		
		
		double Xi_perp=1.;
		double Xi_par=param->m_Kstar/E_Kstar;
	
		double eu=2./3.;
		double ed=-1./3.;
		double eq=-1./3.;
	
		double a1perp=param->a1perp;
		double a2perp=param->a2perp;
		double a1par=param->a1par;
		double a2par=param->a2par;

		a1perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
		a2perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

		a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_perppp,int_perppm,int_perpmp,int_perpmm;
		double complex int_parpp,int_parpm,int_parmp,int_parmm;
		double complex int_perppu,int_parpu,int_parmu;
		double complex Tperppp0,Tperpppf,Tperpppnf,Tperppp;
		double complex Tperppm0,Tperppmf,Tperppmnf,Tperppm;
		double complex Tperpmp0,Tperpmpf,Tperpmpnf,Tperpmp;
		double complex Tperpmm0,Tperpmmf,Tperpmmnf,Tperpmm;
		double complex Tparpp0,Tparppf,Tparppnf,Tparpp;
		double complex Tparpm0,Tparpmf,Tparpmnf,Tparpm;
		double complex Tparmp0,Tparmpf,Tparmpnf,Tparmp;
		double complex Tparmm0,Tparmmf,Tparmmnf,Tparmm;
		double complex Tperppnfu,Tperppu;
		double complex Tparp0u,Tparpnfu,Tparpu;
		double complex Tparm0u,Tparmnfu,Tparmu;
		
		int_perppp=int_perppm=int_perpmp=int_perpmm=0.;
		int_perppu=0.;
		int_parpp=int_parpm=int_parmp=int_parmm=0.;
		int_parpu=int_parmu=0.;


		double lambda_Bp=param->lambda_Bp;
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(param->m_Bd-mb)/3.;
		double complex lambda_Bm=1./(exp(-q2/param->m_Bd/omega0)/omega0*(-Ei(q2/param->m_Bd/omega0)+I*pi));

		double phiKstar_perp,phiKstar_par;
		double complex tperp_mb,tperp_mc,tperp_0;
		double complex tpar_mb,tpar_mc,tpar_0;

		double complex integ3=0;

		double complex Fperp=0.;
		double complex Xperp=0.;
		double x;
		double complex integ4=0.;
		double complex FV;
		double complex integ4u=0.;
		double complex FVu;
		double complex integ5=0.;
		double complex integ5u=0.;

		double zeta3A=param->zeta3A;
		double zeta3V=param->zeta3V;
		double wA10=param->wA10;
		double deltatp=param->deltatp;
		double deltatm=param->deltatm;
		
		int n1=10;
		int n1sav=n1;
		for(ie=0;ie<=n1;ie++)
		{
			u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=(double)(ie-1)/n1;n1*=2;}

			/* Tperp */	
				
			Tperppp0=Tperpmp0=0.;
			Tperpppf=(C7eff+C7effp)*2.*param->m_Bd/E_Kstar/(1.-u);
			Tperpmpf=(C7eff-C7effp)*2.*param->m_Bd/E_Kstar/(1.-u);

			Tperppm0=Tperpmm0=0.;
			Tperppmf=Tperpmmf=0.;
			
			phiKstar_perp=phi_Kstar(u,a1perp,a2perp);
			tperp_mc=tperp_bkll(u,mc,q2,E_Kstar,param);
			tperp_mb=tperp_bkll(u,mb,q2,E_Kstar,param);
			tperp_0=tperp_bkll(u,0.,q2,E_Kstar,param);

			Tperpppnf=Tperpmpnf=-4.*ed*C8eff/(u+(1.-u)*q2/param->m_Bd/param->m_Bd)
			+param->m_Bd/2./mb*(eu*tperp_mc*(C2bar+C4bar-C6bar)
			+ed*tperp_mb*(C3bar+C4bar-C6bar-4.*mb/param->m_Bd*C5bar)
			+ed*tperp_0*C3bar);
				
			Tperppmnf=Tperpmmnf=0.;

			Tperppnfu=param->m_Bd/2./mb*eu*(tperp_mc-tperp_0)*(Cmub[2]-Cmub[1]/6.);
		
			Tperppp=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf+Tperpppnf); 
			Tperppm=Tperppm0+alphas_muf*4./3./4./pi*(Tperppmf+Tperppmnf);
			Tperpmp=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf+Tperpmpnf);
			Tperpmm=Tperpmm0+alphas_muf*4./3./4./pi*(Tperpmmf+Tperpmmnf);
	
			Tperppu=alphas_muf*4./3./4./pi*Tperppnfu;
	
			int_perppp+=phiKstar_perp*Tperppp/n1/lambda_Bp; 
			int_perppm+=phiKstar_perp*Tperppm/n1/lambda_Bm;
			int_perpmp+=phiKstar_perp*Tperpmp/n1/lambda_Bp;
			int_perpmm+=phiKstar_perp*Tperpmm/n1/lambda_Bm;

			int_perppu+=phiKstar_perp*Tperppu/n1/lambda_Bp;


			/* Tpar */		

			phiKstar_par=phi_Kstar(u,a1par,a2par);
			tpar_mc=tpar_bkll(u,mc,q2,E_Kstar,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_Kstar,param);
			tpar_0=tpar_bkll(u,0.,q2,E_Kstar,param);
		
			Tparpp0=Tparmp0=0.;
		
			Tparppf=(C7eff+C7effp)*4.*param->m_Bd/E_Kstar/(1.-u);
			Tparmpf=(C7eff-C7effp)*4.*param->m_Bd/E_Kstar/(1.-u);

			Tparppnf=Tparmpnf=param->m_Bd/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
				
			Tparpnfu=param->m_Bd/mb*eu*(tpar_mc-tpar_0)*(Cmub[2]-Cmub[1]/6.);
	
			Tparpu=alphas_muf*4./3./4./pi*Tparpnfu;
	
	
			Tparpm0=Tparmm0=-eq*4.*param->m_Bd/mb*(C3bar+3.*C4bar);
			
			Tparm0u=-eq*4.*param->m_Bd/mb*(-3*Cmub[2]);
		
			Tparpmf=Tparmmf=0.;

			h_mc=h_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,0.,mu_b);

			Tparpmnf=Tparmmnf=eq*(8.*C8eff/((1.-u)+u*q2/param->m_Bd/param->m_Bd)
			+6.*param->m_Bd/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf);
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);
			Tparmp=Tparmp0+alphas_muf*4./3./4./pi*(Tparmpf+Tparmpnf);
			Tparmm=Tparmm0+alphas_muf*4./3./4./pi*(Tparmmf+Tparmmnf);
			
			Tparmnfu=eq*(6.*param->m_Bd/mb*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.));

			Tparpu=Tparp0u+alphas_muf*4./3./4./pi*Tparpnfu;
	
			Tparmu=Tparm0u+alphas_muf*4./3./4./pi*Tparmnfu;


			int_parpp+=(phiKstar_par*Tparpp/lambda_Bp)/n1;
			int_parpm+=(phiKstar_par*Tparpm/lambda_Bm)/n1;
			
			int_parmp+=(phiKstar_par*Tparmp/lambda_Bp)/n1;
			int_parmm+=(phiKstar_par*Tparmm/lambda_Bm)/n1;
			
			int_parpu+=(phiKstar_par*Tparpu/lambda_Bp)/n1;
			
			int_parmu+=(phiKstar_par*Tparmu/lambda_Bm)/n1;


			integ3+=phiKstar_perp/((1.-u)+u*shat)/n1;

			x=(1.-u)*param->m_B*param->m_B+u*q2;
			
			h_mc=h_bkll(x,mc,mu_b);
			h_mb=h_bkll(x,mbpole,mu_b);
			h_0=h_bkll(x,0.,mu_b);
			
			FV=3./4.*(h_mc*(C2bar+C4bar+C6bar)+h_mb*(C3bar+C4bar+C6bar)+h_0*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));

			FVu=3./4.*(h_mc-h_0)*(Cmub[2]-Cmub[1]/6.);
		

			integ4+=phiKstar_perp/((1.-u)+u*shat)*FV/n1;
			integ4u+=phiKstar_perp/((1.-u)+u*shat)*FVu/n1;
		
			Fperp+=phiKstar_perp/((1.-u)+u*shat)/3./n1;
			Xperp+=(u<=1.-0.5/param->m_B)*phiKstar_perp/pow((1.-u)+u*shat,2.)/3.*(0.5/param->m_B)/n1;
		
			integ5+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;

			integ5u+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FVu/n1;

			if(ie==0||ie==n1sav) n1=n1sav;
		}
		
		/* Tau_perp */		

		double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_B*f_K_perp/param->m_Bd*Xi_perp*(int_perppp+int_perppm);
		double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_B*f_K_perp/param->m_Bd*Xi_perp*(int_perpmp+int_perpmm);
		
		double complex Tauperppu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_perp*Cperppu+pi*pi/3.*param->f_B*f_K_perp/param->m_Bd*Xi_perp*int_perppu);
		double complex Tauperppu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_perp*Cperppu+pi*pi/3.*param->f_B*f_K_perp/param->m_Bd*Xi_perp*int_perppu);
		
		double complex Tauperpmu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*xi_perp*Cperpmu;
		double complex Tauperpmu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*xi_perp*Cperpmu;

		
		/* Tau_par */		

		double complex Tauparp=xi_par*Cparp+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*(int_parpp+int_parpm);
		
		double complex Tauparm=xi_par*Cparm+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*(int_parmp+int_parmm);
		
		double complex Tauparpu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_par*Cparpu+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*int_parpu);
				
		double complex Tauparmu=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(xi_par*Cparmu+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*int_parmu);
				
		double complex Tauparmu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_par*Cparmu+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*int_parmu);
		
		double complex Tauparpu_bar=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(xi_par*Cparpu+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*int_parpu);
		
				
		double complex DeltaTauperpWA=-eq*4.*pi*pi/3.*param->f_B*f_K_perp/mb/param->m_Bd*(Cmub[3]+4./3.*Cmub[4]+4.*Cmub[5]+16./3.*Cmub[6])*integ3
		+eq*2.*pi*pi/3.*param->f_B*f_K_par/mb/param->m_Bd*param->m_Kstar/(1.-shat)/lambda_Bp*(Cmub[3]+4./3.*Cmub[4]+48./3.*Cmub[5]+64./3.*Cmub[6]);

		double complex DeltaTauperpuWA=-eq*2.*pi*pi/3.*param->f_B*f_K_par/mb/param->m_Bd*param->m_Kstar/(1.-shat)/lambda_Bp*3.*Cmub[2];


		double rho=0.;
		double phi=0.;
		Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;


		double complex DeltaTauperpHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/param->m_Bd*(12.*C8eff*mb/param->m_Bd*f_K_perp*Xperp
		+8.*f_K_perp*integ4-4.*param->m_Kstar*f_K_par/(1.-shat)/lambda_Bp*integ5);
		
		double complex DeltaTauperpuHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/param->m_Bd*
		(8.*f_K_perp*integ4u-4.*param->m_Kstar*f_K_par/(1.-shat)/lambda_Bp*integ5u);

		Tauperpp+=DeltaTauperpWA+DeltaTauperpHSA;
		Tauperpm+=DeltaTauperpWA+DeltaTauperpHSA;

		Tauperppu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
		Tauperppu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);	
		
		Tauperpmu+=param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts)*(DeltaTauperpuWA+DeltaTauperpuHSA);
		Tauperpmu_bar+=param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb)*(DeltaTauperpuWA+DeltaTauperpuHSA);		

		double lambda=pow(param->m_Bd,4.)+pow(param->m_Kstar,4.)+q2*q2-2.*(param->m_Bd*param->m_Bd*param->m_Kstar*param->m_Kstar+param->m_Kstar*param->m_Kstar*q2+param->m_Bd*param->m_Bd*q2);
	
		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/param->m_Bd*shat*sqrt(lambda)*beta_l);
		
		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/param->m_Bd*shat*sqrt(lambda)*beta_l);
 
		ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(param->m_Bd+param->m_Kstar)+2.*mb/q2*(Tauperpp+Tauperppu));
		
		ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(param->m_Bd+param->m_Kstar)+2.*mb/q2*(Tauperpp+Tauperppu));
		
		ALpar=-N*sqrt(2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(((C9-C9p)-(C10-C10p))*A1/(param->m_Bd-param->m_Kstar)+4.*mb/param->m_Bd*E_Kstar/q2*(Tauperpm+Tauperpmu));
	
		ARpar=-N*sqrt(2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(((C9-C9p)+(C10-C10p))*A1/(param->m_Bd-param->m_Kstar)+4.*mb/param->m_Bd*E_Kstar/q2*(Tauperpm+Tauperpmu));
	
		AL0=-N/2./param->m_Kstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*((param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*(param->m_Bd+param->m_Kstar)*A1-lambda*A2/(param->m_Bd+param->m_Kstar))+2.*mb*(2.*E_Kstar/param->m_Bd*(param->m_Bd*param->m_Bd+3.*param->m_Kstar*param->m_Kstar-q2)*(Tauperpm+Tauperpmu)-lambda/(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(Tauperpm+Tauparm+Tauperpmu+Tauparmu)));
	
		AR0=-N/2./param->m_Kstar/sqrt(q2)*(((C9-C9p)+(C10-C10p))*((param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*(param->m_Bd+param->m_Kstar)*A1-lambda*A2/(param->m_Bd+param->m_Kstar))+2.*mb*(2.*E_Kstar/param->m_Bd*(param->m_Bd*param->m_Bd+3.*param->m_Kstar*param->m_Kstar-q2)*(Tauperpm+Tauperpmu)-lambda/(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(Tauperpm+Tauparm+Tauperpmu+Tauparmu)));

		At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/param->mass_mu*(CQ2-CQ2p)/(mb_mub+param->mass_s))*E_Kstar/param->m_Kstar*xi_par;
	
		AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*E_Kstar/param->m_Kstar*xi_par;
		
		
		ALperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(param->m_Bd+param->m_Kstar)+2.*mb/q2*(Tauperpp+Tauperppu_bar));
	
		ARperp_bar=Nbar*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(param->m_Bd+param->m_Kstar)+2.*mb/q2*(Tauperpp+Tauperppu_bar));
	
		ALpar_bar=-Nbar*sqrt(2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(((C9-C9p)-(C10-C10p))*A1/(param->m_Bd-param->m_Kstar)+4.*mb/param->m_Bd*E_Kstar/q2*(Tauperpm+Tauperpmu_bar));
	
		ARpar_bar=-Nbar*sqrt(2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(((C9-C9p)+(C10-C10p))*A1/(param->m_Bd-param->m_Kstar)+4.*mb/param->m_Bd*E_Kstar/q2*(Tauperpm+Tauperpmu_bar));
	
		AL0_bar=-Nbar/2./param->m_Kstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*((param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*(param->m_Bd+param->m_Kstar)*A1-lambda*A2/(param->m_Bd+param->m_Kstar))+2.*mb*(2.*E_Kstar/param->m_Bd*(param->m_Bd*param->m_Bd+3.*param->m_Kstar*param->m_Kstar-q2)*(Tauperpm+Tauperpmu_bar)-lambda/(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(Tauperpm+Tauparm+Tauperpmu_bar+Tauparmu_bar)));
		AR0_bar=-Nbar/2./param->m_Kstar/sqrt(q2)*(((C9-C9p)+(C10-C10p))*((param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*(param->m_Bd+param->m_Kstar)*A1-lambda*A2/(param->m_Bd+param->m_Kstar))+2.*mb*(2.*E_Kstar/param->m_Bd*(param->m_Bd*param->m_Bd+3.*param->m_Kstar*param->m_Kstar-q2)*(Tauperpm+Tauperpmu_bar)-lambda/(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(Tauperpm+Tauparm+Tauperpmu_bar+Tauparmu_bar)));

		At_bar=Nbar/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/param->mass_mu*(CQ2-CQ2p)/(mb_mub+param->mass_s))*E_Kstar/param->m_Kstar*xi_par;
	
		AS_bar=-2.*Nbar*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*E_Kstar/param->m_Kstar*xi_par;		
	}
	
	if(q2>6.)
	{
		if(q2>12.)
		{
		 	double tau_plus=pow(param->m_Bd+param->m_Kstar,2.);
		 	double z_q2=(sqrt(tau_plus-q2)-sqrt(tau_plus-12.))/(sqrt(tau_plus-q2)+sqrt(tau_plus-12.));
		 	
		 	double a0V=0.496;
		 	double a1V=-2.03;
		 	double c01V=1.38;
		 	double c01sV=1.066;
		 	
		 	double a0A0=0.469;
		 	double a1A0=-2.11;
		 	double c01A0=2.46;
		 	double c01sA0=0.421;
		
		 	double a0A1=0.286;
		 	double a1A1=0.19;
		 	double c01A1=1.07;
		 	double c01sA1=0.841;
		 	
		 	double a0A12=0.216;
		 	double a1A12=0.32;
		 	double c01A12=-0.12;
		 	double c01sA12=0.151;
		
		 	double PF_V=1.-q2/pow(param->m_Bd+0.135,2.);
		 	double PF_A0=1.-q2/pow(param->m_Bd+0.087,2.);
		 	double PF_A1=1.-q2/pow(param->m_Bd+0.55,2.);
		 	double PF_A12=1.-q2/pow(param->m_Bd+0.55,2.);
			
			double f_pi=0.132;
			double m_pi=0.3134;
			double m_etas=0.7319; 
			double mphys_pi=0.140;
			double mphys_etas=0.686;

			double Deltax=(m_pi*m_pi-mphys_pi*mphys_pi)/pow(4.*pi*f_pi,2.);
			double Deltaxs=(m_etas*m_etas-mphys_etas*mphys_etas)/pow(4.*pi*f_pi,2.);

			V=(a0V*(1.+c01V*Deltax+c01sV*Deltaxs)+a1V*z_q2)/PF_V;
			double A0=(a0A0*(1.+c01A0*Deltax+c01sA0*Deltaxs)+a1A0*z_q2)/PF_A0;
			A1=(a0A1*(1.+c01A1*Deltax+c01sA1*Deltaxs)+a1A1*z_q2)/PF_A1;
			double A12=(a0A12*(1.+c01A12*Deltax+c01sA12*Deltaxs)+a1A12*z_q2)/PF_A12;
			
			double lambda=(pow(param->m_Bd+param->m_Kstar,2.)-q2)*(pow(param->m_Bd-param->m_Kstar,2.)-q2);
			
			A2=(pow(param->m_Bd+param->m_Kstar,2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*A1-16.*param->m_Bd*param->m_Kstar*param->m_Kstar*(param->m_Bd+param->m_Kstar)*A12)/lambda;
					 	
		 	xi_par=param->m_Kstar/E_Kstar*A0;
		 	
		 	xi_perp=param->m_Bd/(param->m_Bd+param->m_Kstar)*V;
		}
		
		double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
	
		double shat=q2/param->m_Bd/param->m_Bd;
	
		double z=4.*mb*mb/q2;
		
		double complex x1=0.5+0.5*I*csqrt(z-1.);
		double complex x2=0.5-0.5*I*csqrt(z-1.);
		double complex x3=0.5+0.5*I/csqrt(z-1.);
		double complex x4=0.5-0.5*I/csqrt(z-1.);
		
		double complex A=
		-104./243.*log(mb*mb/mu_b/mu_b)+4.*shat/27./(1.-shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*log(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(cabs(z-1.)))+9.*shat*shat*shat*log(shat)*log(shat)+18.*pi*I*shat*(1.-2.*shat)*log(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
		
		double complex B=
		8./243./shat*((4.-34.*shat-17.*pi*I*shat)*log(mb*mb/mu_b/mu_b)+8.*shat*pow(log(mb*mb/mu_b/mu_b),2.)+17.*shat*log(shat)*log(mb*mb/mu_b/mu_b))
		+(2.+shat)*csqrt(z-1.)/729./shat*(-48.*log(mb*mb/mu_b/mu_b)*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
		-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
		-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
		-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
		+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
		-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
		+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
		double complex C=-16./81.*log(q2/mu_b/mu_b)+428./243.-64./27.*zeta3+16./81.*pi*I;
		
		double kappa=1.-2.*alphas_mub/3./pi*log(mu_b/mb);
		
		double complex C9eff=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(shat))
		+8.*mc*mc/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts))+2.*Cmub[3]+20.*Cmub[5]);
		
		double complex C9eff_bar=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(shat))
		+8.*mc*mc/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vus*conj(param->Vub)/param->Vts/conj(param->Vtb))+2.*Cmub[3]+20.*Cmub[5]);

		double complex C7eff=Cmub[7]
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(shat,log(mu_b/mb)));
		
		double lambda_hat=1.+shat*shat+pow(param->m_Kstar/param->m_Bd,4.)-2.*(shat+shat*param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd+param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd);
		
		double complex N=param->Vtb*conj(param->Vts)*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)*param->m_Bd*shat*sqrt(lambda_hat));
				
		double complex Nbar=conj(param->Vtb)*param->Vts*sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)*param->m_Bd*shat*sqrt(lambda_hat));
				
		double complex f_perp=N*param->m_Bd*sqrt(2.*lambda_hat)/(1.+param->m_Kstar/param->m_Bd)*V;
		
		double complex f_par=N*param->m_Bd*sqrt(2.)*(1.+param->m_Kstar/param->m_Bd)*A1;
		
		double complex f_0=N*param->m_Bd*((1.-shat-param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd)*pow(1.+param->m_Kstar/param->m_Bd,2.)*A1-lambda_hat*A2)/(2.*param->m_Kstar/param->m_Bd*(1.+param->m_Kstar/param->m_Bd)*sqrt(shat));
		
		double complex f_perp_bar=Nbar*param->m_Bd*sqrt(2.*lambda_hat)/(1.+param->m_Kstar/param->m_Bd)*V;
		
		double complex f_par_bar=Nbar*param->m_Bd*sqrt(2.)*(1.+param->m_Kstar/param->m_Bd)*A1;
		
		double complex f_0_bar=Nbar*param->m_Bd*((1.-shat-param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd)*pow(1.+param->m_Kstar/param->m_Bd,2.)*A1-lambda_hat*A2)/(2.*param->m_Kstar/param->m_Bd*(1.+param->m_Kstar/param->m_Bd)*sqrt(shat));
		
		ALperp_high=((C9eff-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_perp;		
		ARperp_high=((C9eff+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_perp;
		
		ALpar_high=-((C9eff-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_par;
			
		ARpar_high=-((C9eff+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_par;
		
		AL0_high=-((C9eff-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_0;
			
		AR0_high=-((C9eff+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_0;
				
				
		ALperp_bar_high=((C9eff_bar-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_perp_bar;		
		ARperp_bar_high=((C9eff_bar+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_perp_bar;
		
		ALpar_bar_high=-((C9eff_bar-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_par_bar;		
		ARpar_bar_high=-((C9eff_bar+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_par_bar;
		
		AL0_bar_high=-((C9eff_bar-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_0_bar;		
		AR0_bar_high=-((C9eff_bar+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_0_bar;
						
		double C10=Cmub[10];
		double C10p=Cpb[10];
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double lambda=pow(param->m_Bd,4.)+pow(param->m_Kstar,4.)+q2*q2-2.*(param->m_Bd*param->m_Bd*param->m_Kstar*param->m_Kstar+param->m_Kstar*param->m_Kstar*q2+param->m_Bd*param->m_Bd*q2);
		At_high=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/param->mass_mu*(CQ2-CQ2p)/(mb_mub+param->mass_s))*E_Kstar/param->m_Kstar*xi_par;
	
		AS_high=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*E_Kstar/param->m_Kstar*xi_par;
		
		At_bar_high=Nbar/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/param->mass_mu*(CQ2-CQ2p)/(mb_mub+param->mass_s))*E_Kstar/param->m_Kstar*xi_par;
	
		AS_bar_high=-2.*Nbar*sqrt(lambda)*(CQ1-CQ1p)/(mb_mub+param->mass_s)*E_Kstar/param->m_Kstar*xi_par;


		if(q2>14.)
		{
			ALperp=ALperp_high;
			ARperp=ARperp_high;
			ALpar=ALpar_high;
			ARpar=ARpar_high;
			AL0=AL0_high;
			AR0=AR0_high;
			At=At_high;
			AS=AS_high;

			ALperp_bar=ALperp_bar_high;
			ARperp_bar=ARperp_bar_high;
			ALpar_bar=ALpar_bar_high;
			ARpar_bar=ARpar_bar_high;
			AL0_bar=AL0_bar_high;
			AR0_bar=AR0_bar_high;
			At_bar=At_bar_high;
			AS_bar=AS_bar_high;
		}
		else
		{
			ALperp=ALperp*(14.-q2)/8.+ALperp_high*(q2-6.)/8.;
			ARperp=ARperp*(14.-q2)/8.+ARperp_high*(q2-6.)/8.;
			ALpar=ALpar*(14.-q2)/8.+ALpar_high*(q2-6.)/8.;
			ARpar=ARpar*(14.-q2)/8.+ARpar_high*(q2-6.)/8.;
			AL0=AL0*(14.-q2)/8.+AL0_high*(q2-6.)/8.;
			AR0=AR0*(14.-q2)/8.+AR0_high*(q2-6.)/8.;
			At=At*(14.-q2)/8.+At_high*(q2-6.)/8.;
			AS=AS*(14.-q2)/8.+AS_high*(q2-6.)/8.;

			ALperp_bar=ALperp_bar*(14.-q2)/8.+ALperp_bar_high*(q2-6.)/8.;
			ARperp_bar=ARperp_bar*(14.-q2)/8.+ARperp_bar_high*(q2-6.)/8.;
			ALpar_bar=ALpar_bar*(14.-q2)/8.+ALpar_bar_high*(q2-6.)/8.;
			ARpar_bar=ARpar_bar*(14.-q2)/8.+ARpar_bar_high*(q2-6.)/8.;
			AL0_bar=AL0_bar*(14.-q2)/8.+AL0_bar_high*(q2-6.)/8.;
			AR0_bar=AR0_bar*(14.-q2)/8.+AR0_bar_high*(q2-6.)/8.;
			At_bar=At_bar*(14.-q2)/8.+At_bar_high*(q2-6.)/8.;
			AS_bar=AS_bar*(14.-q2)/8.+AS_bar_high*(q2-6.)/8.;
		}
	}
							
	double A02=AL0*conj(AL0)+AR0*conj(AR0);
	double Apar2=ALpar*conj(ALpar)+ARpar*conj(ARpar);
	double Aperp2=ALperp*conj(ALperp)+ARperp*conj(ARperp);
	
	double A02_bar=AL0_bar*conj(AL0_bar)+AR0_bar*conj(AR0_bar);
	double Apar2_bar=ALpar_bar*conj(ALpar_bar)+ARpar_bar*conj(ARpar_bar);
	double Aperp2_bar=ALperp_bar*conj(ALperp_bar)+ARperp_bar*conj(ARperp_bar);
	
	
	double J1s=0.25*(2.+beta_l*beta_l)*(Aperp2 + Apar2) + 4.*param->mass_mu*param->mass_mu/q2*creal(ALperp*conj(ARperp)+ALpar*conj(ARpar));

	double J1c=A02 + 4.*param->mass_mu*param->mass_mu/q2*(At*conj(At)+2.*creal(AL0*conj(AR0)))+beta_l*beta_l*AS*conj(AS);

	double J2s=0.25*beta_l*beta_l*(Aperp2+Apar2);

	double J2c=-beta_l*beta_l*A02;
	
	double J3=0.5*beta_l*beta_l*(Aperp2-Apar2);
	
	double J4=1./sqrt(2.)*beta_l*beta_l*creal(AL0*conj(ALpar)+AR0*conj(ARpar));
	
	double J5=sqrt(2.)*beta_l*(creal(AL0*conj(ALperp)-AR0*conj(ARperp))-param->mass_mu/sqrt(q2)*(creal(ALpar*conj(AS)+ARpar*conj(AS))));
	
	double J6s=2.*beta_l*creal(ALpar*conj(ALperp)-ARpar*conj(ARperp));

	double J6c=4.*beta_l*param->mass_mu/sqrt(q2)*creal(AL0*conj(AS)+AR0*conj(AS));

	double J7=sqrt(2.)*beta_l*(cimag(AL0*conj(ALpar)-AR0*conj(ARpar))+param->mass_mu/sqrt(q2)*(cimag(ALperp*conj(AS)+ARperp*conj(AS))));

	double J8=1./sqrt(2.)*beta_l*beta_l*cimag(AL0*conj(ALperp)+AR0*conj(ARperp));
	
	double J9=beta_l*beta_l*cimag(ALpar*conj(ALperp)+ARpar*conj(ARperp));

	double J1s_bar=0.25*(2.+beta_l*beta_l)*(Aperp2_bar + Apar2_bar) + 4.*param->mass_mu*param->mass_mu/q2*creal(ALperp_bar*conj(ARperp_bar)+ALpar_bar*conj(ARpar_bar));

	double J1c_bar=A02_bar + 4.*param->mass_mu*param->mass_mu/q2*(At_bar*conj(At_bar)+2.*creal(AL0_bar*conj(AR0_bar)))+beta_l*beta_l*AS_bar*conj(AS_bar);

	double J2s_bar=0.25*beta_l*beta_l*(Aperp2_bar+Apar2_bar);

	double J2c_bar=-beta_l*beta_l*A02_bar;
	
	double J3_bar=0.5*beta_l*beta_l*(Aperp2_bar-Apar2_bar);
	
	double J4_bar=1./sqrt(2.)*beta_l*beta_l*creal(AL0_bar*conj(ALpar_bar)+AR0_bar*conj(ARpar_bar));
	
	double J5_bar=sqrt(2.)*beta_l*(creal(AL0_bar*conj(ALperp_bar)-AR0_bar*conj(ARperp_bar))-param->mass_mu/sqrt(q2)*(creal(ALpar_bar*conj(AS_bar)+ARpar_bar*conj(AS_bar))));
	
	double J6s_bar=2.*beta_l*creal(ALpar_bar*conj(ALperp_bar)-ARpar_bar*conj(ARperp_bar));

	double J6c_bar=4.*beta_l*param->mass_mu/sqrt(q2)*creal(AL0_bar*conj(AS_bar)+AR0_bar*conj(AS_bar));

	double J7_bar=sqrt(2.)*beta_l*(cimag(AL0_bar*conj(ALpar_bar)-AR0_bar*conj(ARpar_bar))+param->mass_mu/sqrt(q2)*(cimag(ALperp_bar*conj(AS_bar)+ARperp_bar*conj(AS_bar))));

	double J8_bar=1./sqrt(2.)*beta_l*beta_l*cimag(AL0_bar*conj(ALperp_bar)+AR0_bar*conj(ARperp_bar));
	
	double J9_bar=beta_l*beta_l*cimag(ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar));	
	
	double dGamma_BKstarmumu_dq2=3./4.*(2.*J1s+J1c-(2.*J2s+J2c)/3.);

	double dGamma_bar_BKstarmumu_dq2=3./4.*(2.*J1s_bar+J1c_bar-(2.*J2s_bar+J2c_bar)/3.);

	double AFB[3],FL[3],FT[3],AT1[3],AT2[3],AT3[3],AT4[3],AT5[3],HT1[3],HT2[3],HT3[3],alpha_Kstar[3],AIm[3],P2[3],P3[3],P6[3],P8[3],P4prime[3],P5prime[3],P6prime[3],P8prime[3];

	AFB[0]=-3./4.*(J6s+J6s_bar)/(dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2);
	AFB[1]=-3./4.*(J6s+J6s_bar);
	AFB[2]=dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2;

	FL[0]=-(J2c+J2c_bar)/(dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2);
	FL[1]=-(J2c+J2c_bar);
	FL[2]=dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2;

	FT[0]=4.*(J2s+J2s_bar)/(dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2);
	FT[1]=4.*(J2s+J2s_bar);
	FT[2]=dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2;

	AT1[0]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp)+ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar))/(Apar2+Aperp2+Apar2_bar+Aperp2_bar);
	AT1[1]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp)+ALpar_bar*conj(ALperp_bar)+ARpar_bar*conj(ARperp_bar));
	AT1[2]=Apar2+Aperp2+Apar2_bar+Aperp2_bar;

	AT2[0]=(J3+J3_bar)/2./(J2s+J2s_bar);
	AT2[1]=J3+J3_bar;
	AT2[2]=2.*(J2s+J2s_bar);

	AT3[0]=sqrt(fabs((4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar))/(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar))));
	AT3[1]=sqrt(fabs(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar)));
	AT3[2]=sqrt(fabs(-2.*(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	
	AT4[0]=sqrt((beta_l*beta_l*(J5+J5_bar)*(J5+J5_bar)+4.*(J8+J8_bar)*(J8+J8_bar))/(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar)));
	AT4[1]=sqrt(beta_l*beta_l*(J5+J5_bar)*(J5+J5_bar)+4.*(J8+J8_bar)*(J8+J8_bar));
	AT4[2]=sqrt(4.*(J4+J4_bar)*(J4+J4_bar)+beta_l*beta_l*(J7+J7_bar)*(J7+J7_bar));
	
	AT5[0]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp)+ALperp_bar*conj(ARpar_bar)+ALpar_bar*conj(ARperp_bar))/(Apar2+Aperp2+Apar2_bar+Aperp2_bar);
	AT5[1]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp)+ALperp_bar*conj(ARpar_bar)+ALpar_bar*conj(ARperp_bar));
	AT5[2]=Apar2+Aperp2+Apar2_bar+Aperp2_bar;
	
	HT1[0]=sqrt(2.)*(J4+J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s-J3+2.*J2s_bar-J3_bar)));
	HT1[1]=sqrt(2.)*(J4+J4_bar);
	HT1[2]=2.*J2s-J3+2.*J2s_bar-J3_bar;
	
	HT2[0]=beta_l*(J5+J5_bar)/sqrt(2.)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	HT2[1]=beta_l*(J5+J5_bar)/sqrt(2.);
	HT2[2]=2.*J2s+J3+2.*J2s_bar+J3_bar;
	
	HT3[0]=(J6s+J6s_bar)/2./sqrt(fabs(4.*(J2s+J2s_bar)*(J2s+J2s_bar)-(J3+J3_bar)*(J3+J3_bar)));
	HT3[1]=(J6s+J6s_bar)/2.;
	HT3[2]=J3+J3_bar;

	alpha_Kstar[0]=-(2.*J2s+J2c+2.*J2s_bar+J2c_bar)/2./(J2s+J2s_bar);
	alpha_Kstar[1]=-(2.*J2s+J2c+2.*J2s_bar+J2c_bar);
	alpha_Kstar[2]=2.*(J2s+J2s_bar);

	AIm[0]=(J9+J9_bar)/(dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2);
	AIm[1]=J9+J9_bar;
	AIm[2]=dGamma_BKstarmumu_dq2+dGamma_bar_BKstarmumu_dq2;

	P2[0]=(J6s+J6s_bar)/8./(J2s+J2s_bar);
	P2[1]=J6s+J6s_bar;
	P2[2]=8.*(J2s+J2s_bar);

	P3[0]=-(J9+J9_bar)/4./(J2s+J2s_bar);
	P3[1]=-J9+J9_bar;
	P3[2]=4.*(J2s+J2s_bar);
	
	P6[0]=-beta_l*(J7+J7_bar)/sqrt(2.)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	P6[1]=-beta_l*(J7+J7_bar)/sqrt(2.);
	P6[2]=2.*J2s+J3+2.*J2s_bar+J3_bar;
	
	P4prime[0]=(J4+J4_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P4prime[1]=J4+J4_bar;
	P4prime[2]=-(J2c+J2c_bar);
	
	P5prime[0]=(J5+J5_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P5prime[1]=(J5+J5_bar)/2.;
	P5prime[2]=-(J2c+J2c_bar);
	
	P6prime[0]=-(J7+J7_bar)/2./sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P6prime[1]=-(J7+J7_bar)/2.;
	P6prime[2]=-(J2c+J2c_bar);

	P8[0]=-sqrt(2.)*(J8+J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(2.*J2s+J3+2.*J2s_bar+J3_bar)));
	P8[1]=-sqrt(2.)*(J8+J8_bar);
	P8[2]=2.*J2s+J3+2.*J2s_bar+J3_bar;
	
	P8prime[0]=-(J8+J8_bar)/sqrt(fabs(-(J2c+J2c_bar)*(J2s+J2s_bar)));
	P8prime[1]=-(J8+J8_bar);
	P8prime[2]=-(J2c+J2c_bar);
		
	for(je=0;je<=Nobs_BKsll;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=FL[ie];
		obs[3][ie]=FT[ie];
		obs[4][ie]=AT1[ie];
		obs[5][ie]=AT2[ie]; /* =P1 */
		obs[6][ie]=AT3[ie];
		obs[7][ie]=AT4[ie];
		obs[8][ie]=AT5[ie];
		obs[9][ie]=HT1[ie]; /* =P4 */
		obs[10][ie]=HT2[ie]; /* =P5 */
		obs[11][ie]=HT3[ie];
		obs[12][ie]=alpha_Kstar[ie];
		obs[13][ie]=AIm[ie];
		obs[14][ie]=P2[ie];
		obs[15][ie]=P3[ie];
		obs[16][ie]=P6[ie];
		obs[17][ie]=P4prime[ie];
		obs[18][ie]=P5prime[ie];
		obs[19][ie]=P6prime[ie];
		obs[20][ie]=P8[ie];
		obs[21][ie]=P8prime[ie];
	}

	return dGamma_BKstarmumu_dq2;
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarmumu_dq2_calculator(double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return dGamma_BKstarmumu_dq2(q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu(double smin, double smax, double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=100.;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_BKsll+1],obs_den[Nobs_BKsll+1];
	for(je=0;je<=Nobs_BKsll;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated FL */
	obs[3]=0.; /* integrated FT */
	obs[4]=0.; /* integrated AT1 */
	obs[5]=0.; /* integrated AT2=P1 */
	obs[6]=0.; /* integrated AT3 */
	obs[7]=0.; /* integrated AT4 */
	obs[8]=0.; /* integrated AT5 */
	obs[9]=0.; /* integrated HT1=P4 */
	obs[10]=0.; /* integrated HT2=P5 */
	obs[11]=0.; /* integrated HT3 */
	obs[12]=0.; /* integrated alpha */
	obs[13]=0.; /* integrated AIm */
	obs[14]=0.; /* integrated P2 */
	obs[15]=0.; /* integrated P3 */
	obs[16]=0.; /* integrated P6 */
	obs[17]=0.; /* integrated P4prime */
	obs[18]=0.; /* integrated P5prime */
	obs[19]=0.; /* integrated P6prime */
	obs[20]=0.; /* integrated P8 */
	obs[21]=0.; /* integrated P8prime */
	
	double dobs[Nobs_BKsll+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=0.;
	s0p=1.;
	
	for(ie=1;ie<=nmax;ie++)
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s=smin+(smax-smin)*ie/nmax;
		s0p=s;	
		Gamma+=dGamma_BKstarmumu_dq2(s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		
		for(je=1;je<=Nobs_BKsll;je++) 
		{
			obs_num[je]+=dobs[je][1];
			obs_den[je]+=dobs[je][2];
		}
		
		if((dAFBtmp/dobs[1][0]<0.)&&(ie>1)) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	Gamma*=(smax-smin)/nmax;
	obs[0]=s0;
	for(je=1;je<=Nobs_BKsll;je++) 
	if((je==17)||(je==18)||(je==19)||(je==21))
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[15]/4.));
	}
	else if(je==11)
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[5]*obs_den[5]-obs_den[je]*obs_den[je]));
	}
	else if((je==9)||(je==10)||(je==16)||(je==20))
	{
		obs[je]=obs_num[je]/sqrt(fabs(obs_den[je]*obs_den[19]));
	}
	else obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_BKsll;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	return param->life_Bd/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_lowq2(double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarmumu(1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_highq2(double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarmumu(14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[1];
}

/*----------------------------------------------------------------------*/

double A_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[1];
}

/*----------------------------------------------------------------------*/

double FL_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[2];
}

/*----------------------------------------------------------------------*/

double FL_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[2];
}

/*----------------------------------------------------------------------*/

double FT_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[3];
}

/*----------------------------------------------------------------------*/

double FT_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[3];
}

/*----------------------------------------------------------------------*/

double AT1_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[4];
}

/*----------------------------------------------------------------------*/

double AT1_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[4];
}

/*----------------------------------------------------------------------*/

double AT2_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[5];
}

/*----------------------------------------------------------------------*/

double AT2_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[5];
}

/*----------------------------------------------------------------------*/

double AT3_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[6];
}

/*----------------------------------------------------------------------*/

double AT3_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[6];
}

/*----------------------------------------------------------------------*/

double AT4_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[7];
}

/*----------------------------------------------------------------------*/

double AT4_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[7];
}

/*----------------------------------------------------------------------*/

double AT5_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[8];
}

/*----------------------------------------------------------------------*/

double AT5_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[8];
}

/*----------------------------------------------------------------------*/

double HT1_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[9];
}

/*----------------------------------------------------------------------*/

double HT1_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[9];
}

/*----------------------------------------------------------------------*/

double HT2_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[10];
}

/*----------------------------------------------------------------------*/

double HT2_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[10];
}

/*----------------------------------------------------------------------*/

double HT3_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[11];
}

/*----------------------------------------------------------------------*/

double HT3_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[11];
}

/*----------------------------------------------------------------------*/

double alpha_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[12];
}

/*----------------------------------------------------------------------*/

double alpha_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[12];
}

/*----------------------------------------------------------------------*/

double AIm_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[13];
}

/*----------------------------------------------------------------------*/

double AIm_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[13];
}

/*----------------------------------------------------------------------*/

double P2_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[14];
}

/*----------------------------------------------------------------------*/

double P2_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[14];
}

/*----------------------------------------------------------------------*/

double P3_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[15];
}

/*----------------------------------------------------------------------*/

double P3_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[15];
}

/*----------------------------------------------------------------------*/

double P6_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[16];
}

/*----------------------------------------------------------------------*/

double P6_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[16];
}

/*----------------------------------------------------------------------*/

double P8_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[20];
}

/*----------------------------------------------------------------------*/

double P8_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[20];
}

/*----------------------------------------------------------------------*/

double P4prime_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[17];
}

/*----------------------------------------------------------------------*/

double P4prime_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[17];
}

/*----------------------------------------------------------------------*/

double P5prime_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[18];
}

/*----------------------------------------------------------------------*/

double P5prime_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[18];
}

/*----------------------------------------------------------------------*/

double P6prime_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[19];
}

/*----------------------------------------------------------------------*/

double P6prime_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[19];
}

/*----------------------------------------------------------------------*/

double P8prime_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[21];
}

/*----------------------------------------------------------------------*/

double P8prime_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[21];
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_lowq2_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */

	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_highq2_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BKstarmumu_zero_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	double smin=pow(2.*param.mass_mu,2.);
	double smax=pow(param.m_Bd-param.m_Kstar,2.)*0.999; 
	
	BRBKstarmumu(smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[0];
}

/*----------------------------------------------------------------------*/

double dAI_BKstarmumu_dq2(double q2, double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	double mc=mc_pole(param);

	int ie;
	double shat=q2/param->m_B/param->m_B;
	
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);

	double mu_f=sqrt(mu_b*0.5);
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
	double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi;

	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(param->m_B*param->m_B+param->m_Kstar*param->m_Kstar-q2)/2./param->m_B;

	int nf=5;
	double f_K_perp=param->f_K_perp;
	f_K_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_K_par=param->f_K_par;

	double V=0.;
	double A1=0.;
	double A2=0.;
		
	double xi_par=0.;
	double xi_perp=0.;	

	if(q2<12.)
	{
		double tau_plus=pow(param->m_Bd+param->m_Kstar,2.);
		double tau_minus=pow(param->m_Bd-param->m_Kstar,2.);
		double tau_0=tau_plus-sqrt((tau_plus-tau_minus)*tau_plus);
		double z_q2tau0=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));
		double z_0tau0=(sqrt(tau_plus)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus)+sqrt(tau_plus-tau_0));
		
		double V0=0.36;
		double b1V=-4.8;
		double A00=0.29;
		double b1A0=-18.2*0.3;
		V=V0/(1.-q2/param->m_Bs/param->m_Bs)*(1.+b1V*(z_q2tau0-z_0tau0+0.5*(z_q2tau0*z_q2tau0-z_0tau0*z_0tau0)));
		double A0=A00/(1.-q2/param->m_Bs/param->m_Bs)*(1.+b1A0*(z_q2tau0-z_0tau0+0.5*(z_q2tau0*z_q2tau0-z_0tau0*z_0tau0)));
		A1=2.*param->m_Bd*E_Kstar/pow(param->m_Bd+param->m_Kstar,2.)*V;
		A2=((param->m_Bd+param->m_Kstar)*A1/2./E_Kstar-(param->m_Kstar/E_Kstar)*A0)*param->m_Bd/(param->m_Bd-param->m_Kstar);
		
		xi_par=(param->m_Bd+param->m_Kstar)/2./E_Kstar*A1-(param->m_Bd-param->m_Kstar)/param->m_Bd*A2;
		xi_perp=param->m_Bd/(param->m_Bd+param->m_Kstar)*V;
	}
	else
	{
		double tau_plus=pow(param->m_Bd+param->m_Kstar,2.);
		double z_q2=(sqrt(tau_plus-q2)-sqrt(tau_plus-12.))/(sqrt(tau_plus-q2)+sqrt(tau_plus-12.));
		
		double a0V=0.496;
		double a1V=-2.03;
		double c01V=1.38;
		double c01sV=1.066;
		
		double a0A0=0.469;
		double a1A0=-2.11;
		double c01A0=2.46;
		double c01sA0=0.421;
		
		double a0A1=0.286;
		double a1A1=0.19;
		double c01A1=1.07;
		double c01sA1=0.841;
		
		double a0A12=0.216;
		double a1A12=0.32;
		double c01A12=-0.12;
		double c01sA12=0.151;
		
		double PF_V=1.-q2/pow(param->m_Bd+0.135,2.);
		double PF_A0=1.-q2/pow(param->m_Bd+0.087,2.);
		double PF_A1=1.-q2/pow(param->m_Bd+0.55,2.);
		double PF_A12=1.-q2/pow(param->m_Bd+0.55,2.);
			
		double f_pi=0.132;
		double m_pi=0.3134;
		double m_etas=0.7319; 
		double mphys_pi=0.140;
		double mphys_etas=0.686;

		double Deltax=(m_pi*m_pi-mphys_pi*mphys_pi)/pow(4.*pi*f_pi,2.);
		double Deltaxs=(m_etas*m_etas-mphys_etas*mphys_etas)/pow(4.*pi*f_pi,2.);

		V=(a0V*(1.+c01V*Deltax+c01sV*Deltaxs)+a1V*z_q2)/PF_V;
		double A0=(a0A0*(1.+c01A0*Deltax+c01sA0*Deltaxs)+a1A0*z_q2)/PF_A0;
		A1=(a0A1*(1.+c01A1*Deltax+c01sA1*Deltaxs)+a1A1*z_q2)/PF_A1;
		double A12=(a0A12*(1.+c01A12*Deltax+c01sA12*Deltaxs)+a1A12*z_q2)/PF_A12;
			
		double lambda=(pow(param->m_Bd+param->m_Kstar,2.)-q2)*(pow(param->m_Bd-param->m_Kstar,2.)-q2);
			
		A2=(pow(param->m_Bd+param->m_Kstar,2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*A1-16.*param->m_Bd*param->m_Kstar*param->m_Kstar*(param->m_Bd+param->m_Kstar)*A12)/lambda;
					 	
		xi_par=param->m_Kstar/E_Kstar*A0;
		 	
		xi_perp=param->m_Bd/(param->m_Bd+param->m_Kstar)*V;
	}
	
	double C7eff=Cmub[7];
	double C8eff=Cmub[8];
	double C9=Cmub[9];
	double C10=Cmub[10];
		
	double C1bar=Cmub[1]/2.;
	double C2bar=Cmub[2]-Cmub[1]/6.;
	double C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
	double C4bar=Cmub[4]/2.+8.*Cmub[6];
	double C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
	double C6bar=Cmub[4]/2.+2.*Cmub[6];
	
	double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
	+h_bkll(q2,mc,mu_b)*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
	+h_bkll(q2,param->mass_b_pole,mu_b)*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
	+h_bkll(q2,0.,mu_b)*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

	double complex C90perp=C9+Y+2.*mb*param->m_B/q2*C7eff;
	double complex C90par=C9+Y+2.*mb/param->m_B*C7eff;

	double a1perp=param->a1perp;
	double a2perp=param->a2perp;
	double a1par=param->a1par;
	double a2par=param->a2par;

	a1perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
	a2perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));

	a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
	a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

	double u;

	double lambda_Bp=param->lambda_Bp;
	lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

	double omega0=2.*(param->m_Bd-mb)/3.;
	double complex lambda_Bm=1./(exp(-q2/param->m_Bd/omega0)/omega0*(-Ei(q2/param->m_Bd/omega0)+I*pi));

	double complex integ1,integ2,integ3;
	integ1=integ2=integ3=0.;
	double complex Fperp=0.;
	double complex Xperp=0.;
	double complex Fpar=0.;
	double complex FV;
	double x;
	
	double zeta3A=param->zeta3A;
	double zeta3V=param->zeta3V;
	double wA10=param->wA10;
	double deltatp=param->deltatp;
	double deltatm=param->deltatm;
	
	double phiKstar_perp,phiKstar_par;

	int n1=10;
	int n1sav=n1;
	for(ie=0;ie<=n1;ie++)
	{
		u=(double)ie/n1;
		if(ie==0) n1*=2;
		if(ie==n1){u=(double)(ie-1)/n1;n1*=2;}
		
		phiKstar_perp=phi_Kstar(u,a1perp,a2perp);
		phiKstar_par=phi_Kstar(u,a1par,a2par);
		
		x=(1.-u)*param->m_B*param->m_B+u*q2;
		FV=3./4.*(h_bkll(x,mc,mu_b)*(C2bar+C4bar+C6bar)+h_bkll(x,param->mass_b_pole,mu_b)*(C3bar+C4bar+C6bar)+h_bkll(x,0.,mu_b)*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));
		
		integ1+=phiKstar_perp/((1.-u)+u*shat)*FV/n1;
		
		integ2+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;
	
		integ3+=phiKstar_par*FV/n1;
		
		Fperp+=phiKstar_perp/((1.-u)+u*shat)/3./n1;
		Xperp+=(u<=1.-0.5/param->m_B)*phiKstar_perp/pow((1.-u)+u*shat,2.)/3.*(0.5/param->m_B)/n1;
		Fpar+=2.*phiKstar_par/((1.-u)+u*shat)/n1;

		if(ie==0||ie==n1sav) n1=n1sav;
	}
	double rho=0.;
	double phi=0.;
	Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;
	
	double complex lambda_u=creal(conj(param->Vus)*param->Vub);
	double complex lambda_t=creal(conj(param->Vts)*param->Vtb);

	double complex K1upar_a=-lambda_u/lambda_t*(C1bar/3.+C2bar)+(C4bar+C3bar/3.);
	double complex K1dpar_a=C4bar+C3bar/3.;
	
	double complex K1perp_a=-(C6bar+C5bar/3.)*Fperp;
	
	double complex K2uperp_a=-lambda_u/lambda_t*(C1bar/3.+C2bar)+(C4bar+C3bar/3.);
	double complex K2dperp_a=C4bar+C3bar/3.;
	
	double complex K1perp_b=C8eff*mb/param->m_B*4./9.*alphas_mub/4./pi*Xperp;
	double complex K2perp_b=0.;
	double complex K1par_b=-C8eff*mb/param->m_B*4./9.*alphas_mub/4./pi*Fpar;
			
	double complex K1perp_c=4./9.*alphas_mub/4./pi*2./3.*integ1;
	double complex K2perp_c=-4./9.*alphas_mub/4./pi*1./2.*integ2;
	
	double complex K1par_c=-4./9.*alphas_mub/4./pi*2.*integ3;
	
	double complex K1perp=K1perp_a+K1perp_b+K1perp_c;
	double complex K2uperp=K2uperp_a+K2perp_c;
	double complex K2dperp=K2dperp_a+K2perp_c;
	double complex K1upar=K1upar_a+K1par_b+K1par_c;
	double complex K1dpar=K1dpar_a+K1par_b+K1par_c;
	
	double eu=2./3.;
	double ed=-1./3.;
	
	double complex bd_perp=24.*pi*pi*param->m_B*param->f_B*ed/q2/xi_perp/C90perp*(f_K_perp/param->m_B*K1perp+f_K_par*param->m_Kstar/6./lambda_Bp/param->m_B*K2dperp/(1.-q2/param->m_B/param->m_B));
	
	double complex bu_perp=24.*pi*pi*param->m_B*param->f_B*eu/q2/xi_perp/C90perp*(f_K_perp/param->m_B*K1perp+f_K_par*param->m_Kstar/6./lambda_Bp/param->m_B*K2uperp/(1.-q2/param->m_B/param->m_B));
	
	double complex bd_par=24.*pi*pi*param->f_B*ed*param->m_Kstar/param->m_B/E_Kstar/xi_par/C90par*(f_K_par/3./lambda_Bm*K1dpar);
	
	double complex bu_par=24.*pi*pi*param->f_B*eu*param->m_Kstar/param->m_B/E_Kstar/xi_par/C90par*(f_K_par/3./lambda_Bm*K1upar);	
	
	double dAI_dq2=creal(bd_perp-bu_perp)*pow(cabs(C90perp),2.)/(pow(cabs(C90perp),2.)+C10*C10)*(1.+0.25/q2*pow(E_Kstar*param->m_B/param->m_Kstar*xi_par/xi_perp*cabs(C90par/C90perp),2.)*creal(bd_par-bu_par)/creal(bd_perp-bu_perp))/(1.+0.25/q2*pow(E_Kstar*param->m_B/param->m_Kstar*xi_par/xi_perp,2.)*(pow(cabs(C90par),2.)+C10*C10)/(pow(cabs(C90perp),2.)+C10*C10));
		
	return dAI_dq2;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu(double smin, double smax, double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=100.;
	double AI=0.;
	double s;
		
	for(ie=1;ie<=nmax;ie++)
	{
		s=smin+(smax-smin)*ie/nmax;
		AI+=dAI_BKstarmumu_dq2(s,C0b,C1b,C2b,param,mu_b);
	}
	AI*=(smax-smin)/nmax;

	return AI;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_lowq2(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	return AI_BKstarmumu(1.,6.,C0b,C1b,C2b,param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_highq2(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	return AI_BKstarmumu(14.18,16.,C0b,C1b,C2b,param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_zero(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	double smin=pow(2.*param->mass_mu,2.);
	double smax=pow(param->m_B-param->m_Kstar,2.)*0.999; 

	double stemp;

	while(fabs(1.-smin/smax)>1.e-4)
	{
		stemp=(smax+smin)/2.;
		
		if(dAI_BKstarmumu_dq2(stemp,C0b,C1b,C2b,param,mu_b)>0.) smin=stemp; else smax=stemp;
	}
	return stemp;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating AI(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_lowq2(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating AI(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_highq2(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_zero_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the zero of AI(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_zero(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/
