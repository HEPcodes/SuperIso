#include "include.h"


double D3(double x)
{
	if(x<1.e-3) return 0.;
	if(fabs(x-1.)<1.e-3) return -1.;
	
	return x*log(x)/(1.-x);
}

/*---------------------------------------------------------------------*/

double D2(double x, double y)
{
	if(fabs(x-y)<1.e-3)
	{
		if(fabs(x-1.)<1.e-3) return -0.5;	
		return (1.-x+log(x))/(1.-x)/(1.-x);
	}

	return (D3(x)-D3(y))/(x-y);
}

/*---------------------------------------------------------------------*/

double D1(double x, double y, double z)
{
	return (D3(x)/(x-z)-D3(y)/(y-z))/(x-y) + D3(z)/(z-x)/(z-y);
}

/*---------------------------------------------------------------------*/

double Bplus(double x, double y)
{
	return y/(x-y)*(log(y)/(y-1.)-log(x)/(x-1.));
}

/*---------------------------------------------------------------------*/

void C_SUSY(struct parameters* param, double *Ccount_S, double *Ccount_P, double *Cbox_S, double *Cbox_P, double *Cpeng_S, double *Cpeng_P, double *CHp_S, double *CHp_P)
{

	if(param->SM==1)
	{
		*Ccount_S=*Ccount_P=*Cbox_S=*Cbox_P=*Cpeng_S=*Cpeng_P=*CHp_S=*CHp_P=0.;
		return;
	}
	
	double sw=sin(atan(param->gp/param->g2));
	double MU[4];
	
	MU[1]=param->mass_u;
	MU[2]=param->mass_c;
	MU[3]=param->mtmt;

 	if(param->THDM_model==0)
	{
		*CHp_S=param->mass_mu*pow(param->tan_beta/2./sw/param->mass_W,2.)*D3(param->mass_H*param->mass_H/MU[3]/MU[3])*MU[3]*MU[3]/param->mass_H/param->mass_H;
		*CHp_P=-*CHp_S;
	

#ifdef SM_ChargedHiggs
	*Ccount_S=*Ccount_P=*Cbox_S=*Cbox_P=*Cpeng_S=*Cpeng_P=0.;
	return;
#endif
	}
 	else
	{
		double CSbox=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./sw/param->mass_W,2.)*Bplus(param->mass_H*param->mass_H/param->mass_W/param->mass_W,MU[3]*MU[3]/param->mass_W/param->mass_W);

		double CPbox=-CSbox;
		
		double Pplus=-D2(param->mass_H*param->mass_H/param->mass_W/param->mass_W,MU[3]*MU[3]/param->mass_W/param->mass_W)*MU[3]*MU[3]/param->mass_W/param->mass_W;
	
		double CSpeng1=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./sw,2.)*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0);
	
		double CPpeng1=param->mass_mu*(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./sw,2.)*Pplus/param->mass_A0/param->mass_A0;
		

		double CSpeng2=param->mass_mu*(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./sw,2.)*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0*(param->mass_H*param->mass_H-param->mass_h0*param->mass_h0)/param->mass_W/param->mass_W+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0*(param->mass_H*param->mass_H-param->mass_H0*param->mass_H0)/param->mass_W/param->mass_W);
	
		double CPpeng2=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./sw,2.)*Pplus/param->mass_A0/param->mass_A0*(param->mass_H*param->mass_H-param->mass_A0*param->mass_A0)/param->mass_W/param->mass_W;
 
	
		double CSself=-param->mass_mu*param->lambda_d[3][3]*param->lambda_l[2][2]*pow(1./2./sw,2.)*(param->mass_H*param->mass_H/param->mass_W/param->mass_W+((param->lambda_d[3][3]+param->mass_s/param->mass_b*param->lambda_d[2][2])*param->lambda_u[3][3]))*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0);
	
		double CPself=param->mass_mu*param->lambda_d[3][3]*param->lambda_l[2][2]*pow(1./2./sw,2.)*(param->mass_H*param->mass_H/param->mass_W/param->mass_W+((param->lambda_d[3][3]-param->mass_s/param->mass_b*param->lambda_d[2][2])*param->lambda_u[3][3]))*Pplus/param->mass_A0/param->mass_A0;
		
	
		CSbox *= 1.+param->mass_s/param->mass_b*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_s*param->lambda_u[3][3]);
	
		CPbox *= 1.+param->mass_s/param->mass_b*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3]);
	
		CSpeng1 *= 1.+param->mass_s/param->mass_b*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3]);
	
		CPpeng1 *= 1.+param->mass_s/param->mass_b*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3]);

		CSpeng2 *= 1.+param->mass_s/param->mass_b*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3]);
	
		CPpeng2 *= 1.+param->mass_s/param->mass_b*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_b*param->lambda_u[3][3]);
	
		CSself *= 1.+param->mass_s/param->mass_b*param->lambda_d[2][2]/param->lambda_d[3][3];

		CPself *= 1.+param->mass_s/param->mass_b*param->lambda_d[2][2]/param->lambda_d[3][3];
	
	
		*CHp_S=CSbox+CSpeng1+CSpeng2+CSself;
		*CHp_P=CPbox+CPpeng1+CPpeng2+CPself;
	
		*Ccount_S=*Ccount_P=*Cbox_S=*Cbox_P=*Cpeng_S=*Cpeng_P=0.;
		return;	
	}

	double C,Ctmp,Ca,Cb,Dtmp;
	double m_cha[3], m_squ[7], m_snu[4], GammaULT[4][7], GammaURT[4][7], G_aimn[7][3][4][4], lmn[4][4];
	
	double mA=param->mass_A0;
	
	if((param->mass_A02!=0.)&&(param->mass_H03!=0.)) mA=param->mass_A02;
	
	int ie,ae,me,ne,be,ke,je;
		
	for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) GammaULT[ne][ae]=GammaURT[ne][ae]=0.;
	for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) lmn[me][ne]=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) G_aimn[ae][ie][me][ne]=0.;
		
	m_cha[1]=param->mass_cha1;
	m_cha[2]=param->mass_cha2;
	

	m_squ[1]=param->mass_upl;
	m_squ[2]=param->mass_chl;
	m_squ[3]=param->mass_t1;
	m_squ[4]=param->mass_upr;
	m_squ[5]=param->mass_chr;
	m_squ[6]=param->mass_t2;
	
	m_snu[1]=param->mass_nuel;
	m_snu[2]=param->mass_numl;
	m_snu[3]=param->mass_nutl;
	
	GammaULT[1][1]=1.;
	GammaULT[2][2]=1.;
	GammaULT[3][3]=param->stop_mix[1][1];
	GammaULT[3][6]=param->stop_mix[2][1];
	
	GammaURT[1][4]=1.;
	GammaURT[2][5]=1.;
	GammaURT[3][6]=param->stop_mix[1][1];
	GammaURT[3][3]=param->stop_mix[1][2];
	
	lmn[1][1]=-0.022;
	lmn[2][2]=-1.;
	lmn[3][3]=1.;
	
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++)
	{
		G_aimn[ae][ie][1][1]=0.5/sqrt(2.)/sw/sw*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*GammaULT[1][ae]
		-MU[1]*param->charg_Vmix[ie][2]*GammaURT[1][ae])*lmn[1][1];
	
		G_aimn[ae][ie][2][2]=0.5/sqrt(2.)/sw/sw*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*GammaULT[2][ae]
		-MU[2]*param->charg_Vmix[ie][2]*GammaURT[2][ae])*lmn[2][2];
	
		G_aimn[ae][ie][3][3]=0.5/sqrt(2.)/sw/sw*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*GammaULT[3][ae]
		-MU[3]*param->charg_Vmix[ie][2]*GammaURT[3][ae])*lmn[3][3];
	}
	

	C=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) 
	{
		Ctmp=m_cha[ie]*param->charg_Umix[ie][2]*GammaULT[me][ae]*G_aimn[ae][ie][me][ne];
		if(Ctmp!=0.) Ctmp *= D3(m_squ[ae]*m_squ[ae]/m_cha[ie]/m_cha[ie]);
		
		C += Ctmp;
	}
	
	double cfac=param->mass_mu*pow(param->tan_beta/param->mass_W/mA,2.);
	
	*Ccount_P = cfac*param->tan_beta/sqrt(2.)*C;
	*Ccount_S = - *Ccount_P;
	

	*Cpeng_S=*Cpeng_P=0.;

	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++)  for(be=1;be<=6;be++) for(ke=1;ke<=3;ke++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) 
	{
		Ctmp=G_aimn[ae][ie][me][ne]*GammaULT[me][be]*param->charg_Umix[je][2];
		
		if(Ctmp!=0.)
		{	
			Ca=Cb=0.;
			
			if((ae==be)&&(ke==me))
			{
				Dtmp=D2(m_squ[ae]*m_squ[ae]/m_cha[je]/m_cha[je],m_cha[ie]*m_cha[ie]/m_cha[je]/m_cha[je]);
				Ca+=param->mass_W*m_squ[ae]*m_squ[ae]/m_cha[je]/m_cha[je]*param->charg_Umix[je][2]*param->charg_Vmix[ie][1]*Dtmp;
				Cb+=param->mass_W*m_cha[ie]/m_cha[je]*param->charg_Umix[ie][2]*param->charg_Vmix[je][1]*Dtmp;				
			}
			if(ie==je) 
			{	
				Dtmp=D2(m_squ[ae]*m_squ[ae]/m_cha[ie]/m_cha[ie],m_squ[be]*m_squ[be]/m_cha[ie]/m_cha[ie]);
				Ca-=MU[ke]/sqrt(2.)/m_cha[ie]*param->mu_Q*GammaURT[ke][ae]*GammaULT[ke][be]*Dtmp;
				Cb-=MU[ke]/sqrt(2.)/m_cha[ie]*param->mu_Q*GammaULT[ke][ae]*GammaURT[ke][be]*Dtmp;
			}
		
			*Cpeng_S+=Ctmp*(Ca+Cb);
			*Cpeng_P+=Ctmp*(Cb-Ca);
		}
	}
		
	*Cpeng_S *= cfac;
	*Cpeng_P *= cfac;


	*Cbox_S=*Cbox_P=0.;

	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(ke=1;ke<=3;ke++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) 
	{
		Ctmp=G_aimn[ae][ie][me][ne]*GammaULT[me][ae]*param->charg_Umix[je][2]/m_cha[ie]/m_cha[ie];
		
		if(Ctmp!=0.)
		{	
			Dtmp=D1(m_snu[ke]*m_snu[ke]/m_cha[ie]/m_cha[ie],m_squ[ae]*m_squ[ae]/m_cha[ie]/m_cha[ie],m_cha[je]*m_cha[je]/m_cha[ie]/m_cha[ie]);
			Ca=m_squ[ae]*m_squ[ae]/m_cha[ie]/m_cha[ie]*param->charg_Umix[je][2]*param->charg_Vmix[ie][1]*Dtmp;
			Cb=m_cha[je]/m_cha[ie]*param->charg_Umix[ie][2]*param->charg_Vmix[je][1]*Dtmp;
		
			*Cbox_S-=Ctmp*(Ca+Cb);
			*Cbox_P+=Ctmp*(Ca-Cb);
		}
	}
		
	*Cbox_S *= param->mass_mu*pow(param->tan_beta,2.)/2./param->mass_W;
	*Cbox_P *= param->mass_mu*pow(param->tan_beta,2.)/2./param->mass_W;
		
	return;
	
}

/*---------------------------------------------------------------------*/

double Bsmumu(struct parameters* param)
/* computes the inclusive branching ratio of Bs -> mu+ mu- */
{
	double alpha_em=1./137.036;
	double VtbVts=cabs(param->Vtb*param->Vts);
	
	double Ccount_S,Ccount_P,Cbox_S,Cbox_P,Cpeng_S,Cpeng_P,CHp_S,CHp_P;
	double CA,CS,CP;
	
	C_SUSY(param,&Ccount_S,&Ccount_P,&Cbox_S,&Cbox_P,&Cpeng_S,&Cpeng_P,&CHp_S,&CHp_P);
	
	double epsfac;
	if(param->THDM_model>0) epsfac=1.;
	else epsfac=pow((1.+epsilon_b(param)*param->tan_beta),2.);

	CS=(CHp_S+Ccount_S+Cbox_S+Cpeng_S);
	CP=(CHp_P+Ccount_P+Cbox_P+Cpeng_P);
	
	CS/=epsfac;
	CP/=epsfac;
	
#ifdef DEBUG
	printf("-----------------\n");
	printf("BR(Bs -> mu+ mu-)\n");
	printf("-----------------\n");
	printf("CS=%.5e\t CP=%.5e\n",CS,CP);
#endif
	
	CA=1.033*pow(mt_mt(param)/170.,1.55)/pow(sin(atan(param->gp/param->g2)),2.);
	
	double BRmumu=param->Gfermi*param->Gfermi*alpha_em*alpha_em*pow(param->m_Bs,5.)*param->f_Bs*param->f_Bs*param->life_Bs/hbar/64./pi/pi/pi*VtbVts*VtbVts*sqrt(1.-4.*param->mass_mu*param->mass_mu/param->m_Bs/param->m_Bs)
	*( (1.-4.*param->mass_mu*param->mass_mu/param->m_Bs/param->m_Bs)*CS*CS + pow(CP-2.*CA*param->mass_mu/param->m_Bs/param->m_Bs,2.) );
	
	return BRmumu;
}

/*--------------------------------------------------------------------*/

double Bsmumu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Bs-> mu+ mu-) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Bsmumu(&param);
}
