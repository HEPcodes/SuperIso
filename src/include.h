#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/wait.h>
#include <complex.h>
#include <string.h>
#include <strings.h>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/
/*#define SM_ChargedHiggs*/

/*--------------------------------------------------------------------*/

#define pi    3.1415926535897932385
#define zeta3 1.2020569031595942855
#define hbar  6.58211889e-25 /* in GeV.s */

/*--------------------------------------------------------------------*/

typedef struct parameters
/* structure containing all the scanned parameters from the SLHA file */
{
	int SM;
	int model; /* mSUGRA = 1, GMSB = 2, AMSB = 3 */
	int generator; /* ISAJET = 1, SOFTSUSY = 2 */
	double Q; /* Qmax ; default = M_EWSB = sqrt(m_stop1*mstop2) */
	double m0,m12,tan_beta,sign_mu,A0; /* mSUGRA parameters */
	double Lambda,Mmess,N5,cgrav,m32; /* AMSB, GMSB parameters */
	double mass_Z,mass_W,mass_b,mass_top_pole,mass_tau_pole; /* SM parameters */
	double inv_alpha_em,alphas_MZ,Gfermi,GAUGE_Q; /* SM parameters */
	double charg_Umix[3][3],charg_Vmix[3][3],stop_mix[3][3],sbot_mix[3][3],stau_mix[3][3],neut_mix[6][6],mass_neut[6],alpha; /* mass mixing matrices */
	double Min,M1_Min,M2_Min,M3_Min,At_Min,Ab_Min,Atau_Min,M2H1_Min,M2H2_Min,mu_Min,M2A_Min,tb_Min,mA_Min; /* optional input parameters at scale Min */
	double MeL_Min,MmuL_Min,MtauL_Min,MeR_Min,MmuR_Min,MtauR_Min; /* optional input parameters at scale Min */
	double MqL1_Min,MqL2_Min,MqL3_Min,MuR_Min,McR_Min,MtR_Min,MdR_Min,MsR_Min,MbR_Min; /* optional input parameters at scale Min */
	double N51,N52,N53,M2H1_Q,M2H2_Q; /* optional input parameters (N51...3: GMSB) */
	double mass_d,mass_u,mass_s,mass_c,mass_t,mass_e,mass_nue,mass_mu,mass_num,mass_tau,mass_nut; /* SM masses */
	double mass_gluon,mass_photon,mass_Z0; /* SM masses */
	double mass_h0,mass_H0,mass_A0,mass_H,mass_dnl,mass_upl,mass_stl,mass_chl,mass_b1,mass_t1; /* Higgs & superparticle masses */
	double mass_el,mass_nuel,mass_mul,mass_numl,mass_tau1,mass_nutl,mass_gluino,mass_cha1,mass_cha2; /* superparticle masses */
	double mass_dnr,mass_upr,mass_str,mass_chr,mass_b2,mass_t2,mass_er,mass_mur,mass_tau2; /* superparticle masses */
	double mass_nuer,mass_numr,mass_nutr,mass_graviton,mass_gravitino; /* superparticle masses */
	double gp,g2,g3,YU_Q,yut[4],YD_Q,yub[4],YE_Q,yutau[4]; /* Yukawa couplings */
	double HMIX_Q,mu_Q,tanb_GUT,Higgs_VEV,mA2_Q,MSOFT_Q,M1_Q,M2_Q,M3_Q; /* parameters at scale Q */
	double MeL_Q,MmuL_Q,MtauL_Q,MeR_Q,MmuR_Q,MtauR_Q,MqL1_Q,MqL2_Q,MqL3_Q,MuR_Q,McR_Q,MtR_Q,MdR_Q,MsR_Q,MbR_Q; /* masses at scale Q */
	double AU_Q,A_u,A_c,A_t,AD_Q,A_d,A_s,A_b,AE_Q,A_e,A_mu,A_tau; /* trilinear couplings */
	
	/* SLHA2 */
	int NMSSM,Rparity,CPviolation,Flavor;
	double mass_nutau2,mass_e2,mass_nue2,mass_mu2,mass_numu2,mass_d2,mass_u2,mass_s2,mass_c2;
	double CKM_lambda,CKM_A,CKM_rho,CKM_eta;
	double PMNS_theta12,PMNS_theta23,PMNS_theta13,PMNS_delta13,PMNS_alpha1,PMNS_alpha2;
	double lambdaNMSSM_Min,kappaNMSSM_Min,AlambdaNMSSM_Min,AkappaNMSSM_Min,lambdaSNMSSM_Min,xiFNMSSM_Min,xiSNMSSM_Min,mupNMSSM_Min,mSp2NMSSM_Min,mS2NMSSM_Min,mass_H03,mass_A02,NMSSMRUN_Q,lambdaNMSSM,kappaNMSSM,AlambdaNMSSM,AkappaNMSSM,lambdaSNMSSM,xiFNMSSM,xiSNMSSM,mupNMSSM,mSp2NMSSM,mS2NMSSM; /* NMSSM parameters */
	double PMNSU_Q,CKM_Q,IMCKM_Q,MSE2_Q,MSU2_Q,MSD2_Q,MSL2_Q,MSQ2_Q,TU_Q,TD_Q,TE_Q;
	
	double CKM[4][4],IMCKM[4][4]; /* CKM matrix */
	double H0_mix[4][4],A0_mix[4][4]; /* Higgs mixing matrices */
	double sU_mix[7][7],sD_mix[7][7],sE_mix[7][7], sNU_mix[4][4]; /* mixing matrices */
	double sCKM_msq2[4][4],sCKM_msl2[4][4],sCKM_msd2[4][4],sCKM_msu2[4][4],sCKM_mse2[4][4]; /* super CKM matrices */
	double PMNS_U[4][4]; /* PMNS mixing matrices */
	double TU[4][4],TD[4][4],TE[4][4]; /* trilinear couplings */
	
	/* non-SLHA*/
	double mass_b_1S,mass_b_pole,mtmt;
	double Lambda5; /* Lambda QCD */
	
	/* Flavor constants */
	double f_B,f_Bs,f_Ds,f_D,fK_fpi;
	double m_B,m_Bs,m_pi,m_Ds,m_K,m_Kstar,m_D0,m_D;
	double life_pi,life_K,life_B,life_Bs,life_D,life_Ds;
	
	/* Decay widths */
	double width_h0,width_H0,width_A0,width_H;
	
	/* CKM matrix */
	complex double Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb;

	/* 2HDM */
	int THDM_model;
	double lambda_u[4][4],lambda_d[4][4],lambda_l[4][4];
}
parameters;

/*--------------------------------------------------------------------*/
/* Prototypes */

/* isajet.c */
int isajet_sugra(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, char name[]);
int isajet_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, char name[]);
int isajet_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, char name[]);
int isajet_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, char name[]);

/* softsusy.c */
int softsusy_sugra(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, double mbot, double alphas_mz, char name[]);

/* leshouches.c */ 
int Les_Houches_Reader(char name[], struct parameters* param);
void Init_param(struct parameters* param);
void slha_adjust(struct parameters* param);
int test_slha(char name[]);

/* alphas.c */
double alphas_running(double Q, double mtop, double mbot, struct parameters* param);

/* masses.c */
double running_mass(double quark_mass, double Qinit, double Qfin,  double mtop, double mbot, struct parameters* param);
double mb_pole(struct parameters* param);
double mc_pole(struct parameters* param);
double mb_1S(struct parameters* param);
double mt_mt(struct parameters* param);

/* general.c */
double max(double x, double y);
double min(double x, double y);
double Li2(double x);
double Li3(double x);
complex double CLi2(complex double x);
double Cl2(double x);
double I0(double x);
double I1(double x);
double K0(double x);
double K1(double x);
double K2(double x);
double K0exp(double x,double z);
double K1exp(double x,double z);
double K2exp(double x,double z);
int test_integer(char name[]);

/* wilson.c */
double A0t(double x);
double A1t(double x, double l);
double F0t(double x);
double F1t(double x,double l);
double E0t(double x);
double G1t(double x, double l);
double E1t(double x, double l);
double T(double x);
double Ech(double x);
double F7_1(double x);
double F7_2(double x);
double F7_3(double x);
double F8_1(double x);
double F8_2(double x);
double F8_3(double x);
double H2(double x, double y);
double B(double m1, double m2, double Q);
double G7H(double x, double lu, double ld);
double Delta7H(double x, double lu, double ld);
double EH(double x, double lu);
double G8H(double x, double lu, double ld);
double Delta8H(double x, double lu, double ld);
double C7t2mt(double x);
double C7c2MW(double x);
double C8t2mt(double x);
double C8c2MW(double x);
double epsilon_0(struct parameters* param);
double epsilon_2(struct parameters* param);
double epsilon_b(struct parameters* param);
double epsilon_bp(struct parameters* param);
double epsilon_0p(struct parameters* param);
double epsilon_1p(struct parameters* param);
void CW_calculator(double C0w[], double C1w[], double C2w[], double mu_W, struct parameters* param); 
void C_calculator_base1(double C0w[], double C1w[], double C2w[], double mu_W, double C0b[], double C1b[], double C2b[], double mu, struct parameters* param); 
void C_calculator_base2(double C0w[], double C1w[], double mu_W, double C0b[], double C1b[], double mu, struct parameters* param); 

/* bsgamma.c */
double phi77(double delta);
double phi78(double delta);
double phi88(double delta, double b);
double G1(double t);
double phi22(double delta, double z);
double phi11(double delta, double z);
double phi12(double delta, double z);
double G2(double t);
double phi27(double delta, double z);
double phi17(double delta, double z);
double phi18(double delta, double z);
double phi28(double delta, double z);
double phi47(double delta);
double phi48(double delta);
double F2_nf(double z);
double F2_a(double z);
double F2_na(double z);
double phi77_2beta(double delta,double mu, struct parameters* param);
double phi77_2rem(double delta, struct parameters* param);
double Re_a(double z);
double Re_b(double z);
double bsgamma(double C0[], double C1[], double C2[], double mu, double mu_W, struct parameters* param);
double bsgamma_calculator(char name[]);

/* isospin.c */
double F_orth(double a1_orth, double a2_orth);
double X_orth1(double a1_orth, double a2_orth);
double X_orth2(double a1_orth, double a2_orth);
double complex G(double s, double x);
double complex G_orth(double s, double a1_orth, double a2_orth);
double complex H_orth(double s, double a1_par, double a2_par);
double H8_orth(double a1_orth, double a2_orth);
double complex h(double u,double s);
double complex H2_orth(double s, double a1_orth, double a2_orth);
double delta0(double C0[],double C0_spec[],double C1[],double C1_spec[],struct parameters* param,double mub,double muspec, double lambda_h);
double delta0_calculator(char name[]);

/* excluded_masses.c */
int excluded_mass_calculator(char name[]);
int excluded_masses(struct parameters* param);
int charged_LSP_calculator(char name[]);
int charged_LSP(struct parameters* param);

/* gmuon.c */
double F1N(double x);
double F2N(double x);
double F1C(double x);
double F2C(double x);
double fPS(double x);
double fS(double x);
double fft(double x);
double muonf(double z);
double muong(double z);
double muon_gm2_calculator(char name[]);
double muon_gm2(struct parameters* param);

/* btaunu.c */
double Btaunu(struct parameters* param);
double Btaunu_calculator(char name[]);
double RBtaunu(struct parameters* param);
double RBtaunu_calculator(char name[]);

/* bdtaunu.c */
double GBDlnu(double w);
double tBDlnu(double w, double m_B, double m_D);
double rhoV(double w, double ml, double m_B, double m_D);
double rhoS(double w, double ml, double m_B, double m_D);
double dGammaBDlnu_dw(double w, double ml, struct parameters* param);
double GammaBDlnu(double ml, struct parameters* param);
double BDtaunu_BDenu(struct parameters* param);
double BDtaunu_BDenu_calculator(char name[]);
double BDtaunu(struct parameters* param);
double BDtaunu_calculator(char name[]);

/* bsmumu.c */
double D3(double x);
double D2(double x, double y);
double D1(double x, double y, double z);
void C_SUSY(struct parameters* param, double *Ccount_S, double *Ccount_P, double *Cbox_S, double *Cbox_P, double *Cpeng_S, double *Cpeng_P, double *CHp_S, double *CHp_P);
double Bsmumu(struct parameters* param);
double Bsmumu_calculator(char name[]);

/* kmunu.c */
double Kmunu_pimunu(struct parameters* param);
double Kmunu_pimunu_calculator(char name[]);
double Rl23(struct parameters* param);
double Rl23_calculator(char name[]);

/* dslnu.c */
double Dstaunu(struct parameters* param);
double Dstaunu_calculator(char name[]);
double Dsmunu(struct parameters* param);
double Dsmunu_calculator(char name[]);

/* dmunu.c */
double Dmunu(struct parameters* param);
double Dmunu_calculator(char name[]);

/* 2hdmc.c */
int thdmc_types(double l1, double l2, double l3, double l4, double l5, double l6, double l7, double m12_2, double tanb, int type, char name[]);

