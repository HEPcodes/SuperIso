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
/*#define SMONLY*/
/*#define SM_ChargedHiggs*/

/*--------------------------------------------------------------------*/

#define pi    3.1415926535897932385
#define zeta3 1.2020569031595942855
#define hbar  6.58211889e-25 /* in GeV.s */

/*--------------------------------------------------------------------*/

typedef struct parameters
/* structure containing all the scanned parameters from the SLHA file */
{
	int model; /* mSUGRA = 1, GMSB = 2, AMSB = 3 */
	int generator; /* ISAJET = 1, SOFTSUSY = 2 */
	float Q; /* Qmax ; default = M_EWSB = sqrt(m_stop1*mstop2) */
	float m0,m12,tan_beta,sign_mu,A0; /* mSUGRA parameters */
	float Lambda,Mmess,N5,cgrav,m32; /* AMSB, GMSB parameters */
	float mass_Z,mass_W,mass_b,mass_top_pole,mass_tau_pole; /* SM parameters */
	float inv_alpha_em,alphas_MZ,alpha,Gfermi,GAUGE_Q; /* SM parameters */
	float charg_Umix[3][3],charg_Vmix[3][3],stop_mix[3][3],sbot_mix[3][3],stau_mix[3][3],neut_mix[6][6],mass_neut[6]; /* mass mixing matrices */
	float Min,M1_Min,M2_Min,M3_Min,At_Min,Ab_Min,Atau_Min,M2H1_Min,M2H2_Min,mu_Min,M2A_Min,tb_Min,mA_Min; /* optional input parameters at scale Min */
	float MeL_Min,MmuL_Min,MtauL_Min,MeR_Min,MmuR_Min,MtauR_Min; /* optional input parameters at scale Min */
	float MqL1_Min,MqL2_Min,MqL3_Min,MuR_Min,McR_Min,MtR_Min,MdR_Min,MsR_Min,MbR_Min; /* optional input parameters at scale Min */
	float N51,N52,N53,M2H1_Q,M2H2_Q; /* optional input parameters (N51...3: GMSB)  */
	float mass_d,mass_u,mass_s,mass_c,mass_t,mass_e,mass_nue,mass_mu,mass_num,mass_tau,mass_nut; /* SM masses */
	float mass_gluon,mass_photon,mass_Z0; /* SM masses */
	float mass_h0,mass_H0,mass_A0,mass_H,mass_dnl,mass_upl,mass_stl,mass_chl,mass_b1,mass_t1; /* Higgs & superparticle masses */
	float mass_el,mass_nuel,mass_mul,mass_numl,mass_tau1,mass_nutl,mass_gluino,mass_cha1,mass_cha2; /* superparticle masses */
	float mass_dnr,mass_upr,mass_str,mass_chr,mass_b2,mass_t2,mass_er,mass_mur,mass_tau2; /* superparticle masses */
	float mass_nuer,mass_numr,mass_nutr,mass_graviton,mass_gravitino; /* superparticle masses */
	float gp,g2,g3,YU_Q,yut[4],YD_Q,yub[4],YE_Q,yutau[4]; /* Yukawa couplings */
	float HMIX_Q,mu_Q,tanb_GUT,Higgs_VEV,mA2_Q,MSOFT_Q,M1_Q,M2_Q,M3_Q; /* parameters at scale Q */
	float MeL_Q,MmuL_Q,MtauL_Q,MeR_Q,MmuR_Q,MtauR_Q,MqL1_Q,MqL2_Q,MqL3_Q,MuR_Q,McR_Q,MtR_Q,MdR_Q,MsR_Q,MbR_Q; /* masses at scale Q */
	float AU_Q,A_u,A_c,A_t,AD_Q,A_d,A_s,A_b,AE_Q,A_e,A_mu,A_tau; /* trilinear couplings */
	
	/* SLHA2 */
	int NMSSM,Rparity,CPviolation,Flavor;
	float mass_nutau2,mass_e2,mass_nue2,mass_mu2,mass_numu2,mass_d2,mass_u2,mass_s2,mass_c2;
	float CKM_lambda,CKM_A,CKM_rho,CKM_eta;
	float PMNS_theta12,PMNS_theta23,PMNS_theta13,PMNS_delta13,PMNS_alpha1,PMNS_alpha2;
	float lambdaNMSSM_Min,kappaNMSSM_Min,AlambdaNMSSM_Min,AkappaNMSSM_Min,lambdaSNMSSM_Min,xiFNMSSM_Min,xiSNMSSM_Min,mupNMSSM_Min,mSp2NMSSM_Min,mS2NMSSM_Min,mass_H03,mass_A02,NMSSMRUN_Q,lambdaNMSSM,kappaNMSSM,AlambdaNMSSM,AkappaNMSSM,lambdaSNMSSM,xiFNMSSM,xiSNMSSM,mupNMSSM,mSp2NMSSM,mS2NMSSM; /* NMSSM parameters */
	float PMNSU_Q,CKM_Q,MSE2_Q,MSU2_Q,MSD2_Q,MSL2_Q,MSQ2_Q,TU_Q,TD_Q,TE_Q;
	
	float CKM[4][4]; /* CKM matrix */
	float H0_mix[4][4],A0_mix[4][3]; /* Higgs mixing matrices */
	float sU_mix[7][7],sD_mix[7][7],sE_mix[7][7], sNU_mix[4][4]; /* mixing matrices */
	float sCKM_msq2[4][4],sCKM_msl2[4][4],sCKM_msd2[4][4],sCKM_msu2[4][4],sCKM_mse2[4][4]; /* super CKM matrices */
	float PMNS_U[4][4]; /* PMNS mixing matrices */
	float TU[4][4],TD[4][4],TE[4][4]; /* trilinear couplings */
	
	/* non-SLHA*/
	float mass_b_1S,mass_b_pole,mtmt;
	float Lambda5; /* Lambda QCD */
	
	/* Flavor constants */
	float f_B,f_Bs,f_Ds,m_B,m_Bs,m_Ds,m_K,m_Kstar,m_D,life_B,life_Bs,life_Ds;
}
parameters;

/*--------------------------------------------------------------------*/
/* Prototypes */

/* isajet.c */
int isajet_sugra(float m0, float m12, float tanb, float A0, float sgnmu, float mtop, char name[]);

/* softsusy.c */
int softsusy_sugra(float m0, float m12, float tanb, float A0, float sgnmu, float mtop, float mbot, float alphas_mz, char name[]);
int softsusy_gmsb(float Lambda, float Mmess, float tanb, int N5, float cGrav, float sgnmu, float mtop, float mbot, float alphas_mz, char name[]);
int softsusy_amsb(float m0, float m32, float tanb, float sgnmu, float mtop, float mbot, float alphas_mz, char name[]);
int softsusy_nuhm(float m0, float m12, float tanb, float A0, float mu, float mA, float mtop, float mbot, float alphas_mz, char name[]);

/* leshouches.c */ 
int Les_Houches_Reader(char name[], struct parameters* param);
void Init_param(struct parameters* param);
void slha_adjust(struct parameters* param);
int test_slha(char name[]);

/* alphas.c */
float alphas_running(float Q, float mtop, float mbot, struct parameters* param);

/* masses.c */
float running_mass(float quark_mass, float Qinit, float Qfin,  float mtop, float mbot, struct parameters* param);
float mb_pole(struct parameters* param);
float mc_pole(struct parameters* param);
float mb_1S(struct parameters* param);
float mt_mt(struct parameters* param);

/* general.c */
float max(float x, float y);
float min(float x, float y);
float Li2(float x);
float Li3(float x);
complex float CLi2(complex float x);
float Cl2(float x);

/* wilson.c */
float A0t(float x);
float A1t(float x, float l);
float F0t(float x);
float F1t(float x,float l);
float E0t(float x);
float G1t(float x, float l);
float E1t(float x, float l);
float T(float x);
float Ech(float x);
float F7_1(float x);
float F7_2(float x);
float F7_3(float x);
float F8_1(float x);
float F8_2(float x);
float F8_3(float x);
float H2(float x, float y);
float B(float m1, float m2, float Q);
float G7H(float x, float tb);
float Delta7H(float x, float tb);
float EH(float x, float tb);
float G8H(float x, float tb);
float Delta8H(float x, float tb);
float C7t2mt(float x);
float C7c2MW(float x);
float C8t2mt(float x);
float C8c2MW(float x);
float epsilon_0(struct parameters* param);
float epsilon_2(struct parameters* param);
float epsilon_b(struct parameters* param);
float epsilon_bp(struct parameters* param);
float epsilon_0p(struct parameters* param);
float epsilon_1p(struct parameters* param);
void CW_calculator(float C0w[], float C1w[], float C2w[], float mu_W, struct parameters* param); 
void C_calculator_base1(float C0w[], float C1w[], float C2w[], float mu_W, float C0b[], float C1b[], float C2b[], float mu, struct parameters* param); 
void C_calculator_base2(float C0w[], float C1w[], float mu_W, float C0b[], float C1b[], float mu, struct parameters* param); 

/* bsgamma.c */
float phi77(float delta);
float phi78(float delta);
float phi88(float delta, float b);
float G1(float t);
float phi22(float delta, float z);
float phi11(float delta, float z);
float phi12(float delta, float z);
float G2(float t);
float phi27(float delta, float z);
float phi17(float delta, float z);
float phi18(float delta, float z);
float phi28(float delta, float z);
float phi47(float delta);
float phi48(float delta);
float F2_nf(float z);
float F2_a(float z);
float F2_na(float z);
float phi77_2beta(float delta,float mu, struct parameters* param);
float phi77_2rem(float delta, struct parameters* param);
float Re_a(float z);
float Re_b(float z);
float bsgamma(float C0[], float C1[], float C2[], float mu, float mu_W, struct parameters* param);
float bsgamma_calculator(char name[]);

/* isospin.c */
float F_orth(float a1_orth, float a2_orth);
float X_orth1(float a1_orth, float a2_orth);
float X_orth2(float a1_orth, float a2_orth);
float complex G(float s, float x);
float complex G_orth(float s, float a1_orth, float a2_orth);
float complex H_orth(float s, float a1_par, float a2_par);
float H8_orth(float a1_orth, float a2_orth);
float complex h(float u,float s);
float complex H2_orth(float s, float a1_orth, float a2_orth);
float delta0(float C0[],float C0_spec[],float C1[],float C1_spec[],struct parameters* param,float mub,float muspec, float lambda_h);
float delta0_calculator(char name[]);

/* excluded_masses.c */
int excluded_mass_calculator(char name[]);
int excluded_masses(struct parameters* param);
int charged_LSP_calculator(char name[]);
int charged_LSP(struct parameters* param);

/* gmuon.c */
float F1N(float x);
float F2N(float x);
float F1C(float x);
float F2C(float x);
float muon_gm2_calculator(char name[]);
float muon_gm2(struct parameters* param);

/* btaunu.c */
float Btaunu(struct parameters* param);
float Btaunu_calculator(char name[]);
float RBtaunu(struct parameters* param);
float RBtaunu_calculator(char name[]);

/* bdtaunu.c */
float GBDlnu(float w);
float tBDlnu(float w, float m_B, float m_D);
float rhoV(float w, float ml, float m_B, float m_D);
float rhoS(float w, float ml, float m_B, float m_D);
float dGammaBDlnu_dw(float w, float ml, struct parameters* param);
float GammaBDlnu(float ml, struct parameters* param);
float BDtaunu_BDenu(struct parameters* param);
float BDtaunu_BDenu_calculator(char name[]);
float BDtaunu(struct parameters* param);
float BDtaunu_calculator(char name[]);

/* bsmumu.c */
float D3(float x);
float D2(float x, float y);
float D1(float x, float y, float z);
void C_SUSY(struct parameters* param, float *Ccount_S, float *Ccount_P, float *Cbox_S, float *Cbox_P, float *Cpeng_S, float *Cpeng_P, float *CHp_S, float *CHp_P);
float Bsmumu(struct parameters* param);
float Bsmumu_calculator(char name[]);

/* kmunu.c */
float Kmunu_pimunu(struct parameters* param);
float Kmunu_pimunu_calculator(char name[]);
float Rl23(struct parameters* param);
float Rl23_calculator(char name[]);

/* dslnu.c */
float Dstaunu(struct parameters* param);
float Dstaunu_calculator(char name[]);
float Dsmunu(struct parameters* param);
float Dsmunu_calculator(char name[]);

