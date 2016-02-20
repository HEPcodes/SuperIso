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

/*--------------------------------------------------------------------*/

typedef struct parameters
/* structure containing all the scanned parameters from the SLHA file */
{
	int model; /* mSUGRA = 1, GMSB = 2, AMSB = 3 */
	int generator; /* ISAJET = 1, SOFTSUSY = 2 */
	float Q; /* Qmax ; default = M_EWSB = sqrt(m_stop1*mstop2) */
	float m0,m12,tan_beta,sign_mu,A0,mass_W; /* mSUGRA parameters */
	float Lambda,Mmess,N5,cgrav,m32; /* AMSB, GMSB parameters */
	float mass_Z,mass_b,mass_top_pole,mass_tau_pole; /* SM parameters */
	float inv_alpha_em,alphas_MZ,alpha,Gfermi,GAUGE_Q; /* SM parameters */
	float charg_Umix[3][3],charg_Vmix[3][3],stop_mix[3][3],sbot_mix[3][3],stau_mix[3][3],neut_mix[6][6],mass_neut[6]; /* mass mixing matrices */
	float Min,M1_Min,M2_Min,M3_Min,At_Min,Ab_Min,Atau_Min,M2H1_Min,M2H2_Min,mu_Min,M2A_Min,tb_Min,mA_Min; /* optional input parameters at scale Min */
	float MeL_Min,MmuL_Min,MtauL_Min,MeR_Min,MmuR_Min,MtauR_Min; /* optional input parameters at scale Min */
	float MqL1_Min,MqL2_Min,MqL3_Min,MuR_Min,McR_Min,MtR_Min,MdR_Min,MsR_Min,MbR_Min; /* optional input parameters at scale Min */
	float N51,N52,N53,M2H1_Q,M2H2_Q; /* optional input parameters (N51...3: GMSB)  */
	float mass_d,mass_u,mass_s,mass_c,mass_t,mass_e,mass_nue,mass_mu,mass_num,mass_tau,mass_nut; /* SM masses */
	float mass_gluon,mass_photon,mass_Z0; /* SM masses */
	float mass_h0,mass_H0,mass_A0,mass_H,mass_dnl,mass_upl,mass_stl,mass_chl,mass_b1,mass_t1; /* Higgs & superparticle masses */
	float mass_el,mass_nuel,mass_mul,mass_numl,mass_tau1,mass_nutl,mass_gluino,mass_cha1,mass_cha2; /* Higgs & superparticle masses */
	float mass_dnr,mass_upr,mass_str,mass_chr,mass_b2,mass_t2,mass_er,mass_mur,mass_tau2; /* superparticle masses */
	float mass_nuer,mass_numr,mass_nutr,mass_graviton,mass_gravitino; /* superparticle masses */
	float gp,g2,g3,YU_Q,yut,YD_Q,yub,YE_Q,yutau; /* Yukawa couplings */
	float HMIX_Q,mu_Q,tanb_GUT,Higgs_VEV,mA2_Q,MSOFT_Q,M1_Q,M2_Q,M3_Q; /* parameters at scale Q */
	float MeL_Q,MmuL_Q,MtauL_Q,MeR_Q,MmuR_Q,MtauR_Q,MqL1_Q,MqL2_Q,MqL3_Q,MuR_Q,McR_Q,MtR_Q,MdR_Q,MsR_Q,MbR_Q; /* masses at scale Q */
	float AU_Q,A_u,A_c,A_t,AD_Q,A_d,A_s,A_b,AE_Q,A_e,A_mu,A_tau; /* trilinear couplings */
	
	/* SLHA2 */
	int NMSSM,Rparity,CPviolation,Flavor;
	float mass_nutau2,mass_e2,mass_nue2,mass_mu2,mass_numu2,mass_d2,mass_u2,mass_s2,mass_c2,CKM_lambda,CKM_A,CKM_rho,CKM_eta,PMNS_theta12,PMNS_theta23,PMNS_theta13,PMNS_delta13,PMNS_alpha1,PMNS_alpha2,lambdaNMSSM_Min,kappaNMSSM_Min,AlambdaNMSSM_Min,AkappaNMSSM_Min,lambdaSNMSSM_Min,xiFNMSSM_Min,xiSNMSSM_Min,mupNMSSM_Min,mSp2NMSSM_Min,mS2NMSSM_Min,mass_H03,mass_A02,NMSSMRUN_Q,lambdaNMSSM,kappaNMSSM,AlambdaNMSSM,AkappaNMSSM,lambdaSNMSSM,xiFNMSSM,xiSNMSSM,mupNMSSM,mSp2NMSSM,mS2NMSSM;
	float H0_mix[4][4],A0_mix[4][3];
	
	/* non-SLHA*/
	float mass_b_1S,mass_b_pole;
	float Lambda5;
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
float epsilon_b(struct parameters* param);
float epsilon_bp(struct parameters* param);
float epsilon_tp(struct parameters* param);
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
float BRbsgamma_calculator(char name[]);

/* isospin.c */
float F_orth(float a1_orth, float a2_orth);
float X_orth1(float a1_orth, float a2_orth);
float X_orth2(float a1_orth, float a2_orth);
float complex G(float s, float x);
float complex G_orth(float s, float a1_orth, float a2_orth);
float complex H_orth(float s, float a1_par, float a2_par);
float H8_orth(float a1_orth, float a2_orth);
float complex h(float u,float s);
float complex H1_orth(float s, float a1_orth, float a2_orth);
float delta0m(float C0[],float C0_spec[],float C1[],float C1_spec[],struct parameters* param,float mub,float muspec, float lambda_h);
float delta0_calculator(char name[]);

/* excluded_masses.c */
int excluded_mass_calculator(char name[]);
int excluded_masses(struct parameters* param);
int charged_LSP_calculator(char name[]);
int charged_LSP(struct parameters* param);

/* muon.c */
float F1N(float x);
float F2N(float x);
float F1C(float x);
float F2C(float x);
float muon_gm2_calculator(char name[]);
float muon_gm2(struct parameters* param);

/* btaunu.c */
float btaunu(struct parameters* param);
float btaunu_calculator(char name[]);

/* bdtaunu.c */
float GBDlnu(float w);
float tBDlnu(float w);
float rhoV(float w, float ml);
float rhoS(float w, float ml);
float dGammaBDlnu_dw(float w, float ml, struct parameters* param);
float GammaBDlnu(float ml, struct parameters* param);
float Bbdtaunu_Bbdenu(struct parameters* param);
float Bbdtaunu_Bbdenu_calculator(char name[]);
float Bbdtaunu(struct parameters* param);
float Bbdtaunu_calculator(char name[]);

/* kmunu.c */
float Bkmunu_Bpimunu(struct parameters* param);
float Bkmunu_Bpimunu_calculator(char name[]);
float Rl23(struct parameters* param);
float Rl23_calculator(char name[]);

