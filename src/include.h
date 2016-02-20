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

typedef struct parameters
/* structure containing all the scanned parameters from the SLHA file */
{
	int model; /* mSUGRA = 1, GMSB = 2, AMSB = 3 */
	int generator; /* ISAJET = 1, SOFTSUSY = 2 */
	float Q; /* Qmax ; default = M_EWSB = sqrt(m_stop1*mstop2) */
	float m0,m1_2,tan_beta,sign_mu,A0,mass_W; /* mSUGRA parameters */
	float Lambda,Mmess,N5,cgrav,m3_2; /* AMSB, GMSB parameters */
	float mass_Z,mass_b,mass_top_pole,mass_tau_pole; /* SM parameters */
	float inv_alpha_em,alpha_s_MZ,alpha,Gfermi,GAUGE_Q; /* SM parameters */
	float charg_Umix[3][3],charg_Vmix[3][3],stop_mix[3][3],sbot_mix[3][3],stau_mix[3][3],neut_mix[5][5],mass_neut[5]; /* mass mixing matrices */
	float Min,M1_Min,M2_Min,M3_Min,At_Min,Ab_Min,Atau_Min,M2H1_Min,M2H2_Min,mu_Min,M2A_Min,tb_Min,mA_Min; /* optional input parameters at scale Min */
	float MeL_Min,MmuL_Min,MtauL_Min,MeR_Min,MmuR_Min,MtauR_Min; /* optional input parameters at scale Min */
	float MqL1_Min,MqL2_Min,MqL3_Min,MuR_Min,McR_Min,MtR_Min,MdR_Min,MsR_Min,MbR_Min; /* optional input parameters at scale Min */
	float N51,N52,N53,M2H1_Q,M2H2_Q; /* optional input parameters (N51...3: GMSB)  */
	float mass_d,mass_u,mass_s,mass_c,mass_t,mass_e,mass_nue,mass_mu,mass_num,mass_tau,mass_nut; /* SM masses */
	float mass_gluon,mass_photon,mass_Z0; /* SM masses */
	float mass_h0,mass_cH0,mass_A0,mass_H,mass_dnl,mass_upl,mass_stl,mass_chl,mass_b1,mass_t1; /* Higgs & superparticle masses */
	float mass_el,mass_nuel,mass_mul,mass_numl,mass_tau1,mass_nutl,mass_gluino,mass_cha1,mass_cha2; /* Higgs & superparticle masses */
	float mass_dnr,mass_upr,mass_str,mass_chr,mass_b2,mass_t2,mass_er,mass_mur,mass_tau2; /* superparticle masses */
	float mass_nuer,mass_numr,mass_nutr,mass_graviton,mass_gravitino; /* superparticle masses */
	float gp,g2,g3,YU_Q,yut,YD_Q,yub,YE_Q,yutau; /* Yukawa couplings */
	float HMIX_Q,mu_Q,tanb_GUT,Higgs_VEV,mA2_Q,MSOFT_Q,M1_Q,M2_Q,M3_Q; /* parameters at scale Q */
	float MeL_Q,MmuL_Q,MtauL_Q,MeR_Q,MmuR_Q,MtauR_Q,MqL1_Q,MqL2_Q,MqL3_Q,MuR_Q,McR_Q,MtR_Q,MdR_Q,MsR_Q,MbR_Q; /* masses at scale Q */
	float AU_Q,A_u,A_c,A_t,AD_Q,A_d,A_s,A_b,AE_Q,A_e,A_mu,A_tau; /* trilinear couplings */
}
parameters;

/*--------------------------------------------------------------------*/
/* Prototypes */
int softsusy_sugra(float m0, float m12, float tanb, float A0, float sgnmu, float mtop, float mbot, float alphas_mz, char name[]);
int isajet_sugra(float m0, float m12, float tanb, float A0, float sgnmu, float mtop, char name[]);
int softsusy_gmsb(float Lambda, float Mmess, float tanb, int N5, float cGrav, float sgnmu, float mtop, float mbot, float alphas_mz, char name[]);
int softsusy_amsb(float m0, float m32, float tanb, float sgnmu, float mtop, float mbot, float alphas_mz, char name[]);
int softsusy_nuhm(float m0, float m12, float tanb, float A0, float mu, float mA, float mtop, float mbot, float alphas_mz, char name[]);
void Init_param(struct parameters* param);
float alpha_s_running(float Q, float mtop, float mbot, float alphas_MZ, float MZ);
float running_mass(float quark_mass, float Qinit, float Qfin,  float mtop, float mbot, float alphas_MZ, float MZ);
float Li2(float x);
float E(float x);
float F4H(float x);
float F4ch(float x);
float F7_1(float x);
float F7_2(float x);
float F7_3(float x);
float F8_1(float x);
float F8_2(float x);
float F8_3(float x);
float G7(float x);
float G8(float x);
float H2(float x, float y);
float G7H(float x, float tb);
float Delta7H(float x, float tb);
float EH(float x, float tb);
float G8H(float x, float tb);
float Delta8H(float x, float tb);
float max(float x, float y);
float min(float x, float y);
float f(float x);
float f77(float x);
float f78(float x);
float f88(float x, float b);
float f22(float delta, float z);
float f27(float delta, float z);
int Les_Houches_Reader(char name[], struct parameters* param);
void C_calculator(float C0[], float C1[],float mu, struct parameters* param);
void Cem_calculator(float Cem[], float mu, struct parameters* param);
float BR_bsgamma(float C0[], float C1[], float Cem[], float mu, struct parameters* param);
float BRbsgamma_calculator(char name[]);
float delta0m(float C0[],float C0_spec[],float C1[],float C1_spec[],struct parameters* param,float mub,float muspec, float lambda_h);
float delta0_calculator(char name[]);
int excluded_mass_calculator(char name[]);
int excluded_masses(struct parameters* param);
int charged_LSP_calculator(char name[]);
int charged_LSP(struct parameters* param);
