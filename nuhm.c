#include "src/include.h"

/*#define USE_ISAJET*/
#define USE_SOFTSUSY

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to a NUHM point generated with SOFTSUSY */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	float m0,m12,tanb,A0,mu,mA,mtop,mbot,alphas_mz,delta0;

	if(argc<7) 
	{ 
		printf(" This program needs 6 parameters:\n"
		"   m0      common scalar mass at GUT scale\n"
           	"   m12     common gaugino mass at GUT scale\n"
		"   A0      trilinear soft breaking parameter at GUT scale\n"
		"   tanb      tan(beta) \n"
		"   mu      \n"
		"   mA      \n");
		printf(" Auxiliary parameters are:\n"
		"   mtop     top quark pole mass\n"
		"   mbot     Mb(Mb) scale independent b-quark mass (SOFTSUSY only)\n"
		"   alphas_mz  strong coupling at MZ (SOFTSUSY only)\n");                        
		exit(1); 
	} 
	else  
	{
		sscanf(argv[1],"%f",&m0);
		sscanf(argv[2],"%f",&m12);
		sscanf(argv[3],"%f",&A0);
		sscanf(argv[4],"%f",&tanb);
		sscanf(argv[5],"%f",&mu);
		sscanf(argv[6],"%f",&mA);
		if(argc>7) sscanf(argv[7],"%f",&mtop); else mtop=172.5;   
		if(argc>8) sscanf(argv[8],"%f",&mbot); else mbot=4.2;
		if(argc>9) sscanf(argv[9],"%f",&alphas_mz); else alphas_mz=0.1172;
	}	
	sprintf(name,"nuhm.lha");

#ifdef USE_SOFTSUSY
	softsusy_nuhm(m0, m12, tanb, A0, mu, mA, mtop, mbot, alphas_mz, name);
	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("delta0=%f\n",delta0);
       		printf("BR=%f\n",BRbsgamma_calculator(name));
		printf("excluded_mass=%d\n\n",excluded_mass_calculator(name));
	}
	else printf("Invalid point\n\n");
#endif

	sprintf(name,"rm nuhm.lha");
	system(name);		
	return 1;
}
