#include "src/include.h"


#define USE_ISAJET /* to be commented if ISAJET is unavailable */
#define USE_SOFTSUSY /* to be commented if SOFTSUSY is unavailable */

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to a mSUGRA point generated with SOFTSUSY and/or ISAJET */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	float m0,m12,tanb,A0,sgnmu,mtop,mbot,alphas_mz,delta0;

  	if(argc<5) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   m0      universal scalar mass at GUT scale\n"
           	"   m12     universal gaugino mass at GUT scale\n"
           	"   A0      trilinear soft breaking parameter at GUT scale\n"
           	"   tanb      tan(beta) \n");
    		printf(" Auxiliary parameters are:\n"
           	"   sgnmu    +/-1,  sign of Higgsino mass term (default 1)\n"    
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
     		if(argc>5) sscanf(argv[5],"%f",&sgnmu); else sgnmu=1;
     		if(argc>6) sscanf(argv[6],"%f",&mtop); else mtop=172.5;   
     		if(argc>7) sscanf(argv[7],"%f",&mbot); else mbot=4.2;
     		if(argc>8) sscanf(argv[8],"%f",&alphas_mz); else alphas_mz=0.1172;
  	}	
	sprintf(name,"msugra.lha");

#ifdef USE_SOFTSUSY
	softsusy_sugra(m0, m12, tanb, A0, sgnmu, mtop, mbot, alphas_mz, name);

	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("delta0_softsusy=%f\n",delta0_calculator(name));
       		printf("BR_softsusy=%f\n",BRbsgamma_calculator(name));
       		printf("charged_LSP_softsusy=%d\n",charged_LSP_calculator(name));
       		printf("excluded_masses_softsusy=%d\n\n",excluded_mass_calculator(name));
	}
	else printf("Invalid point\n\n");
#endif

#ifdef USE_ISAJET
	isajet_sugra(m0, m12, tanb, A0, sgnmu, mtop, name);

	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("delta0_isajet=%f\n",delta0);
       		printf("BR_isajet=%f\n",BRbsgamma_calculator(name));
       		printf("charged_LSP_isajet=%d\n",charged_LSP_calculator(name));
		printf("excluded_masses_isajet=%d\n\n",excluded_mass_calculator(name));
	}
	else printf("Invalid point\n\n");
#endif	

	sprintf(name,"rm msugra.lha");
	system(name);		
	return 1;
}
