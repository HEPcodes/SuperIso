#include "src/include.h"
#include "src/higgsbounds.h"
#include "src/hdecay.h"


/* #define USE_HIGGSBOUNDS */ /* to be commented if HIGGSBOUNDS or HDECAY is unavailable */

/*--------------------------------------------------------*/
/* Calculation of the observables using a given SLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	int test;

  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the SLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}


	int filesOK=1;
#ifdef USE_HIGGSBOUNDS
	if(!test_file(HIGGSBOUNDS)) 
	{
		printf("\"%s\" absent. Please check the HIGGSBOUNDS path or comment \"#define USE_HIGGSBOUNDS\" in slha.c\n",HIGGSBOUNDS);
		filesOK=0;
	}
	if(!test_file(HDECAY)) 
	{
		printf("\"%s\" absent. Please check the HDECAY path or comment \"#define USE_HIGGSBOUNDS\" in slha.c\n",HDECAY);
		filesOK=0;
	}
#endif
	if(!filesOK) return 1;

 
	test=test_slha(name);
	
	if(test>0)
	{
		if(test==2) printf("WARNING: only tested in the MFV scenario!\n");
		printf("delta0=%.3e\n",delta0_calculator(name));
       		printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
		printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
      		printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      		printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      		printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      		printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
     		printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     		printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
     		printf("BR_Dmunu=%.3e\n",Dmunu_calculator(name));
      		printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     		printf("Rl23=%.3e\n",Rl23_calculator(name));
		printf("a_muon=%.3e\n",muon_gm2_calculator(name));

#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds=%d\n",higgsbounds_calculator(name));
#else
 		if(test!=3) printf("excluded_Higgs_mass=%d\n",excluded_Higgs_mass_calculator(name));
#endif
 		if(test==3)
		{ 	printf("excluded_collider_NMSSMTools=%d\n",NMSSM_collider_excluded(name));
			printf("theory_excluded=%d\n",NMSSM_theory_excluded(name));
		}
		else printf("excluded_SUSY_mass=%d\n",excluded_SUSY_mass_calculator(name));	

		flha_generator(name,"output.flha");
		printf("output.flha generated\n\n");	
	}
	else if(test==-1) printf("Invalid point\n\n");
	else if(test==-2) printf("Model not yet implemented\n\n");
	else if(test==-3) printf("Invalid SLHA file\n\n");
	else if(test==-4) printf("SLHA file absent\n\n");
	
	return 1;
}
