#include "src/include.h"
#include "src/isajet.h"
#include "src/softsusy.h"
#include "src/suspect.h"
#include "src/spheno.h"
#include "src/higgsbounds.h"
#include "src/hdecay.h"


#define USE_ISAJET /* to be commented if ISAJET (>=7.80) is unavailable */
#define USE_SOFTSUSY /* to be commented if SOFTSUSY is unavailable */
#define USE_SUSPECT /* to be commented if SUSPECT is unavailable */
#define USE_SPHENO /* to be commented if SPHENO is unavailable */

#define USE_HIGGSBOUNDS /* to be commented if HIGGSBOUNDS or HDECAY is unavailable */

/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to a GMSB point generated by SOFTSUSY and/or ISAJET and/or SUSPECT and/or SPHENO */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	double Lambda,Mmess,tanb,cGrav,sgnmu,mtop,mbot,alphas_mz,delta0;
	int N5;

  	if(argc<5) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   Lambda  scale of the SUSY breaking (10000-100000 GeV)\n"
           	"   Mmess   messenger mass scale > Lambda\n"
	   	"   N5      equivalent number of 5+5bar messenger fields\n"
           	"   tanb    tan(beta) \n");
    		printf(" Auxiliary parameters are:\n"
           	"   sgnmu   +/-1,  sign of Higgsino mass term (default 1)\n" 
	   	"   cGrav   >=1, ratio of the gravitino mass to its value for a breaking scale of F_m (not available for SUSPECT) \n"
           	"   mtop     top quark pole mass\n"
           	"   mbot     Mb(Mb) scale independent b-quark mass (not available for ISAJET)\n"
           	"   alphas_mz  strong coupling at MZ (not available for ISAJET)\n");                        
      		exit(1); 
  	} 
	else  
  	{  
		sscanf(argv[1],"%lf",&Lambda);
     		sscanf(argv[2],"%lf",&Mmess);
     		sscanf(argv[3],"%d",&N5);
     		sscanf(argv[4],"%lf",&tanb);
     		if(argc>5) sscanf(argv[5],"%lf",&sgnmu); else sgnmu=1.;
     		if(argc>6) sscanf(argv[6],"%lf",&cGrav); else cGrav=1.;   
     		if(argc>7) sscanf(argv[7],"%lf",&mtop); else mtop=173.3;   
     		if(argc>8) sscanf(argv[8],"%lf",&mbot); else mbot=4.19;
     		if(argc>9) sscanf(argv[9],"%lf",&alphas_mz); else alphas_mz=0.1176;
  	}

	
	int filesOK=1;
#ifdef USE_ISAJET
	if(!test_file(ISAJET)) 
	{
		printf("\"%s\" absent. Please check the ISAJET path or comment \"#define USE_ISAJET\" in gmsb.c\n",ISAJET);
		filesOK=0;
	};
#endif
#ifdef USE_SOFTSUSY
	if(!test_file(SOFTSUSY)) 
	{
		printf("\"%s\" absent. Please check the SOFTSUSY path or comment \"#define USE_SOFTSUSY\" in gmsb.c\n",SOFTSUSY);
		filesOK=0;
	};
#endif
#ifdef USE_SUSPECT
	if(!test_file(SUSPECT)) 
	{
		printf("\"%s\" absent. Please check the SUSPECT path or comment \"#define USE_SUSPECT\" in gmsb.c\n",SUSPECT);
		filesOK=0;
	};
#endif
#ifdef USE_SPHENO
	if(!test_file(SPHENO)) 
	{
		printf("\"%s\" absent. Please check the SPHENO path or comment \"#define USE_SPHENO\" in gmsb.c\n",SPHENO);
		filesOK=0;
	};
#endif
#ifdef USE_HIGGSBOUNDS
	if(!test_file(HIGGSBOUNDS)) 
	{
		printf("\"%s\" absent. Please check the HIGGSBOUNDS path or comment \"#define USE_HIGGSBOUNDS\" in gmsb.c\n",HIGGSBOUNDS);
		filesOK=0;
	};
	if(!test_file(HDECAY)) 
	{
		printf("\"%s\" absent. Please check the HDECAY path or comment \"#define USE_HIGGSBOUNDS\" in gmsb.c\n",HDECAY);
		filesOK=0;
	};
#endif
	if(!filesOK) return 1;


	if(Lambda>Mmess) printf("Lambda=%.0f must be smaller than Mmess=%.0f\n\n",Lambda,Mmess);
	else
	{

		if(!test_file("tmp")) system("mkdir tmp");
		chdir("tmp");

#ifdef USE_SOFTSUSY
		sprintf(name,"gmsb_softsusy%d.tmplha",getpid());
		softsusy_gmsb(Lambda, Mmess, tanb, N5, cGrav, sgnmu, mtop, mbot, alphas_mz, name);
		printf("SLHA file generated by SOFTSUSY\n");
		
		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{
			printf("delta0=%.3e\n",delta0);
       			printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
      			printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      			printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      			printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     			printf("Rl23=%.3e\n",Rl23_calculator(name));
      			printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      			printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
			printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
     			printf("BR_Dmunu=%.3e\n",Dmunu_calculator(name));
     			printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     			printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
			printf("a_muon=%.3e\n",muon_gm2_calculator(name));
	#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds=%d\n",higgsbounds_calculator(name));
#else
 		printf("excluded_Higgs_mass=%d\n",excluded_Higgs_mass_calculator(name));
#endif
 		printf("excluded_SUSY_mass=%d\n",excluded_SUSY_mass_calculator(name));
		flha_generator(name,"../output1.flha");
		printf("output1.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_softsusy%d.tmplha",getpid());
		system(name);		
#endif

#ifdef USE_ISAJET
		sprintf(name,"gmsb_isajet%d.tmplha",getpid());
		isajet_gmsb(Lambda, Mmess, tanb, N5, cGrav, sgnmu, mtop, name);
		printf("SLHA file generated by ISAJET\n");
		
		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{
			printf("delta0=%.3e\n",delta0);
       			printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
      			printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      			printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      			printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     			printf("Rl23=%.3e\n",Rl23_calculator(name));
      			printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      			printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
			printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
     			printf("BR_Dmunu=%.3e\n",Dmunu_calculator(name));
     			printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     			printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
			printf("a_muon=%.3e\n",muon_gm2_calculator(name));
	#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds=%d\n",higgsbounds_calculator(name));
#else
 		printf("excluded_Higgs_mass=%d\n",excluded_Higgs_mass_calculator(name));
#endif
 		printf("excluded_SUSY_mass=%d\n",excluded_SUSY_mass_calculator(name));
		flha_generator(name,"../output2.flha");
		printf("output2.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_isajet%d.tmplha",getpid());
		system(name);		
#endif
	
#ifdef USE_SUSPECT
		sprintf(name,"gmsb_suspect%d.tmplha",getpid());
		suspect_gmsb(Lambda, Mmess, tanb, N5, sgnmu, mtop, mbot, alphas_mz, name);
		printf("SLHA file generated by SUSPECT\n");
		
		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{
			printf("delta0=%.3e\n",delta0);
       			printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
      			printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      			printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      			printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     			printf("Rl23=%.3e\n",Rl23_calculator(name));
      			printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      			printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
			printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
     			printf("BR_Dmunu=%.3e\n",Dmunu_calculator(name));
     			printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     			printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
			printf("a_muon=%.3e\n",muon_gm2_calculator(name));
	#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds=%d\n",higgsbounds_calculator(name));
#else
 		printf("excluded_Higgs_mass=%d\n",excluded_Higgs_mass_calculator(name));
#endif
 		printf("excluded_SUSY_mass=%d\n",excluded_SUSY_mass_calculator(name));
		flha_generator(name,"../output3.flha");
		printf("output3.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_suspect%d.tmplha",getpid());
		system(name);		
#endif

#ifdef USE_SPHENO
		sprintf(name,"gmsb_spheno%d.tmplha",getpid());
		spheno_gmsb(Lambda, Mmess, tanb, N5, sgnmu, mtop, mbot, alphas_mz, name);
		printf("SLHA file generated by SPHENO\n");
		
		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{
			printf("delta0=%.3e\n",delta0);
       			printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
      			printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      			printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      			printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     			printf("Rl23=%.3e\n",Rl23_calculator(name));
      			printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      			printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
			printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
     			printf("BR_Dmunu=%.3e\n",Dmunu_calculator(name));
     			printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     			printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
			printf("a_muon=%.3e\n",muon_gm2_calculator(name));
#ifdef USE_HIGGSBOUNDS
			printf("excluded_HiggsBounds=%d\n",higgsbounds_calculator(name));
#else
 			printf("excluded_Higgs_mass=%d\n",excluded_Higgs_mass_calculator(name));
#endif
 			printf("excluded_SUSY_mass=%d\n",excluded_SUSY_mass_calculator(name));
		flha_generator(name,"../output4.flha");
		printf("output4.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_spheno%d.tmplha",getpid());
		system(name);		
#endif
	}
	return 1;
}
