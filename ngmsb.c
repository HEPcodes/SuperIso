#include "src/include.h"
#include "src/nmssmtools.h"
#include "src/higgsbounds.h"
#include "src/hdecay.h"

#define USE_NMSSMTOOLS /* to be commented if NMSSMTOOLS is unavailable */
#define USE_HIGGSBOUNDS /* to be commented if HIGGSBOUNDS or HDECAY is unavailable */


/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to an NGMSB point generated by NMSSMTools */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[200];
	double Lambda,Mmess,tanb,N5,lambda,AK,sgnmu,del_h,mtop,mbot,alphas_mz,delta0;

  	if(argc<6) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   Lambda  scale of the SUSY breaking (10000-100000 GeV)\n"
           	"   Mmess   messenger mass scale > Lambda\n"
	   	"   N5      equivalent number of 5+5bar messenger fields\n"
           	"   tanb    tan(beta) \n"
          	"   lambda    Yukawa coupling\n");
    		printf(" Auxiliary parameters are:\n"
          	"   sgnmu    +/-1,  sign of Higgsino mass term (default 1)\n"   
           	"   AK       (if non-zero)\n"
           	"   Del_h    (if non-zero)\n"
             	"   mtop     top quark pole mass\n"
           	"   mbot     Mb(Mb) scale independent b-quark mass\n"
           	"   alphas_mz  strong coupling at MZ\n");                        
      		exit(1); 
  	} 
	else  
  	{  
		sscanf(argv[1],"%lf",&Lambda);
     		sscanf(argv[2],"%lf",&Mmess);
     		sscanf(argv[3],"%lf",&N5);
     		sscanf(argv[4],"%lf",&tanb);
     		sscanf(argv[5],"%lf",&lambda);
     		if(argc>6) sscanf(argv[6],"%lf",&sgnmu); else sgnmu=1;
     		if(argc>7) sscanf(argv[7],"%lf",&AK); else AK=0.;
     		if(argc>8) sscanf(argv[8],"%lf",&del_h); else del_h=0.;
     		if(argc>8) sscanf(argv[8],"%lf",&mtop); else mtop=173.3;   
     		if(argc>9) sscanf(argv[9],"%lf",&mbot); else mbot=4.19;
     		if(argc>10) sscanf(argv[10],"%lf",&alphas_mz); else alphas_mz=0.1176;
  	}	

	int filesOK=1;
#ifdef USE_NMSSMTOOLS
	if(!test_file(NMSSMTools)) 
	{
		printf("\"%s\" absent. Please check the NMSSMTOOLS path or comment \"#define USE_NMSSMTOOLS\" in ngmsb.c\n",NMSSMTools);
		filesOK=0;
	};
#endif
#ifdef USE_HIGGSBOUNDS
	if(!test_file(HIGGSBOUNDS)) 
	{
		printf("\"%s\" absent. Please check the HIGGSBOUNDS path or comment \"#define USE_HIGGSBOUNDS\" in ngmsb.c\n",HIGGSBOUNDS);
		filesOK=0;
	};
	if(!test_file(HDECAY)) 
	{
		printf("\"%s\" absent. Please check the HDECAY path or comment \"#define USE_HIGGSBOUNDS\" in ngmsb.c\n",HDECAY);
		filesOK=0;
	};
#endif
	if(!filesOK) return 1;
	
	if(!test_file("tmp")) system("mkdir tmp");
	chdir("tmp");

#ifdef USE_NMSSMTOOLS	
	sprintf(name,"ngmsb_nmssmtools%d.tmplha",getpid());
	nmssmtools_ngmsb(Lambda,Mmess,tanb,(int)N5,lambda,AK,del_h,sgnmu,mtop,mbot,alphas_mz,name);
	printf("SLHA file generated by NMSSMTOOLS\n");
	
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
       		printf("charged_LSP=%d\n",charged_LSP_calculator(name));

#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds=%d\n",higgsbounds_calculator(name));
#else
 		printf("excluded_Higgs_mass=%d\n",NMSSM_collider_excluded(name));
#endif
 		printf("excluded_SUSY_mass=%d\n",excluded_SUSY_mass_calculator(name));
		printf("theory_excluded=%d\n",NMSSM_theory_excluded(name));
		flha_generator(name,"../output.flha");
		printf("output.flha generated\n\n");	
	}
	else printf("Invalid point\n\n");
	sprintf(name,"rm ngmsb_nmssmtools%d.tmplha",getpid());
	system(name);		
#endif
	return 1;
}

