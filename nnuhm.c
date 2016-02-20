#include "src/include.h"
#include "src/nmssmtools.h"

#define USE_NMSSMTOOLS /* to be commented if NMSSMTOOLS is unavailable */

/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to an NNUHM point generated by NMSSMTools */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[200];
	double m0,m12,tanb,A0,lambda,AK,sgnmu,MHDGUT,MHUGUT,mtop,mbot,alphas_mz,delta0;

  	if(argc<8) 
  	{ 
    		printf(" This program needs 7 parameters:\n"
           	"   m0      universal scalar mass at GUT scale\n"
           	"   m12     universal gaugino mass at GUT scale\n"
           	"   A0      trilinear soft breaking parameter at GUT scale\n"
           	"   tanb      tan(beta) \n"
          	"   lambda    Yukawa coupling\n"
		"   MHDGUT    Down Higgs mass at GUT scale\n"
		"   MHUGUT    Up Higgs mass at GUT scale\n");
    		printf(" Auxiliary parameters are:\n"
         	"   sgnmu    +/-1,  sign of Higgsino mass term (default 1)\n"               	"   AK  (if different from A0)\n"
            	"   mtop     top quark pole mass\n"
           	"   mbot     Mb(Mb) scale independent b-quark mass\n"
           	"   alphas_mz  strong coupling at MZ\n");                        
      		exit(1); 
  	} 
	else  
  	{  
		sscanf(argv[1],"%lf",&m0);
     		sscanf(argv[2],"%lf",&m12);
     		sscanf(argv[3],"%lf",&A0);
     		sscanf(argv[4],"%lf",&tanb);
     		sscanf(argv[5],"%lf",&lambda);
     		sscanf(argv[6],"%lf",&MHDGUT);
     		sscanf(argv[7],"%lf",&MHUGUT);	
     		if(argc>8) sscanf(argv[8],"%lf",&sgnmu); else sgnmu=1;
     		if(argc>9) sscanf(argv[9],"%lf",&AK); else AK=A0;
     		if(argc>10) sscanf(argv[10],"%lf",&mtop); else mtop=173.3;   
     		if(argc>11) sscanf(argv[11],"%lf",&mbot); else mbot=4.19;
     		if(argc>12) sscanf(argv[12],"%lf",&alphas_mz); else alphas_mz=0.1176;
  	}	

	int filesOK=1;
#ifdef USE_NMSSMTOOLS
	if(!test_file(NMSSMTools)) 
	{
		printf("\"%s\" absent. Please check the NMSSMTOOLS path or comment \"#define USE_NMSSMTOOLS\" in nnuhm.c\n",NMSSMTools);
		filesOK=0;
	}
#endif
	if(!filesOK) return 1;
	
	if(!test_file("tmp")) system("mkdir tmp");
	chdir("tmp");
	
#ifdef USE_NMSSMTOOLS	
	sprintf(name,"nnuhm_nmssmtools%d.tmplha",getpid());
	nmssmtools_nnuhm(m0, m12, tanb, A0, MHDGUT, MHUGUT, lambda, AK, sgnmu, mtop, mbot, alphas_mz,name);
	printf("SLHA file generated by NMSSMTOOLS\n");

	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("delta0=%.3e\n",delta0);
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
       		printf("charged_LSP=%d\n",charged_LSP_calculator(name));

 		printf("excluded_collider_NMSSMTools=%d\n",NMSSM_collider_excluded(name));
		printf("theory_excluded=%d\n",NMSSM_theory_excluded(name));
		flha_generator(name,"../output.flha");
		printf("output.flha generated\n\n");	
	}
	else printf("Invalid point\n\n");
	sprintf(name,"rm nnuhm_nmssmtools%d.tmplha",getpid());
	system(name);		
#endif
	return 1;
}
