#include "src/include.h"
#include "src/isajet.h"
#include "src/softsusy.h"
#include "src/higgsbounds.h"


#define USE_SOFTSUSY /* to be commented if SOFTSUSY is unavailable */
//#define USE_ISAJET /* to be commented if ISAJET (>=7.80) is unavailable */

//#define USE_HIGGSBOUNDS /* to be commented if HIGGSBOUNDS is unavailable */

/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to a NUHM point generated by SOFTSUSY and/or ISAJET */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	double m0,m12,tanb,A0,mu,mA,mtop,mbot,alphas_mz,delta0;
	double obs[Nobs_BKsll+1];

	if(argc<7) 
	{ 
		printf(" This program needs 6 parameters:\n"
		"   m0      common scalar mass at GUT scale\n"
           	"   m12     common gaugino mass at GUT scale\n"
		"   A0      trilinear soft breaking parameter at GUT scale\n"
		"   tanb    tan(beta) \n"
		"   mu      mu parameter\n"
		"   mA      CP-odd Higgs mass\n");
		printf(" Auxiliary parameters are:\n"
		"   mtop    top quark pole mass\n"
		"   mbot    Mb(Mb) scale independent b-quark mass (SOFTSUSY only)\n"
		"   alphas_mz  strong coupling at MZ (SOFTSUSY only)\n");                        
		exit(1); 
	} 
	else  
	{
		sscanf(argv[1],"%lf",&m0);
		sscanf(argv[2],"%lf",&m12);
		sscanf(argv[3],"%lf",&A0);
		sscanf(argv[4],"%lf",&tanb);
		sscanf(argv[5],"%lf",&mu);
		sscanf(argv[6],"%lf",&mA);
		if(argc>7) sscanf(argv[7],"%lf",&mtop); else mtop=173.34;   
		if(argc>8) sscanf(argv[8],"%lf",&mbot); else mbot=4.18;
		if(argc>9) sscanf(argv[9],"%lf",&alphas_mz); else alphas_mz=0.1184;
	}	

	int filesOK=1;
#ifdef USE_ISAJET
	if(!test_file(ISAJET)) 
	{
		printf("\"%s\" absent. Please check the ISAJET path or comment \"#define USE_ISAJET\" in nuhm.c\n",ISAJET);
		filesOK=0;
	}
#endif
#ifdef USE_SOFTSUSY
	if(!test_file(SOFTSUSY)) 
	{
		printf("\"%s\" absent. Please check the SOFTSUSY path or comment \"#define USE_SOFTSUSY\" in nuhm.c\n",SOFTSUSY);
		filesOK=0;
	}
#endif
#ifdef USE_HIGGSBOUNDS
	if(!test_file(HBwithFH)) 
	{
		printf("\"%s\" absent. Please check the HBwithFH path or comment \"#define USE_HIGGSBOUNDS\" in nuhm.c\n",HBwithFH);
		filesOK=0;
	}
#endif
	if(!filesOK) return 1;


	if(!test_file("tmp")) system("mkdir tmp");
	chdir("tmp");

	printf("\n");
	
	printf("SuperIso v3.4 - F. Mahmoudi\n\n");

#ifdef USE_SOFTSUSY
	sprintf(name,"nuhm_softsusy%d.tmplha",getpid());
	softsusy_nuhm(m0, m12, tanb, A0, mu, mA, mtop, mbot, alphas_mz, name);
	printf("---------------------------------------\n\n");
	printf("NUHM - SLHA file generated by SOFTSUSY\n\n");

	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("Observable\t\t\tValue\n\n");

		printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
		printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

		printf("BR(Bs->mu mu)\t\t\t%.3e\n",Bsmumu_calculator(name));
		printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag_calculator(name));
		printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu_calculator(name));
	
		printf("BR(B->K* mu mu)_low\t\t%.3e\n",BRobs_BKstarmumu_lowq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_low\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_low\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_low\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_low\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_low\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_low\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_low\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_low\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_low\t\t%.3e\n\n",AI_BKstarmumu_lowq2_calculator(name));
	
		printf("BR(B->K* mu mu)_high\t\t%.3e\n",BRobs_BKstarmumu_highq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_high\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_high\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_high\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_high\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_high\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_high\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_high\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_high\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_high\t\t%.3e\n\n",AI_BKstarmumu_highq2_calculator(name));

		printf("q0^2(AFB(B->K* mu mu))\t\t%.3e\n",A_BKstarmumu_zero_calculator(name));
		printf("q0^2(AI(B->K* mu mu))\t\t%.3e\n\n",AI_BKstarmumu_zero_calculator(name));


		printf("BR(B->Xs mu mu)_low\t\t%.3e\n",BRBXsmumu_lowq2_calculator(name));
		printf("BR(B->Xs mu mu)_high\t\t%.3e\n",BRBXsmumu_highq2_calculator(name));
		printf("q0^2(AFB(B->Xs mu mu)\t\t%.3e\n",A_BXsmumu_zero_calculator(name));
		printf("BR(B->Xs tau tau)_high\t\t%.3e\n\n",BRBXstautau_highq2_calculator(name));
	
		printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu_calculator(name));
      		printf("R(B->tau nu)\t\t\t%.3e\n",RBtaunu_calculator(name));
      		printf("BR(B->D tau nu)\t\t\t%.3e\n",BDtaunu_calculator(name));
      		printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n",BDtaunu_BDenu_calculator(name));
     		printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu_calculator(name));
     		printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu_calculator(name));
     		printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu_calculator(name));
      		printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu_calculator(name));
     		printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23_calculator(name));

		printf("a_muon\t\t\t\t%.3e\n\n",muon_gm2_calculator(name));

#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds\t\t%d\n",(higgsbounds_calculator(name)>1.));
#endif
 		printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));
		
		printf("charged_LSP\t\t\t%d\n\n",charged_LSP_calculator(name));
		
		flha_generator(name,"../output1.flha");
		printf("output1.flha generated\n\n");	
	}
	else printf("Invalid point\n\n");
	sprintf(name,"rm nuhm_softsusy%d.tmplha",getpid());
	system(name);		
#endif

#ifdef USE_ISAJET
	sprintf(name,"nuhm_isajet%d.tmplha",getpid());
	isajet_nuhm(m0, m12, tanb, A0, mu, mA, mtop, name);
	printf("---------------------------------------\n\n");
	printf("NUHM - SLHA file generated by ISAJET\n\n");

	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("Observable\t\t\tValue\n\n");

		printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
		printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

		printf("BR(Bs->mu mu)\t\t\t%.3e\n",Bsmumu_calculator(name));
		printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag_calculator(name));
		printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu_calculator(name));
	
		printf("BR(B->K* mu mu)_low\t\t%.3e\n",BRobs_BKstarmumu_lowq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_low\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_low\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_low\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_low\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_low\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_low\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_low\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_low\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_low\t\t%.3e\n\n",AI_BKstarmumu_lowq2_calculator(name));
	
		printf("BR(B->K* mu mu)_high\t\t%.3e\n",BRobs_BKstarmumu_highq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_high\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_high\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_high\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_high\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_high\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_high\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_high\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_high\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_high\t\t%.3e\n\n",AI_BKstarmumu_highq2_calculator(name));

		printf("q0^2(AFB(B->K* mu mu))\t\t%.3e\n",A_BKstarmumu_zero_calculator(name));
		printf("q0^2(AI(B->K* mu mu))\t\t%.3e\n\n",AI_BKstarmumu_zero_calculator(name));


		printf("BR(B->Xs mu mu)_low\t\t%.3e\n",BRBXsmumu_lowq2_calculator(name));
		printf("BR(B->Xs mu mu)_high\t\t%.3e\n",BRBXsmumu_highq2_calculator(name));
		printf("q0^2(AFB(B->Xs mu mu)\t\t%.3e\n",A_BXsmumu_zero_calculator(name));
		printf("BR(B->Xs tau tau)_high\t\t%.3e\n\n",BRBXstautau_highq2_calculator(name));
	
		printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu_calculator(name));
      		printf("R(B->tau nu)\t\t\t%.3e\n",RBtaunu_calculator(name));
      		printf("BR(B->D tau nu)\t\t\t%.3e\n",BDtaunu_calculator(name));
      		printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n",BDtaunu_BDenu_calculator(name));
     		printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu_calculator(name));
     		printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu_calculator(name));
     		printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu_calculator(name));
      		printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu_calculator(name));
     		printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23_calculator(name));

		printf("a_muon\t\t\t\t%.3e\n\n",muon_gm2_calculator(name));

#ifdef USE_HIGGSBOUNDS
		printf("excluded_HiggsBounds\t\t%d\n",(higgsbounds_calculator(name)>1.));
#endif
 		printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));
		
		printf("charged_LSP\t\t\t%d\n\n",charged_LSP_calculator(name));
		
		flha_generator(name,"../output2.flha");
		printf("output2.flha generated\n\n");	
	}
	else printf("Invalid point\n\n");
	sprintf(name,"rm nuhm_isajet%d.tmplha",getpid());
	system(name);	
#endif

	return 1;
}
