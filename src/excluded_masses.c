#include "include.h"

int excluded_masses(struct parameters* param)
/* tests whether the SUSY point is excluded by the collider contraints */
/* if excluded, return 1, otherwise 0 */
{
	int excluded=0;

#ifdef SMONLY	
	return excluded;
#endif

	if(param->mass_h0!=0.) excluded=(excluded||(fabs(param->mass_h0)<111.)); /* Higgs */
	if(param->mass_H!=0.) excluded=(excluded||(fabs(param->mass_H)<79.3)); /* charged Higgs */
	if(param->mass_A0!=0.) excluded=(excluded||(fabs(param->mass_A0)<93.4)); /* CP-odd Higgs */

#ifdef SM_ChargedHiggs	
	return excluded;
#endif

	if(param->mass_neut[1]!=0.) excluded=(excluded||(fabs(param->mass_neut[1])<46.)); /* neutralino 1 */
	if(param->mass_neut[2]!=0.) if(param->tan_beta<40.) excluded=(excluded||(fabs(param->mass_neut[2])<62.4)); /* neutralino 2 */
	if(param->mass_neut[3]!=0.) if(param->tan_beta<40.) excluded=(excluded||(fabs(param->mass_neut[3])<99.9)); /* neutralino 3 */
	if(param->mass_neut[4]!=0.) if(param->tan_beta<40.) excluded=(excluded||(fabs(param->mass_neut[4])<116.)); /* neutralino 4 */
	if(param->mass_cha1!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_cha1)-fabs(param->mass_neut[1])>3.)) excluded=(excluded||(fabs(param->mass_cha1)<94.)); /* chargino */
	if(param->mass_er!=0.) excluded=(excluded||(fabs(param->mass_er)<73.)); /* slepton R */
	if(param->mass_mur!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_mur)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_mur)<94.)); /* slepton R */
	if(param->mass_nuel!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_er)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_nuel)<94.)); /* sneutrino */
	if(param->mass_numl!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_mur)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_numl)<94.)); /* sneutrino */
	if(param->mass_nutl!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_tau1)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_nutl)<94.)); /* sneutrino */
	if(param->mass_dnr!=0.) excluded=(excluded||(fabs(param->mass_dnr)<379.)); /* squark R */
	if(param->mass_upr!=0.) excluded=(excluded||(fabs(param->mass_upr)<379.)); /* squark R */
	if(param->mass_str!=0.) excluded=(excluded||(fabs(param->mass_str)<379.)); /* squark R */
	if(param->mass_chr!=0.) excluded=(excluded||(fabs(param->mass_chr)<379.)); /* squark R */
	if(param->mass_t1!=0.) if(fabs(param->mass_t1)-fabs(param->mass_neut[1])>10.) excluded=(excluded||(fabs(param->mass_t1)<95.7)); /* stop */	
	if(param->mass_gluino!=0.) excluded=(excluded||(fabs(param->mass_gluino)<308.)); /* gluino */
	if(param->mass_b1!=0.)  if(fabs(param->mass_b1)-fabs(param->mass_neut[1])>8.) excluded=(excluded||(fabs(param->mass_b1)<89.)); /* sbottom */
	if(param->mass_tau1!=0.) if(fabs(param->mass_tau1)-fabs(param->mass_neut[1])>15.) excluded=(excluded||(fabs(param->mass_tau1)<81.9)); /* stau */
	
	return excluded;
}

/*--------------------------------------------------------------------*/

int excluded_mass_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point is excluded by the mass contraints */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return excluded_masses(&param);
}

/*--------------------------------------------------------------------*/

int charged_LSP(struct parameters* param)
/* tests whether the SUSY point corresponds to a charged LSP (NLSP if the LSP is a gravitino) */
/* if the LSP is charged, return 1, otherwise 0 */
{
	double mass_non_chargee;
	int charged_LSP=0;

#ifdef SMONLY	
	return charged_LSP;
#endif

#ifdef SM_ChargedHiggs	
	return charged_LSP;
#endif

	mass_non_chargee=fabs(param->mass_neut[1]);
	
	if(param->mass_nuel!=0.) mass_non_chargee=min(fabs(param->mass_nuel),mass_non_chargee);
	if(param->mass_numl!=0.) mass_non_chargee=min(fabs(param->mass_numl),mass_non_chargee);
	if(param->mass_nutl!=0.) mass_non_chargee=min(fabs(param->mass_nutl),mass_non_chargee);
	if(param->mass_gluino!=0.) mass_non_chargee=min(fabs(param->mass_gluino),mass_non_chargee);
	
	if(param->mass_tau1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_tau1)<mass_non_chargee));
	if(param->mass_t1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_t1)<mass_non_chargee));
	if(param->mass_cha1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_cha1)<mass_non_chargee));
	
	return charged_LSP;
}

/*--------------------------------------------------------------------*/

int charged_LSP_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point corresponds to a charged LSP */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return charged_LSP(&param);
}

