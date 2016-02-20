#include "include.h"

void Init_param(struct parameters* param)
/* initializes the "param" structure by initializing the parameters with either 0 or a value from the PDG2008 */
{
	int ie,je;
	
	param->model=-3; /* this parameter is used to test whether the scan of the SLHA file succeeds. */
	param->generator=0;
	param->Q=0.;
	param->m0=0.;
	param->m12=0.;
	param->tan_beta=0.;
	param->sign_mu=0.;
	param->A0=0.;
	param->mass_W=0.;
	param->Lambda=0.;
	param->Mmess=0.;
	param->N5=0.;
	param->cgrav=0.;
	param->m32=0.;
	param->mass_Z=0.;
	param->mass_b=0.;
	param->mass_top_pole=0.;
	param->mass_tau_pole=0.;
	param->inv_alpha_em=0.;
	param->alphas_MZ=0.;
	param->alpha=0.;
	param->Gfermi=0.;
	param->GAUGE_Q=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++)
	{
		param->charg_Umix[ie][je]=0.;
		param->charg_Vmix[ie][je]=0.;
		param->stop_mix[ie][je]=0.;
		param->sbot_mix[ie][je]=0.;
		param->stau_mix[ie][je]=0.;
	}
	for(ie=1;ie<=5;ie++) for(je=1;je<=5;je++) param->neut_mix[ie][je]=0.;
	for(ie=1;ie<=5;ie++) param->mass_neut[ie]=0.;
	for(ie=1;ie<=3;ie++) param->yub[ie]=param->yut[ie]=param->yutau[ie]=0.;
		
	param->Min=0.;
	param->M1_Min=0.;
	param->M2_Min=0.;
	param->M3_Min=0.;
	param->At_Min=0.;
	param->Ab_Min=0.;
	param->Atau_Min=0.;
	param->M2H1_Min=0.;
	param->M2H2_Min=0.;
	param->mu_Min=0.;
	param->M2A_Min=0.;
	param->tb_Min=0.;
	param->mA_Min=0.;
	param->MeL_Min=0.;
	param->MmuL_Min=0.;
	param->MtauL_Min=0.;
	param->MeR_Min=0.;
	param->MmuR_Min=0.;
	param->MtauR_Min=0.;
	param->MqL1_Min=0.;
	param->MqL2_Min=0.;
	param->MqL3_Min=0.;
	param->MuR_Min=0.;
	param->McR_Min=0.;
	param->MtR_Min=0.;
	param->MdR_Min=0.;
	param->MsR_Min=0.;
	param->MbR_Min=0.;
	param->N51=0.;
	param->N52=0.;
	param->N53=0.;
	param->M2H1_Q=0.;
	param->M2H2_Q=0.;
	param->mass_h0=0.;
	param->mass_H0=0.;
	param->mass_A0=0.;
	param->mass_H=0.;
	param->mass_dnl=0.;
	param->mass_upl=0.;
	param->mass_stl=0.;
	param->mass_chl=0.;
	param->mass_b1=0.;
	param->mass_t1=0.;
	param->mass_el=0.;
	param->mass_nuel=0.;
	param->mass_mul=0.;
	param->mass_numl=0.;
	param->mass_tau1=0.;
	param->mass_nutl=0.;
	param->mass_gluino=0.;
	param->mass_cha1=0.;
	param->mass_cha2=0.;
	param->mass_dnr=0.;
	param->mass_upr=0.;
	param->mass_str=0.;
	param->mass_chr=0.;
	param->mass_b2=0.;
	param->mass_t2=0.;
	param->mass_er=0.;
	param->mass_mur=0.;
	param->mass_tau2=0.;
	param->gp=0.;
	param->g2=0.;
	param->g3=0.;
	param->YU_Q=0.;
	param->YD_Q=0.;
	param->YE_Q=0.;
	param->HMIX_Q=0.;
	param->mu_Q=0.;
	param->tanb_GUT=0.;
	param->Higgs_VEV=0.;
	param->mA2_Q=0.;
	param->MSOFT_Q=0.;
	param->M1_Q=0.;
	param->M2_Q=0.;
	param->M3_Q=0.;
	param->MeL_Q=0.;
	param->MmuL_Q=0.;
	param->MtauL_Q=0.;
	param->MeR_Q=0.;
	param->MmuR_Q=0.;
	param->MtauR_Q=0.;
	param->MqL1_Q=0.;
	param->MqL2_Q=0.;
	param->MqL3_Q=0.;
	param->MuR_Q=0.;
	param->McR_Q=0.;
	param->MtR_Q=0.;
	param->MdR_Q=0.;
	param->MsR_Q=0.;
	param->MbR_Q=0.;
	param->AU_Q=0.;
	param->A_u=0.;
	param->A_c=0.;
	param->A_t=0.;
	param->AD_Q=0.;
	param->A_d=0.;
	param->A_s=0.;
	param->A_b=0.;
	param->AE_Q=0.;
	param->A_e=0.;
	param->A_mu=0.;
	param->A_tau=0.;
	param->mass_graviton=0.;
	param->mass_gravitino=0.;
	param->mass_nuer=0.;
	param->mass_numr=0.;
	param->mass_nutr=0.;
	param->mass_t=0.;
	param->mass_tau=0.;
	param->mass_gluon=0.;
	param->mass_nue=0.;
	param->mass_num=0.;
	param->mass_nut=0.;
	param->mass_photon=0.;
	param->mass_Z0=0.;

	/* SLHA2 */
	param->NMSSM=0;
	param->Rparity=0;
	param->CPviolation=0;
	param->Flavor=0;
	param->mass_nutau2=0.;
	param->mass_e2=0.;
	param->mass_nue2=0.;
	param->mass_mu2=0.;
	param->mass_numu2=0.;
	param->mass_d2=0.;
	param->mass_u2=0.;
	param->mass_s2=0.;
	param->mass_c2=0.;
	param->CKM_lambda=0.;
	param->CKM_A=0.;
	param->CKM_rho=0.;
	param->CKM_eta=0.;
	param->PMNS_theta12=0.;
	param->PMNS_theta23=0.;
	param->PMNS_theta13=0.;
	param->PMNS_delta13=0.;
	param->PMNS_alpha1=0.;
	param->PMNS_alpha2=0.;
	param->lambdaNMSSM_Min=0.;
	param->kappaNMSSM_Min=0.;
	param->AlambdaNMSSM_Min=0.;
	param->AkappaNMSSM_Min=0.;
	param->lambdaSNMSSM_Min=0.;
	param->xiFNMSSM_Min=0.;
	param->xiSNMSSM_Min=0.;
	param->mupNMSSM_Min=0.;
	param->mSp2NMSSM_Min=0.;
	param->mS2NMSSM_Min=0.;
	param->mass_H03=0.;
	param->mass_A02=0.;
	param->NMSSMRUN_Q=0.;
	param->lambdaNMSSM=0.;
	param->kappaNMSSM=0.;
	param->AlambdaNMSSM=0.;
	param->AkappaNMSSM=0.;
	param->lambdaSNMSSM=0.;
	param->xiFNMSSM=0.;
	param->xiSNMSSM=0.;
	param->mupNMSSM=0.;
	param->mSp2NMSSM=0.;
	param->mS2NMSSM=0.;
	param->PMNSU_Q=0.;
	param->CKM_Q=0.;
	param->MSE2_Q=0.;
	param->MSU2_Q=0.;
	param->MSD2_Q=0.;
	param->MSL2_Q=0.;
	param->MSQ2_Q=0.;
	param->TU_Q=0.;
	param->TD_Q=0.;
	param->TE_Q=0.;
	
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		param->sU_mix[ie][je]=0.;
		param->sD_mix[ie][je]=0.;
		param->sE_mix[ie][je]=0.;
	}
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		param->H0_mix[ie][je]=0.;
		param->sNU_mix[ie][je]=0.;
		param->sCKM_msq2[ie][je]=0.;
		param->sCKM_msl2[ie][je]=0.;
		param->sCKM_msd2[ie][je]=0.;
		param->sCKM_msu2[ie][je]=0.;
		param->sCKM_mse2[ie][je]=0.;
		param->PMNS_U[ie][je]=0.;
		param->CKM[ie][je]=0.;
		param->TU[ie][je]=0.;
		param->TD[ie][je]=0.;
		param->TE[ie][je]=0.;
	}
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=2;je++) param->A0_mix[ie][je]=0.;
	
	/* non-SLHA*/
	param->mass_b_1S=0.;
	param->mass_b_pole=0.;
	param->mtmt=0.;
	param->Lambda5=0.;
	
	/* Flavor physics */
	param->f_B=0.200;
	param->f_Bs=0.245;
	param->f_Ds=0.241;
	param->m_B=5.2795;
	param->m_Bs=5.3663;
	param->m_K=0.4937;
	param->m_Kstar=0.8917;
	param->m_D=1.86484;
	param->m_Ds=1.96849;
	param->life_B=1.638e-12;
	param->life_Bs=1.425e-12;
	param->life_Ds=5.e-13;
	
	/* masses and coupling from PDG 2008 */
	param->mass_u = 2.55e-3;
	param->mass_d = 5.04e-3;
	param->mass_s = 0.104;
	param->mass_c = 1.27;
	param->mass_b = 4.20;
	param->mass_top_pole = 172.4; /* from arXiv:0808.1089 */
	
	param->mass_e = 0.511e-3;
	param->mass_mu= 0.106;
	param->mass_tau_pole=1.7769;
	param->mass_tau=param->mass_tau_pole;
	
	param->mass_Z=91.1876;
	param->alphas_MZ=0.1176;
	param->mass_W=80.403;
	
	return;
}

/*--------------------------------------------------------------------*/

int Les_Houches_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[500];
	float dum;
	int ie,je;
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MODSEL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")))
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1:	fscanf(lecture,"%d",&param->model); break;
					case 3:	fscanf(lecture,"%d",&param->NMSSM); break;
					case 4:	fscanf(lecture,"%d",&param->Rparity); break;
					case 5:	fscanf(lecture,"%d",&param->CPviolation); break;
					case 6:	fscanf(lecture,"%d",&param->Flavor); break;
					case 12: fscanf(lecture,"%f",&param->Q); break;
				}	
			}
		}
	}
	fclose(lecture);		

	if(param->NMSSM != 0) param->model=-2;
	if(param->Rparity != 0) param->model=-2;
	if(param->CPviolation != 0) param->model=-2;
	if(param->Flavor != 0) param->model=-2;
	
	if(param->model<0) return 0;	
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"SPINFO"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 	fscanf(lecture,"%s",dummy); 
							if(!strncasecmp(dummy,"ISA",3)) param->generator=1; 
							if(!strncasecmp(dummy,"SOFTSUSY",8)) param->generator=2; 
							break;
					case 4: param->model=-1; fclose(lecture); return 0;
				}
			}	
		}
		else if(!strcasecmp(dummy,"SMINPUTS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->inv_alpha_em); break;
					case 2: fscanf(lecture,"%f",&param->Gfermi); break;
					case 3: fscanf(lecture,"%f",&param->alphas_MZ); break;
					case 4: fscanf(lecture,"%f",&param->mass_Z); break;
					case 5: fscanf(lecture,"%f",&param->mass_b); break;
					case 6: fscanf(lecture,"%f",&param->mass_top_pole); break;
					case 7: fscanf(lecture,"%f",&param->mass_tau_pole); break;
					case 8: fscanf(lecture,"%f",&param->mass_nutau2); break;
					case 11: fscanf(lecture,"%f",&param->mass_e2); break;
					case 12: fscanf(lecture,"%f",&param->mass_nue2); break;
					case 13: fscanf(lecture,"%f",&param->mass_mu2); break;
					case 14: fscanf(lecture,"%f",&param->mass_numu2); break;
					case 21: fscanf(lecture,"%f",&param->mass_d2); break;
					case 22: fscanf(lecture,"%f",&param->mass_u2); break;
					case 23: fscanf(lecture,"%f",&param->mass_s2); break;
					case 24: fscanf(lecture,"%f",&param->mass_c2); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"VCKMIN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->CKM_lambda); break;
					case 2: fscanf(lecture,"%f",&param->CKM_A); break;
					case 3: fscanf(lecture,"%f",&param->CKM_rho); break;
					case 4: fscanf(lecture,"%f",&param->CKM_eta); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"UPMNSIN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->PMNS_theta12); break;
					case 2: fscanf(lecture,"%f",&param->PMNS_theta23); break;
					case 3: fscanf(lecture,"%f",&param->PMNS_theta13); break;
					case 4: fscanf(lecture,"%f",&param->PMNS_delta13); break;
					case 5: fscanf(lecture,"%f",&param->PMNS_alpha1); break;
					case 6: fscanf(lecture,"%f",&param->PMNS_alpha2); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"MINPAR"))
		{
			switch(param->model)
			{
				case 1:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
							switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
							{
								case 1: fscanf(lecture,"%f",&param->m0); break;
								case 2: fscanf(lecture,"%f",&param->m12); break;
								case 3: fscanf(lecture,"%f",&param->tan_beta); break;
								case 4: fscanf(lecture,"%f",&param->sign_mu); break;
								case 5: fscanf(lecture,"%f",&param->A0); break;
							}
					}
					break;
				}
	
				case 2:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->Lambda); break;
							case 2: fscanf(lecture,"%f",&param->Mmess); break;
							case 3: fscanf(lecture,"%f",&param->tan_beta); break;
							case 4: fscanf(lecture,"%f",&param->sign_mu); break;
							case 5: fscanf(lecture,"%f",&param->N5); break;
							case 6: fscanf(lecture,"%f",&param->cgrav); break;
						}
					}
					break;
				}
	
				case 3:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->m32); break;
							case 2: fscanf(lecture,"%f",&param->m0); break;
							case 3: fscanf(lecture,"%f",&param->tan_beta); break;
							case 4: fscanf(lecture,"%f",&param->sign_mu); break;
						}
					}
					break;
				}

				default:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->tan_beta); break;
						}
					}
					break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"EXTPAR"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 0: fscanf(lecture,"%f",&param->Min); break;
					case 1: fscanf(lecture,"%f",&param->M1_Min); break;
					case 2: fscanf(lecture,"%f",&param->M2_Min); break;
					case 3: fscanf(lecture,"%f",&param->M3_Min); break;	
					case 11: fscanf(lecture,"%f",&param->At_Min); break;
					case 12: fscanf(lecture,"%f",&param->Ab_Min); break;
					case 13: fscanf(lecture,"%f",&param->Atau_Min); break;
					case 21: fscanf(lecture,"%f",&param->M2H1_Min); break;
					case 22: fscanf(lecture,"%f",&param->M2H2_Min); break;
					case 23: fscanf(lecture,"%f",&param->mu_Min); break;
					case 24: fscanf(lecture,"%f",&param->M2A_Min); break;
					case 25: fscanf(lecture,"%f",&param->tb_Min); break;
					case 26: fscanf(lecture,"%f",&param->mA_Min); break;
					case 31: fscanf(lecture,"%f",&param->MeL_Min); break;
					case 32: fscanf(lecture,"%f",&param->MmuL_Min); break;
					case 33: fscanf(lecture,"%f",&param->MtauL_Min); break;
					case 34: fscanf(lecture,"%f",&param->MeR_Min); break;
					case 35: fscanf(lecture,"%f",&param->MmuR_Min); break;
					case 36: fscanf(lecture,"%f",&param->MtauR_Min); break;
					case 41: fscanf(lecture,"%f",&param->MqL1_Min); break;
					case 42: fscanf(lecture,"%f",&param->MqL2_Min); break;
					case 43: fscanf(lecture,"%f",&param->MqL3_Min); break;
					case 44: fscanf(lecture,"%f",&param->MuR_Min); break;
					case 45: fscanf(lecture,"%f",&param->McR_Min); break;
					case 46: fscanf(lecture,"%f",&param->MtR_Min); break;
					case 47: fscanf(lecture,"%f",&param->MdR_Min); break;
					case 48: fscanf(lecture,"%f",&param->MsR_Min); break;
					case 49: fscanf(lecture,"%f",&param->MbR_Min); break;
					case 51: fscanf(lecture,"%f",&param->N51); break;
					case 52: fscanf(lecture,"%f",&param->N52); break;
					case 53: fscanf(lecture,"%f",&param->N53); break;
					case 61: fscanf(lecture,"%f",&param->lambdaNMSSM_Min); break;
					case 62: fscanf(lecture,"%f",&param->kappaNMSSM_Min); break;
					case 63: fscanf(lecture,"%f",&param->AlambdaNMSSM_Min); break;
					case 64: fscanf(lecture,"%f",&param->AkappaNMSSM_Min); break;
					case 65: fscanf(lecture,"%f",&param->lambdaSNMSSM_Min); break;
					case 66: fscanf(lecture,"%f",&param->xiFNMSSM_Min); break;
					case 67: fscanf(lecture,"%f",&param->xiSNMSSM_Min); break;
					case 68: fscanf(lecture,"%f",&param->mupNMSSM_Min); break;
					case 69: fscanf(lecture,"%f",&param->mSp2NMSSM_Min); break;
					case 70: fscanf(lecture,"%f",&param->mS2NMSSM_Min); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->mass_d); break;
					case 2: fscanf(lecture,"%f",&param->mass_u); break;
					case 3: fscanf(lecture,"%f",&param->mass_s); break;
					case 4: fscanf(lecture,"%f",&param->mass_c); break;
					case 6: fscanf(lecture,"%f",&param->mass_t); break;
					case 11: fscanf(lecture,"%f",&param->mass_e); break;
					case 12: fscanf(lecture,"%f",&param->mass_nue); break;
					case 13: fscanf(lecture,"%f",&param->mass_mu); break;
					case 14: fscanf(lecture,"%f",&param->mass_num); break;
					case 15: fscanf(lecture,"%f",&param->mass_tau); break;
					case 16: fscanf(lecture,"%f",&param->mass_nut); break;
					case 21: fscanf(lecture,"%f",&param->mass_gluon); break;
					case 22: fscanf(lecture,"%f",&param->mass_photon); break;
					case 23: fscanf(lecture,"%f",&param->mass_Z0); break;
					case 24: fscanf(lecture,"%f",&param->mass_W); break;
					case 25: fscanf(lecture,"%f",&param->mass_h0); break;
					case 35: fscanf(lecture,"%f",&param->mass_H0); break;
					case 36: fscanf(lecture,"%f",&param->mass_A0); break;
					case 37: fscanf(lecture,"%f",&param->mass_H); break;
					case 39: fscanf(lecture,"%f",&param->mass_graviton); break;
					case 45: fscanf(lecture,"%f",&param->mass_H03); break;
					case 46: fscanf(lecture,"%f",&param->mass_A02); break;
					case 1000001: fscanf(lecture,"%f",&param->mass_dnl); break;
					case 1000002: fscanf(lecture,"%f",&param->mass_upl); break;
					case 1000003: fscanf(lecture,"%f",&param->mass_stl); break;
					case 1000004: fscanf(lecture,"%f",&param->mass_chl); break;
					case 1000005: fscanf(lecture,"%f",&param->mass_b1); break;
					case 1000006: fscanf(lecture,"%f",&param->mass_t1); break;
					case 1000011: fscanf(lecture,"%f",&param->mass_el); break;
					case 1000012: fscanf(lecture,"%f",&param->mass_nuel); break;
					case 1000013: fscanf(lecture,"%f",&param->mass_mul); break;
					case 1000014: fscanf(lecture,"%f",&param->mass_numl); break;
					case 1000015: fscanf(lecture,"%f",&param->mass_tau1); break;
					case 1000016: fscanf(lecture,"%f",&param->mass_nutl); break;
					case 1000021: fscanf(lecture,"%f",&param->mass_gluino); break;
					case 1000022: fscanf(lecture,"%f",&param->mass_neut[1]); break;
					case 1000023: fscanf(lecture,"%f",&param->mass_neut[2]); break;
					case 1000024: fscanf(lecture,"%f",&param->mass_cha1); break;
					case 1000025: fscanf(lecture,"%f",&param->mass_neut[3]); break;
					case 1000035: fscanf(lecture,"%f",&param->mass_neut[4]); break;
					case 1000037: fscanf(lecture,"%f",&param->mass_cha2); break;
					case 1000039: fscanf(lecture,"%f",&param->mass_gravitino); break;
					case 1000045: fscanf(lecture,"%f",&param->mass_neut[5]); break;
					case 2000001: fscanf(lecture,"%f",&param->mass_dnr); break;
					case 2000002: fscanf(lecture,"%f",&param->mass_upr); break;
					case 2000003: fscanf(lecture,"%f",&param->mass_str); break;
					case 2000004: fscanf(lecture,"%f",&param->mass_chr); break;
					case 2000005: fscanf(lecture,"%f",&param->mass_b2); break;
					case 2000006: fscanf(lecture,"%f",&param->mass_t2); break;
					case 2000011: fscanf(lecture,"%f",&param->mass_er); break;
					case 2000012: fscanf(lecture,"%f",&param->mass_nuer); break;
					case 2000013: fscanf(lecture,"%f",&param->mass_mur); break;
					case 2000014: fscanf(lecture,"%f",&param->mass_numr); break;
					case 2000015: fscanf(lecture,"%f",&param->mass_tau2); break;
					case 2000016: fscanf(lecture,"%f",&param->mass_nutr); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
		}
		else if(!strcasecmp(dummy,"STOPMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->stop_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"SBOTMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sbot_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"STAUMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->stau_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"NMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->neut_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"NMNMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->neut_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"UMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));

				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->charg_Umix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"VMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->charg_Vmix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"GAUGE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->GAUGE_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->gp); break;
					case 2: fscanf(lecture,"%f",&param->g2); break;
					case 3: fscanf(lecture,"%f",&param->g3); break;	
				}
			}
		}
		else if(!strcasecmp(dummy,"YU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->YU_Q);
				else 				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%f",&param->yut[ie]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"YD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->YD_Q);
				else
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%f",&param->yub[ie]);
					}
				}			
			}		
		}
		else if(!strcasecmp(dummy,"YE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->YE_Q);
				else
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%f",&param->yutau[ie]);
					}
				}			
			}		
		}
		else if(!strcasecmp(dummy,"HMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->HMIX_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->mu_Q); break;
					case 2: fscanf(lecture,"%f",&param->tanb_GUT); break;
					case 3: fscanf(lecture,"%f",&param->Higgs_VEV); break;	
					case 4: fscanf(lecture,"%f",&param->mA2_Q); break;
				}
			}
		}
		else if(!strcasecmp(dummy,"NMHMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->H0_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"NMAMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->A0_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"MSOFT"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->MSOFT_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->M1_Q); break;
					case 2: fscanf(lecture,"%f",&param->M2_Q); break;
					case 3: fscanf(lecture,"%f",&param->M3_Q); break;	
					case 21: fscanf(lecture,"%f",&param->M2H1_Q); break;
					case 22: fscanf(lecture,"%f",&param->M2H2_Q); break;
					case 31: fscanf(lecture,"%f",&param->MeL_Q); break;
					case 32: fscanf(lecture,"%f",&param->MmuL_Q); break;
					case 33: fscanf(lecture,"%f",&param->MtauL_Q); break;
					case 34: fscanf(lecture,"%f",&param->MeR_Q); break;
					case 35: fscanf(lecture,"%f",&param->MmuR_Q); break;
					case 36: fscanf(lecture,"%f",&param->MtauR_Q); break;
					case 41: fscanf(lecture,"%f",&param->MqL1_Q); break;
					case 42: fscanf(lecture,"%f",&param->MqL2_Q); break;
					case 43: fscanf(lecture,"%f",&param->MqL3_Q); break;
					case 44: fscanf(lecture,"%f",&param->MuR_Q); break;
					case 45: fscanf(lecture,"%f",&param->McR_Q); break;
					case 46: fscanf(lecture,"%f",&param->MtR_Q); break;
					case 47: fscanf(lecture,"%f",&param->MdR_Q); break;
					case 48: fscanf(lecture,"%f",&param->MsR_Q); break;
					case 49: fscanf(lecture,"%f",&param->MbR_Q); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"AU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->AU_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->A_u); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 2: fscanf(lecture,"%f",&param->A_c); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->A_t); break;	
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"AD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->AD_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->A_d); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 2: fscanf(lecture,"%f",&param->A_s); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->A_b); break;	
						}
						break;
					}
				}
			}
		}
		else if(!strcasecmp(dummy,"AE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->AE_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->A_e); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 2: fscanf(lecture,"%f",&param->A_mu); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->A_tau); break;	
						}
						break;
					}
				}
			}	
		}		
		else if(!strcasecmp(dummy,"NMSSMRUN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->NMSSMRUN_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: fscanf(lecture,"%f",&param->lambdaNMSSM); break;
					case 2: fscanf(lecture,"%f",&param->kappaNMSSM); break;
					case 3: fscanf(lecture,"%f",&param->AlambdaNMSSM); break;
					case 4: fscanf(lecture,"%f",&param->AkappaNMSSM); break;
					case 5: fscanf(lecture,"%f",&param->lambdaSNMSSM); break;
					case 6: fscanf(lecture,"%f",&param->xiFNMSSM); break;
					case 7: fscanf(lecture,"%f",&param->xiSNMSSM); break;
					case 8: fscanf(lecture,"%f",&param->mupNMSSM); break;
					case 9: fscanf(lecture,"%f",&param->mSp2NMSSM); break;
					case 10: fscanf(lecture,"%f",&param->mS2NMSSM); break;
				}
			}
		}
		else if(!strcasecmp(dummy,"USQMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sU_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"DSQMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sD_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"SELMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sE_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"SELMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sNU_mix[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"MSQ2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->MSQ2_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sCKM_msq2[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"MSL2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->MSL2_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sCKM_msl2[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"MSD2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->MSD2_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sCKM_msd2[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"MSU2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->MSU2_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sCKM_msu2[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"MSE2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->MSE2_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->sCKM_mse2[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"VCKM"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->CKM_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->CKM[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"UPMNS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->PMNSU_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->PMNS_U[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"TU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->TU_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->TU[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"TD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->TD_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->TD[ie][je]);
					}
				}			
			}	
		}
		else if(!strcasecmp(dummy,"TE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=MGUT=")) break;
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->TE_Q);
				else 

				if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%f",&param->TE[ie][je]);
					}
				}			
			}	
		}
	}
	fclose(lecture);
	if(param->model<0) return 0;	
 			
	slha_adjust(param);
	
	return 1;	
}

/*--------------------------------------------------------------------*/

void slha_adjust(struct parameters* param)
{
	float dum;
	float mass[7],invmat[7][7];
	int ie,je,iemax;

	if(param->mass_Z0==0.) param->mass_Z0=param->mass_Z;

	if(param->MSOFT_Q==0.) param->MSOFT_Q=max(param->TU_Q,max(param->TD_Q,max(param->TE_Q,max(param->PMNSU_Q,max(param->CKM_Q,max(param->MSE2_Q,max(param->MSU2_Q,max(param->MSD2_Q,max(param->MSL2_Q,max(param->MSQ2_Q,max(param->NMSSMRUN_Q,max(param->YU_Q,max(param->YD_Q,max(param->YE_Q,max(param->HMIX_Q,max(param->GAUGE_Q,max(param->AU_Q,max(param->AD_Q,param->AE_Q))))))))))))))))));
	
	if(param->Q==0.) param->Q=sqrt(param->mass_t1*param->mass_t2);
	
	if(param->MSOFT_Q==0.) param->TU_Q=param->TD_Q=param->TE_Q=param->PMNSU_Q=param->CKM_Q=param->MSE2_Q=param->MSU2_Q=param->MSD2_Q=param->MSL2_Q=param->MSQ2_Q=param->MSOFT_Q=param->NMSSMRUN_Q=param->YU_Q=param->YD_Q=param->YE_Q=param->HMIX_Q=param->GAUGE_Q=param->AU_Q=param->AD_Q=param->AE_Q=param->Q;
			
	if(param->tan_beta==0.) param->tan_beta=param->tanb_GUT;
	
	if((param->tan_beta*param->MSOFT_Q)==0.) param->model=-3;
	
	if(param->MeL_Q==0.) param->MeL_Q=param->mass_el;
	if(param->MmuL_Q==0.) param->MmuL_Q=param->mass_mul;
	if(param->MtauL_Q==0.) param->MtauL_Q=param->mass_tau1;
	if(param->MeR_Q==0.) param->MeR_Q=param->mass_er;
	if(param->MmuR_Q==0.) param->MmuR_Q=param->mass_mur;
	if(param->MtauR_Q==0.) param->MtauR_Q=param->mass_tau2;
	if(param->MqL1_Q==0.) param->MqL1_Q=param->mass_dnl;
	if(param->MqL2_Q==0.) param->MqL2_Q=param->mass_stl;
	if(param->MqL3_Q==0.) param->MqL3_Q=param->mass_b1;
	if(param->MuR_Q==0.) param->MuR_Q=param->mass_upr;
	if(param->McR_Q==0.) param->McR_Q=param->mass_chr;
	if(param->MtR_Q==0.) param->MtR_Q=param->mass_t2;
	if(param->MdR_Q==0.) param->MdR_Q=param->mass_dnr;
	if(param->MsR_Q==0.) param->MsR_Q=param->mass_str;
	if(param->MbR_Q==0.) param->MbR_Q=param->mass_b2;
		
	if(param->A_tau==0.) param->A_tau=param->TE[3][3];
	if(param->A_b==0.) param->A_b=param->TD[3][3];
	if(param->A_t==0.) param->A_t=param->TU[3][3];
	
	if(param->stop_mix[1][1]==0.) 
	{
		mass[1]=param->mass_upl;
		mass[2]=param->mass_chl;
		mass[3]=param->mass_t1;
		mass[4]=param->mass_upr;
		mass[5]=param->mass_chr;
		mass[6]=param->mass_t2;

		for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) invmat[ie][je]=param->sU_mix[je][ie];
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[1][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[1][ie]);
		}
		param->mass_upl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[2][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[2][ie]);
		}
		param->mass_chl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[3][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[3][ie]);
		}
		param->mass_t2=mass[iemax];		
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[4][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[4][ie]);
		}
		param->mass_upr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[5][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[5][ie]);
		}
		param->mass_chr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[6][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[6][ie]);
		}
		param->mass_t1=mass[iemax];

		param->stop_mix[1][1]=param->sU_mix[1][6];
		param->stop_mix[1][2]=param->sU_mix[1][3];
	}	
			
	dum=atan(param->stop_mix[1][2]/param->stop_mix[1][1]);
	
	if(param->generator==1) dum=atan(param->stop_mix[2][1]/param->stop_mix[1][1]);
	
	param->stop_mix[1][1]=cos(dum);
	param->stop_mix[2][1]=-sin(dum);
	param->stop_mix[1][2]=sin(dum);
	param->stop_mix[2][2]=cos(dum);
		
	if(param->sbot_mix[1][1]==0.)
	{
		mass[1]=param->mass_dnl;
		mass[2]=param->mass_stl;
		mass[3]=param->mass_b1;
		mass[4]=param->mass_dnr;
		mass[5]=param->mass_str;
		mass[6]=param->mass_b2;

		for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) invmat[ie][je]=param->sD_mix[je][ie];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[1][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[1][ie]);
		}
		param->mass_dnl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[2][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[2][ie]);
		}
		param->mass_stl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[3][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[3][ie]);
		}
		param->mass_b1=mass[iemax];		

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[4][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[4][ie]);
		}
		param->mass_dnr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[5][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[5][ie]);
		}
		param->mass_str=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[6][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[6][ie]);
		}
		param->mass_b2=mass[iemax];

		param->sbot_mix[1][1]=param->sD_mix[1][6];
		param->sbot_mix[1][2]=param->sD_mix[1][3];
	}
		
	dum=atan(param->sbot_mix[1][2]/param->sbot_mix[1][1]);

	if(param->generator==1) dum=atan(param->sbot_mix[2][1]/param->sbot_mix[1][1]);

	param->sbot_mix[1][1]=cos(dum);
	param->sbot_mix[2][1]=-sin(dum);
	param->sbot_mix[1][2]=sin(dum);
	param->sbot_mix[2][2]=cos(dum);
		
	if(param->stau_mix[1][1]==0.)
	{
		mass[1]=param->mass_el;
		mass[2]=param->mass_mul;
		mass[3]=param->mass_tau1;
		mass[4]=param->mass_er;
		mass[5]=param->mass_mur;
		mass[6]=param->mass_tau2;

		for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) invmat[ie][je]=param->sE_mix[je][ie];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[1][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[1][ie]);
		}
		param->mass_el=mass[iemax];

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[2][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[2][ie]);
		}
		param->mass_mul=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[3][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[3][ie]);
		}
		param->mass_tau2=mass[iemax];		
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[4][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[4][ie]);
		}
		param->mass_er=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[5][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[5][ie]);
		}
		param->mass_mur=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[6][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[6][ie]);
		}
		param->mass_tau1=mass[iemax];

		param->stau_mix[1][1]=param->sE_mix[1][6];
		param->stau_mix[1][2]=param->sE_mix[1][3];
	}
		
	dum=atan(param->stau_mix[1][2]/param->stau_mix[1][1]);

	if(param->generator==1) dum=atan(param->stau_mix[2][1]/param->stau_mix[1][1]);

	param->stau_mix[1][1]=cos(dum);
	param->stau_mix[2][1]=-sin(dum);
	param->stau_mix[1][2]=sin(dum);
	param->stau_mix[2][2]=cos(dum);

	
 	param->mass_b_pole=mb_pole(param);
 	param->mass_b_1S=mb_1S(param);
	
	param->mtmt=mt_mt(param);
	
	return;
}

/*--------------------------------------------------------------------*/

int test_slha(char name[])
/* "container" function scanning the SLHA file "name" and checking if it is valid. A negative value indicates a problem with the SLHA file. */
{
	struct parameters param;
		
	Init_param(&param);
	if(Les_Houches_Reader(name,&param)) return 1;

	return param.model;
}
