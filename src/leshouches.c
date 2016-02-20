#include "include.h"

void Init_param(struct parameters* param)
/* initializes the "param" structure by initializing the parameters with either 0 or a value from the PDG2006 */
{
	int ie,je;
	
	param->model=-1; /* this parameter is used to test if the scan of the SLHA file succeeds. */
	param->generator=0;
	param->Q=0.;
	param->m0=0.;
	param->m1_2=0.;
	param->tan_beta=0.;
	param->sign_mu=0.;
	param->A0=0.;
	param->mass_W=0.;
	param->Lambda=0.;
	param->Mmess=0.;
	param->N5=0.;
	param->cgrav=0.;
	param->m3_2=0.;
	param->mass_Z=0.;
	param->mass_b=0.;
	param->mass_top_pole=0.;
	param->mass_tau_pole=0.;
	param->inv_alpha_em=0.;
	param->alpha_s_MZ=0.;
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
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) param->neut_mix[ie][je]=0.;
	for(ie=1;ie<=4;ie++) param->mass_neut[ie]=0.;
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
	param->mass_cH0=0.;
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
	param->yut=0.;
	param->YD_Q=0.;
	param->yub=0.;
	param->YE_Q=0.;
	param->yutau=0.;
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
	param->mass_gluon=0.;
	param->mass_nue=0.;
	param->mass_num=0.;
	param->mass_nut=0.;
	param->mass_photon=0.;
	param->mass_Z0=0.;
	param->mass_b_1S = 0.;

	/* masses and coupling from PDG 2006 */
	param->mass_u = 2.e-3;
	param->mass_d = 5.e-3;
	param->mass_s = 0.1;
	param->mass_c = 1.2;
	param->mass_b = 4.20;
	param->mass_top_pole = 172.5;
	
	param->mass_e = 0.511e-3;
	param->mass_mu= 0.106;
	param->mass_tau_pole=1776.9;
	
	param->mass_Z=91.1876;
	param->alpha_s_MZ=0.1172;
	
	return;
}

/*--------------------------------------------------------------------*/

int Les_Houches_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[50];
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
					case 12: fscanf(lecture,"%f",&param->Q); break;
				}	
			}
		}
	}
	fclose(lecture);	
	
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
							if(!strncasecmp(dummy,"ISAJET",6)) param->generator=1; 
							if(!strncasecmp(dummy,"SOFTSUSY",8)) param->generator=2; 
							break;
					case 4: param->model==-1; fclose(lecture); return 0;
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
					case 3: fscanf(lecture,"%f",&param->alpha_s_MZ); break;
					case 4: fscanf(lecture,"%f",&param->mass_Z); break;
					case 5: fscanf(lecture,"%f",&param->mass_b); break;
					case 6: fscanf(lecture,"%f",&param->mass_top_pole); break;
					case 7: fscanf(lecture,"%f",&param->mass_tau_pole); break;
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
								case 2: fscanf(lecture,"%f",&param->m1_2); break;
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
							case 1: fscanf(lecture,"%f",&param->m3_2); break;
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
					case 35: fscanf(lecture,"%f",&param->mass_cH0); break;
					case 36: fscanf(lecture,"%f",&param->mass_A0); break;
					case 37: fscanf(lecture,"%f",&param->mass_H); break;
					case 39: fscanf(lecture,"%f",&param->mass_graviton); break;
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
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->stop_mix[1][1]); break;
							case 2: fscanf(lecture,"%f",&param->stop_mix[1][2]); break;
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->stop_mix[2][1]); break;
							case 2: fscanf(lecture,"%f",&param->stop_mix[2][2]); break;
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"SBOTMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->sbot_mix[1][1]); break;
							case 2: fscanf(lecture,"%f",&param->sbot_mix[1][2]); break;
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->sbot_mix[2][1]); break;
							case 2: fscanf(lecture,"%f",&param->sbot_mix[2][2]); break;
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"STAUMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->stau_mix[1][1]); break;
							case 2: fscanf(lecture,"%f",&param->stau_mix[1][2]); break;
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->stau_mix[2][1]); break;
							case 2: fscanf(lecture,"%f",&param->stau_mix[2][2]); break;
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"NMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->neut_mix[1][1]); break;
							case 2: fscanf(lecture,"%f",&param->neut_mix[1][2]); break;
							case 3: fscanf(lecture,"%f",&param->neut_mix[1][3]); break;
							case 4: fscanf(lecture,"%f",&param->neut_mix[1][4]); break;
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->neut_mix[2][1]); break;
							case 2: fscanf(lecture,"%f",&param->neut_mix[2][2]); break;
							case 3: fscanf(lecture,"%f",&param->neut_mix[2][3]); break;
							case 4: fscanf(lecture,"%f",&param->neut_mix[2][4]); break;
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->neut_mix[3][1]); break;
							case 2: fscanf(lecture,"%f",&param->neut_mix[3][2]); break;
							case 3: fscanf(lecture,"%f",&param->neut_mix[3][3]); break;
							case 4: fscanf(lecture,"%f",&param->neut_mix[3][4]); break;
						}
						break;
					}
					case 4: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->neut_mix[4][1]); break;
							case 2: fscanf(lecture,"%f",&param->neut_mix[4][2]); break;
							case 3: fscanf(lecture,"%f",&param->neut_mix[4][3]); break;
							case 4: fscanf(lecture,"%f",&param->neut_mix[4][4]); break;
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"UMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->charg_Umix[1][1]); break;
							case 2: fscanf(lecture,"%f",&param->charg_Umix[1][2]); break;
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->charg_Umix[2][1]); break;
							case 2: fscanf(lecture,"%f",&param->charg_Umix[2][2]); break;
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"VMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->charg_Vmix[1][1]); break;
							case 2: fscanf(lecture,"%f",&param->charg_Vmix[1][2]); break;
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 1: fscanf(lecture,"%f",&param->charg_Vmix[2][1]); break;
							case 2: fscanf(lecture,"%f",&param->charg_Vmix[2][2]); break;
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"GAUGE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
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
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->YU_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->yut); break;	
						}
						break;
					}
				}
			}	
		}
		else if(!strcasecmp(dummy,"YD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->YD_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->yub); break;	
						}
						break;
					}
				}
			}		
		}
		else if(!strcasecmp(dummy,"YE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%f",&param->YE_Q);
				else switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
						{
							case 3: fscanf(lecture,"%f",&param->yutau); break;	
						}
						break;
					}
				}
			}		
		}
		else if(!strcasecmp(dummy,"HMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
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
		else if(!strcasecmp(dummy,"MSOFT"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
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
	}
	fclose(lecture);

	if(param->MSOFT_Q==0.) param->MSOFT_Q=max(param->YU_Q,max(param->YD_Q,max(param->YE_Q,max(param->HMIX_Q,max(param->GAUGE_Q,max(param->AU_Q,max(param->AD_Q,param->AE_Q)))))));
	
	if(param->Q==0.) param->Q=sqrt(param->mass_t1*param->mass_t2);
	
	if(param->MSOFT_Q==0.) param->MSOFT_Q=param->YU_Q=param->YD_Q=param->YE_Q=param->HMIX_Q=param->GAUGE_Q=param->AU_Q=param->AD_Q=param->AE_Q=param->Q;
		
	if(param->tan_beta==0.) param->tan_beta=param->tanb_GUT;
	
	if((param->tan_beta*param->MSOFT_Q)==0.) param->model=-1;
	
	dum=acos(param->stop_mix[1][1]);
	param->stop_mix[2][1]=-sin(dum);
	param->stop_mix[1][2]=sin(dum);
	param->stop_mix[2][2]=param->stop_mix[1][1];
		
	dum=acos(param->sbot_mix[1][1]);
	param->sbot_mix[2][1]=-sin(dum);
	param->sbot_mix[1][2]=sin(dum);
	param->sbot_mix[2][2]=param->sbot_mix[1][1];		

	dum=acos(param->stau_mix[1][1]);
	param->stau_mix[2][1]=-sin(dum);
	param->stau_mix[1][2]=sin(dum);
	param->stau_mix[2][2]=param->stau_mix[1][1];

	dum=acos(param->stau_mix[1][1]);
	param->stau_mix[2][1]=-sin(dum);
	param->stau_mix[1][2]=sin(dum);
	param->stau_mix[2][2]=param->stau_mix[1][1];
	
 	param->mass_b_1S=b_mass_1S(param);
 	/* param->mass_b_1S=4.68; */
			
	if(param->model==-1) return 0;
	
	return 1;	
}
