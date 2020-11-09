//
//  SMR_data_p1.h
//  SMR_CPP
//
//  Created by Futao Zhang on 10/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__SMR_data_p1__
#define __SMR_CPP__SMR_data_p1__

#include "SMR_data.h"

namespace SMRDATA
{
    typedef struct{
        long lineNum;
        vector<string> _Expo_id;
        vector<int> _Expo_chr;
        vector<string> _Expo_gene;
        vector<int> _Expo_bp;
        vector<string> _Outco_id;
        vector<int> _Outco_chr;
        vector<string> _Outco_gene;
        vector<int> _Outco_bp;
        vector<string> _snp_rs;
        vector<int> _snp_chr;
        vector<int> _snp_bp;
        vector<string> _snp_a1;
        vector<string> _snp_a2;
        vector<float> _snp_frq;
        vector<double> _b;
        vector<double> _se;
        vector<double> _p_smr;
        vector<double> _p_heidi;
         vector<int> _nsnp;
        vector<uint32_t> _include;
    } eSMRrlt;
    


    
    void est_effect_splsize(char* eqtlsmaslstName, char* eqtlFileName, char* snplstName,char* problstName,char* snplst2exclde, char* problst2exclde,float thres);
    void make_cojo(char* outFileName,char* eqtlFileName, char* snplstName,char* snplst2exclde, char* problstName, char* problst2exclde, char* genelistName, bool bFlag);
   // void standardization(char* outFileName, char* eqtlFileName,bool bFlag,char* freqName, char* vpFileName);
    
    void lookup(char* outFileName,char* eqtlFileName, char* snplstName, char* problstName,char* genelistName, double plookup,bool bFlag, int chr, int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl,char* snpproblstName);
    
    void plot(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, char* geneAnnoName);
    
    void smr_multipleSNP(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, int expanWind, double ld_min,double threshpsmrest, bool sampleoverlap, double pmecs, int minCor, double ld_top_multi,double afthresh,double percenthresh);
    void smr_ivw_gwas_gwas(char* outFileName, char* bFileName,char* gwasFileName, char* gwasFileName2, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, double p_y, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, int expanWind, double ld_min,double threshpsmrest, bool sampleoverlap, double pmecs, int minCor, double ld_top_multi,double afthresh,double percenthresh);
    void smr_ivw_multipleSNP(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, double p_y, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, int expanWind, double ld_min,double threshpsmrest, bool sampleoverlap, double pmecs, int minCor, double ld_top_multi,double afthresh,double percenthresh);
    void smr_rev_multipleSNP(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_gwas, double p_y, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, int expanWind, double ld_min,double threshpsmrest, bool sampleoverlap, double pmecs, int minCor, double ld_top_multi,double afthresh,double percenthresh);
    void smr_e2e_mediation(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag, double ld_top,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool opt, double ld_min, bool cis2all, bool sampleoverlap, double pmecs, int minCor,double ld_top_multi, char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh, bool multi_corr_snp_flg, bool multi_exp_priority_flg, bool get_snp_effects_flg, int min_snp);
    void smr_e2multiple_mediation(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag, double ld_top,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool opt, double ld_min, bool cis2all, bool sampleoverlap, double pmecs, int minCor,double ld_top_multi, char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh, bool multi_corr_snp_flg, bool multi_exp_priority_flg, bool get_snp_effects_flg, int min_snp);
    void init_smr_wk(SMRWK* smrwk);
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl, bool heidioffFlag);
    void meta(char* besdlistFileName, char* outFileName, int meta_mth, double pthresh, bool cis_flag, int cis_itvl, int nmecs,char* problstName, char* problst2exclde, char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename, char* snplstName, char* snplst2exclde, int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb, bool smpwindFlag);
    void update_esifile(char* eqtlFileName,char* s_esiFileName);
    void update_epifile(char* eqtlFileName,char* s_epiFileName);
   
}

#endif /* defined(__SMR_CPP__SMR_data_p1__) */
