/*
 * SMR_data_p4.cpp
 * Implementations of univariable and multivariable MR-IVW implementations
 *
 * 2020 by Marie Sadler
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  
 */

#ifndef __SMR_CPP__SMR_data_p4__
#define __SMR_CPP__SMR_data_p4__

#include "SMR_data.h"
#include "SMR_data_p1.h"
#include <numeric>
#include <Eigen_unsupported/Eigen/KroneckerProduct>

namespace SMRDATA
{
    typedef struct{
        int cur_chr;
        int cur_prbidx;
        vector<double> bxz, sexz,freq,byz,seyz,bmz,semz;
        vector<double> pyz,zxz,zmz;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
    } SMRWKM;

    typedef struct{
        vector<string> Med_ID, Gene_ID;
        vector<int> Med_Chr, Med_bp;
        vector<double> b_xm, se_xm, p_xm, n_xm, b_my, se_my, p_my, n_my;
    } MED_info;

    typedef struct{
        int cur_chr;
        int cur_prbidx;
        int medNum_cis;
        vector<double> bxz, sexz,freq,byz,seyz;
        vector<double> pyz,zxz;
        MatrixXd Xm_b, Xm_se;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
    } SMRWKMULT;

    typedef struct{
        int cur_chr;
        int cur_prbidx;
        int medNum_trans;
        bool medFocal;
        bool medFocal_core;
        vector<double> bxz, sexz, freq, byz, seyz, bmz, semz;
        vector<double> pyz, zxz, zmz;
        MatrixXd Xm_b, Xm_se, Xm_z;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs, allele1, allele2;
    } SMRWKMULTTRANS;

    void lookup_dense(char* outFileName,char* eqtlFileName, char* snplstName, char* problstName,char* genelistName, double plookup, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl,char* snpproblstName);
    
    void allele_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata);
    void allele_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata, bool incmp_med);
    void allele_check(bInfo* bdata, gwasData* gdata, gwasData* gdata2);
    double freq_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata, double &freqthresh, double &percenthresh);
    double freq_check(bInfo* bdata, gwasData* gdata, gwasData* gdata2, double &freqthresh, double &percenthresh);
    void update_geIndx(bInfo* bdata, gwasData* gdata, gwasData* gdata2);
    void update_geIndx(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata);
    void init_smr_wk(SMRWKM* smrwk);
    void init_med_info(MED_info* med_info);
    void init_smr_wk(SMRWKMULT* smrwk);
    void init_smr_wk(SMRWKMULTTRANS* smrwk);
    void fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, gwasData* meddata, SMRWKM* smrwk,int cis_itvl);
    long fill_smr_wk_trans(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk);
    void fill_smr_wk_trans(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, gwasData* meddata, SMRWKM* smrwk);

    void sbat_calcu_lambda(MatrixXd &X, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx,vector<double> &zxz4smr);
    vector<int> sort_indexes(vector<double> &v);
    double qnorm5(double p, int lower_tail = 0);
    int min_id(vector<double> &pyz);
    int min_id(VectorXd &pyz);
    
    void mr_ivw(VectorXd &bzx, VectorXd &bzy, VectorXd &sezx, VectorXd &sezy, double &beta_ivw, double &se_ivw, double &p_ivw);
    void mr_ivw_LD(VectorXd &bzx, VectorXd &bzy, VectorXd &sezx, VectorXd &sezy, MatrixXd &C, double &beta_ivw, double &se_ivw, double &p_ivw);
    void mvmr_ivw_LD(MatrixXd &X, VectorXd &y, VectorXd &SEs, MatrixXd &C, VectorXd &betas, VectorXd &vars);
    void mr_mediation_varsub(MatrixXd &X, VectorXd &y, MatrixXd &X_se, VectorXd &betas, VectorXd &vars);
    
    int smr_ivw_test(bInfo* bdata, vector<uint32_t> &slctId, vector<string> &slct_snpName, vector<string> &slct_a1, vector<string> &slct_a2, vector<double> &slct_bxz,vector<double> &slct_sexz, vector<double> &slct_pxz, vector<double> &slct_byz,vector<double> &slct_seyz, double &beta_ivw, double &se_ivw, double &p_ivw, double &h2cis, double p_smr, double ld_top_multi, bool ldmatrix, int N_ref, double trev, vector<string> &snp4msmr, bool p_expo_provided, bool get_snp_effects_flg, int min_snp, FILE* snpfile, string snpfilename, string exponame);
    void smr_e2e_prbmatch(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName,char* eproblstName,char* mproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* eproblst2exclde,char* mproblst2exclde,double p_smr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp);
    void smr_e2e_cis(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool multi_corr_snp_flg, bool get_snp_effects_flg, int min_snp);
    void smr_ivw_gwas_gwas(char* outFileName, char* bFileName,char* gwasFileName, char* gwasFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp);
    void smr_ivw_analysis(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp);
    void smr_ivw_trans_analysis(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp);
    void smr_rev_ivw_analysis(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_gwas, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp);
    void smr_mediation_prbmatch(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool multi_corr_snp_flg, bool get_snp_effects_flg, int min_snp);
    void smr_multi_mediation_cis(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool multi_corr_snp_flg, bool get_snp_effects_flg, bool get_snp_effects_top_flg, int min_snp, bool incmp_expo, double med_R_thresh, bool uncorr_med, bool exclude_top_SNP, double p_shrinkage, double p_expo_med);
    void smr_multi_mediation_trans(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp, bool incmp_med, char* core_med_lst, char* qtl_med_lst, char* med_corr_snps_lst, int num_med, double med_R_thresh, bool uncorr_med, double p_shrinkage, double p_expo_med);

    void calcul_cormat(char* outFileName, char* eqtlFileName, char* eqtlFileName2, char* snplstName, char* problstName, char* genelistName, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl);
    void calcul_cormat(char* outFileName, char* eqtlFileName, char* snplstName, char* problstName, char* genelistName, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl);

}

#endif /* defined(__SMR_CPP__SMR_data_p4__) */