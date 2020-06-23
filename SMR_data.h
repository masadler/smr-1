//
//  SMR_data.h
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#ifndef __SRM_CPP__SMR_data__
#define __SRM_CPP__SMR_data__

#include <stdlib.h>

#include "bfile.hpp"

#define MAX_LEN_PRBNAME 40
#define MAX_LEN_GENENAME 40
#define MAX_LEN_PATH 512



extern int thread_num;
extern bool mute;
extern int xh; //using when testing
extern bool forcefrqck;
extern FILE* techeQTLfile;
extern char* outFileName;

namespace SMRDATA
{
   
    

    typedef struct{
        long snpNum;
        vector<string> snpName;
        vector<int> snpBp;
        vector<string> allele_1;
        vector<string> allele_2;
        vector<double> freq;
        vector<double> byz;
        vector<double> seyz;
        vector<double> pvalue;
        vector<uint32_t> splSize;
        vector<int> _include;
        map<string,int> _snp_name_map;
    } gwasData;
    
    typedef struct{
        vector<int> _esi_chr;
        vector<string> _esi_rs;
        vector<int> _esi_gd;
        vector<int> _esi_bp;
        vector<string> _esi_allele1;
        vector<string> _esi_allele2;
		vector<int> _esi_include; // initialized in the readesi
        map<string,int> _snp_name_map;
        vector<float> _esi_freq;
        
        vector<int> _epi_chr;
        vector<string> _epi_prbID;
        vector<int> _epi_gd;
        vector<int> _epi_bp;
        vector<string> _epi_gene;
        vector<char> _epi_orien;
        vector<int> _include; // initialized in the readepi
        map<string,int> _probe_name_map;
        vector<double> _epi_var;
        vector<int> _epi_start; // if no probe sequence region input, its size should be 0. for the probe not for probe sequence file, the value should be set as -9, no technical eQTL would be removed from this probe.
        vector<int> _epi_end;
        
        //for sparse
        vector<uint64_t> _cols;
        vector<uint32_t> _rowid;
        vector<float> _val;
        // for dense
        vector< vector<float> > _bxz; // first dimension is probe, second is snp
        vector< vector<float> > _sexz;
        
        uint64_t _probNum;
        uint64_t _snpNum;
        uint64_t _valNum;
        
    } eqtlInfo;
    
    typedef struct{
        char* probeId;
        char* genename;
        char* esdpath;
        char* bfilepath;
        int probechr;
        int gd;
        int bp;
        char orien;
        int start;
        int end;
    } probeinfolst;
    
    typedef struct{
        char* snprs;
        char* a1;
        char* a2;
        int snpchr;
        int gd;
        int bp;
        float beta;
        float se;
        float freq;
        float estn;
    } snpinfolst;

    typedef struct{
        vector<string> besdpath;        
        char* probeId;
        char* genename;
        int probechr;
        int gd;
        int bp;
        char orien;
        uint64_t vnum;
        uint32_t* rowid;
        float* beta_se;
        snpinfolst* sinfo;
        
    } probeinfolst2;
    
    typedef struct{
        char* snprs;
        char* a1;
        char* a2;
        int snpchr;
        int gd; // to store curID
        int bp;
        float beta;
        float se;
        float freq;
        float byz;
        float seyz;
        float pyz;
    } info4trans;
    
    typedef struct{
        int cur_chr;
        int cur_prbidx;
        vector<double> bxz, sexz,freq,byz,seyz;
        vector<double> pyz,zxz;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
    } SMRWK;

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
        string ProbeID;
        int ProbeChr;
        string Gene;
        int Probe_bp;
        string SNP;
        int SNP_Chr;
        int SNP_bp;
        string A1;
        string A2;
        float Freq;
        float b_GWAS;
        float se_GWAS;
        double p_GWAS;
        float b_eQTL;
        float se_eQTL;
        double p_eQTL;
        float b_SMR;
        float se_SMR;
        double p_SMR;
        double p_HET;
        int nsnp;
        char Orien;
        double p_SSMR;
    } SMRRLT;

    typedef struct{
        int* ptr;
        char* probeId;
        char* genename;
        char* esdpath;
        char* bfilepath;
        int probechr;
        int gd;
        int bp;
        char orien;
    } smr_probeinfo;
    
    typedef struct{
        int* rstr;
        bool* revs;
        char* snprs;
        char* a1;
        char* a2;
        int snpchr;
        int gd;
        int bp;
        float beta;
        float se;
        float freq;
        float estn;
    } smr_snpinfo;
    

    
    void read_bimfile(bInfo* bdata,string bimfile);
    void read_famfile(bInfo* bdata,string famfile);
    void read_bedfile(bInfo* bdata,string bedfile);
    void read_gwas_data(gwasData* gdata, char* gwasFileName);
    void read_esifile(eqtlInfo* eqtlinfo, string esifile, bool prtscr=true);
    void read_epifile(eqtlInfo* eqtlinfo, string epifile , bool prtscr=true);
    void read_besdfile(eqtlInfo* eqtlinfo, string besdfile,bool prtscr=true);
   
    bool has_prefix(const std::string &str, const std::string &prefix);
    bool has_suffix(const std::string &str, const std::string &suffix);
    void get_square_idxes(vector<int> &sn_ids,VectorXd &zsxz,double threshold);
	void est_cov_bxy(MatrixXd &covbxy, VectorXd &_zsxz, VectorXf &_bxy, VectorXd &_seyz, VectorXd &_bxz, MatrixXd &_LD_heidi);
	float bxy_hetero3(VectorXd &_byz, VectorXd &_bxz, VectorXd &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num);
    float bxy_mltheter_so(VectorXd &_byz, VectorXd &_bxz, VectorXd &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num, double theta);
    void allele_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata);
    void allele_check_opt(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata);
    void update_geIndx(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata);
    void allele_check(gwasData* gdata, eqtlInfo* esdata);
    void update_gwas(gwasData* gdata);
    
    bool make_XMat(bInfo* bdata, MatrixXd &X);
    void makex_eigenVector(bInfo* bdata,int j, VectorXd &x, bool resize, bool minus_2p);
    // inline functions
    template<typename ElemType>
    void makex(bInfo* bdata,int j, vector<ElemType> &x, bool minus_2p = false) {
        int i = 0;
        x.resize(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++) {
            if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) x[i] = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                else x[i] = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
            } else x[i] = bdata->_mu[bdata->_include[j]];
            if (minus_2p) x[i] -= bdata->_mu[bdata->_include[j]];
        }
    }
 
    void makeptrx(bInfo* bdata,int bsnpid,int cursnpid, float* x, bool minus_2p = false);
    void keep_indi(bInfo* bdata,string indi_list_file);
    void extract_snp(bInfo* bdata,string snplistfile);
    void extract_snp(bInfo* bdata, int chr);
    
    void extract_prob(eqtlInfo* eqtlinfo,string problstName);
	void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName);
    void exclude_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName);
    void exclude_prob(eqtlInfo* eqtlinfo,string problstName);
    
   
    
    void free_gwas_data(gwasData* gdata);  
  
    int file_read_check(ifstream* in_file, const char* filename);
    void init_smr_wk(SMRWK* smrwk);
    void init_smr_wk(SMRWKM* smrwk);
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP, int lowerbp,int upperbp,bool heidioffFlag);
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    void fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, gwasData* meddata, SMRWKM* smrwk,int cis_itvl);
    long fill_smr_wk_trans(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,bool heidioffFlag);
    double heidi_test(bInfo* bdata,SMRWK* smrwk, long maxid,double ld_top, double threshold, int m_hetero,long &nsnp );
    double heidi_test_new(bInfo* bdata,SMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, double theta);
    void update_snidx(SMRWK* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat);
    void extract_smrwk(SMRWK* smrwk,vector<int> &sn_ids,SMRWK* smrwk2);
    void rm_cor_sbat(MatrixXd &R, double R_cutoff, int m, vector<int> &rm_ID1);
    int maxabsid(vector<double> &zsxz, vector<int> &ids);
    void rm_cor_sbat(MatrixXd &R, double R_cutoff, int m, vector<int> &rm_ID1,vector<double> &zxz4smr);
    void update_smrwk_x(SMRWK* smrwk,vector<int> &sn_ids,MatrixXd &X);    
    void smr_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_het_mtd, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg);
    void sbat_calcu_lambda(MatrixXd &X, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx,vector<double> &zxz4smr, vector<double> &zyz4smr);
    int smr_setbased_test_ivw(bInfo* bdata, vector<uint32_t> &slctId, vector<double> &slct_bxz,vector<double> &slct_sexz,vector<double> &slct_byz,vector<double> &slct_seyz, double &beta_ivw, double &se_ivw, double &p_ivw, double p_smr, double ld_top_multi, vector<string> &snp4msmr);
    void smr(char* outFileName, char* bFileName,char* bldFileName, char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , int opt_hetero,char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest,bool new_het_mth, bool opt, char* prbseqregion, double ld_min, bool sampleoverlap, double pmecs, int minCor, char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh);
    
    void smr_trans(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,char* refSNP, bool heidioffFlag,int cis_itvl,int trans_itvl,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest, bool new_het_mth, double p_smr, double ld_min);
    void smr_trans_region(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,char* refSNP, bool heidioffFlag,int cis_itvl,int trans_itvl,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest, bool new_het_mth, double p_smr, bool opt, double ld_min,double afthresh,double percenthresh);
    
    void make_full_besd(char* outFileName, char* eqtlFileName, char* snplstName,char* problstName,bool bFlag,bool make_besd_flag, char* snplst2exclde, char* problst2exclde, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,int addn);
    
    void smr_e2e(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName,char* eproblstName,char* mproblstName,bool bFlag,double p_hetero,double ld_prune,int m_hetero, int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* eproblst2exclde,char* mproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm,double prepsmrthres, bool new_het_mth, bool opt, double ld_min, bool cis2all,bool sampleoverlap, double pmecs, int minCor,bool ssmrflg,int expanWind, double ld_top_multi,char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh);
    void smr_e2e_prbmatch(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName,char* eproblstName,char* mproblstName,bool bFlag,double p_hetero,double ld_prune,int m_hetero, int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* eproblst2exclde,char* mproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm,double prepsmrthres, bool new_het_mth, bool opt, double ld_min, bool cis2all,bool sampleoverlap, double pmecs, int minCor,int expanWind, double ld_top_multi,char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh);
    void ssmr_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, double ld_min,int opt_hetero, int expanWind,bool sampleoverlap, double pmecs, int minCor, double ld_top_multi);
    void meta(char* outFileName,char* eqtlFileName, char* eqtlFileName2);
    int comp_esi(const void *a,const void *b);
    int comp(const void *a,const void *b);
    int comp2(const void *a,const void *b);
    int comp_estn(const void *a,const void *b);
    int comp_i4tran(const void *a,const void *b);
    
    void read_smaslist(vector<string> &smasNames, string eqtlsmaslstName);
    void remove_indi(bInfo* bdata, string indi_list_file);
    void extract_snp(bInfo* bdata,string snplistfile);
    void exclude_snp(bInfo* bdata,string snplistfile);
    void allele_check(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata);
    void allele_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata);
    void update_geIndx(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata);
    void update_geIndx(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata);
    void progress_print(float progress);
    void make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool minus_2p = false);
    void get_square_ldpruning_idxes(vector<int> &sn_ids, VectorXd &zsxz, double threshold, VectorXd &ld_v, long maxid, double ld_top);
    void cor_calc(MatrixXd &LD, MatrixXd &X);
    void extract_prob_by_gene(eqtlInfo* eqtlinfo, string genelistName);
    void update_freq(char* eqtlFileName, char* frqfile);
    void write_besd(string outFileName, eqtlInfo* eqtlinfo);
    void extract_region_bp(bInfo* bdata, int chr, int fromkb, int tokb);
    void extract_eqtl_single_snp(eqtlInfo* eqtlinfo, string snprs);
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, string fromsnprs, string tosnprs);
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, int chr, int fromsnpkb, int tosnpkb);
    void extract_prob(eqtlInfo* eqtlinfo, string prbname, int prbWind);
    void extract_eqtl_single_probe(eqtlInfo* eqtlinfo, string prbname, bool prtscr=true);
    void extract_eqtl_prob(eqtlInfo* eqtlinfo, string fromprbname, string toprbname);
    void extract_eqtl_prob(eqtlInfo* eqtlinfo, int chr, int fromprbkb, int toprbkb);
    void extract_prob_by_single_gene(eqtlInfo* eqtlinfo, string genename);
    void extract_epi_by_chr(eqtlInfo* eqtlinfo, int prbchr);
    void extract_eqtl_by_chr(eqtlInfo* eqtlinfo, int snpchr);
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snporprb, int Wind, string msg);
    void epi_man(eqtlInfo* eqtlinfo,char* problstName,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename);
    void esi_man(eqtlInfo* eqtlinfo,char* snplstName,int chr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, int cis_itvl,const char* prbname);
    void exclude_eqtl_single_probe(eqtlInfo* eqtlinfo, string prbname);
    void slct_sparse_per_prb(vector<int> &slct_idx, probeinfolst* prbifo, vector<snpinfolst> &snpinfo, long cis_itvl, long trans_itvl,double transThres,double restThres,FILE* logfile, bool extract_cis_only, bool techHit=false);
    void read_gene_anno(char* geneAnnoName,vector<int> &chr, vector<string> &genename,vector<int> &start,vector<int> &end);
    void read_gene_anno_strand(char* geneAnnoName,vector<int> &chr, vector<string> &genename,vector<int> &start,vector<int> &end,vector<string> &strand);
    void rm_unmatched_snp(gwasData* gdata, eqtlInfo* esdata);
    void rm_unmatched_snp(eqtlInfo* etrait, eqtlInfo* esdata);
    void filter_snp_null(eqtlInfo* eqtlinfo);
    int max_zsmr_id(SMRWK *smrwk , double p_smr);
    double heidi_test_ref_new(bInfo* bdata,SMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp, int refid, double ld_min,int opt_hetero, bool sampleoverlap, double theta);
    double heidi_test_ref_new(ldInfo* ldinfo,FILE* ldfptr, SMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp, int refid , double ld_min, int opt_hetero , bool sampleoverlap, double theta);
    void free_snplist(vector<snpinfolst> &a);
    void free_probelist(vector<probeinfolst> &a);
    void free_snplist(vector<info4trans> &a);
    void free_probelist(vector<probeinfolst2> &a);
    void read_epistartend(eqtlInfo* eqtlinfo,char* prbseqregion);
    int shown(string besdfile);
    double freq_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata, double &freqthresh, double &percenthresh);
    double freq_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata, double &freqthresh, double &percenthresh);
    void slct_trans_per_prb(vector<int> &slct_idx, vector<int> &regionChr, vector<long> &snpNumPerRegion,vector<long> &leftbound, vector<long> &rightbound, probeinfolst* prbifo, vector<info4trans> &snpinfo, long cis_itvl, long trans_itvl,double transThres);
    void extract_targets(eqtlInfo* eqtlinfo, string snpprblistfile, map<string, string> &prb_snp);
}
#endif /* defined(__SRM_CPP__SMR_data__) */
