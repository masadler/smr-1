//
//  SMR.cpp
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//


#include "StrFunc.h"
#include "SMR_data_p1.h"
#include "SMR_data_p2.h"
#include "SMR_data_p3.h"
#include "SMR_data_p4.h"
#include "SMR_plot.h"
#include "bfile.hpp"
#include "SMR.h"

using namespace SMRDATA;
using namespace StatFunc;
int thread_num;
bool mute = false;
int xh=0;
bool forcefrqck = true;
char* outFileName=NULL;


int main(int argc, char** argv) {
 
    cout << "*******************************************************************" << endl;
    cout << "* SMR-IVW: Extension of Summary-data-based Mendelian Randomization (SMR)" << endl;
    cout << "* version 1.0"<< endl;
    cout << "* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* Extensions implemented by Marie Sadler, University of Lausanne"<< endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;
   
    FLAGS_VALID_CK(argc,argv);
    long int time_used = 0, start = time(NULL);
    string months[] = {"Jan", "Feb", "Mar", "Apr", "May",
        "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    string weeks[] = {"Sun","Mon", "Tue", "Wed", "Thu", "Fri",
        "Sat" };
    time_t t = time(0);
    tm* now = localtime(&t);
    cout << "Analysis started: "<<now->tm_hour<<":"<<now->tm_min<<":"<< now->tm_sec<<","<<weeks[now->tm_wday]<<" "<<months[(now->tm_mon)]<<" "<<now->tm_mday<<","<<(now->tm_year + 1900)<<endl;
    
    cout << "\nOptions:" << endl;
    try {
        option(argc, argv);
    } catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    } catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }    
    t = time(0);
    now = localtime(&t);
    cout << "\nAnalysis completed: "<<now->tm_hour<<":"<<now->tm_min<<":"<< now->tm_sec<<","<<weeks[now->tm_wday]<<" "<<months[(now->tm_mon)]<<" "<<now->tm_mday<<","<<(now->tm_year + 1900)<<endl;
    time_used = time(NULL) - start;
    cout << "Computational time: " << time_used / 3600 << ":" << (time_used % 3600) / 60 << ":" << time_used % 60 << endl;
    
    return 0;
}

void option(int option_num, char* option_str[])
{
    char* bFileName=NULL;
    char* gwasFileName=NULL;
    char* eqtlFileName=NULL;
    char* indilstName=NULL;
    char* snplstName=NULL;
    char* indilst2remove=NULL;
    char* snplst2exclde=NULL;
    char* problst2exclde=NULL;
    char* eproblst2exclde=NULL;
     char* oproblst2exclde=NULL;
    bool bFlag=false;// for binary file
    double maf=0.0;
    double p_hetero=1.5654e-3;
    double p_smr=5.0e-8;
    double p_medsmr=5.0e-8;
    double p_gwas=5.0e-8;
    double ld_prune=0.9;
    double ld_min=0.05;
    double ld_prune_multi=0.1;
    unsigned int m_hetero=3;
    int opt_hetero=20;
    bool smr_flag=true;
    bool smr_trans_flag=false;
    
    // for data management
    bool make_besd_flag=false;
    char* problstName=NULL;
    char* eproblstName=NULL;
    char* oproblstName=NULL;
    char* targetsnpproblstName=NULL;
    char* snpproblstName=NULL;
    
    char* eFileName=NULL;
    bool eremlFlag=false;
    double pexsnp = -9;
    double pinsnp = -9;
    
    // for filtering
    bool cis_flag=false;
    int cis_itvl=2000;
    int trans_itvl=1000;
    double transThres=5.0e-8;
    double cisThres=5.0e-8;
    double restThres=1.0e-5;
    
    // for lookup
    double plookup=5e-8;
    bool lookup_flag=false;
    char* genelistName=NULL;
    int chr=0;
    int prbchr=0;
    int snpchr=0;
    char* snprs=NULL;
    char* snp2rm = NULL;
    char* prbname=NULL;
    char* fromsnprs=NULL;
    char* tosnprs=NULL;
    char* fromprbname=NULL;
    char* toprbname=NULL;
    int snpWind=50;   //Kb
    int prbWind=1000; //Kb
    char* genename=NULL;
    int fromsnpkb=-9;
    int tosnpkb=-9;
    int fromprbkb=-9;
    int toprbkb=-9;
    bool snpwindFlag=false;
    bool prbwindFlag=false;
    bool trans_dense=false;
    
    
    char* refSNP=NULL;
    bool heidioffFlag = false;
    
    thread_num = 1;
    
    bool combineFlg = false;
    char* eqtlsmaslstName = NULL;
    //
    char* gwasFileName2=NULL;
    char* eqtlFileName2=NULL; // for the outcome trait
    char* traitlstName=NULL;
    //
    bool plotflg=false;
    //
    char* syllabusName=NULL;
    bool gctaflag=false;
    bool plinkflag=false;
    bool gemmaflag=false;
    bool merlinflag=false;
    bool boltflag=false;
    bool mateqtlflag = false;
    bool fastqtlnflag = false;
    bool fastqtlpflag = false;
    bool qtltoolsnflag = false;
    bool qtltoolspflag = false;
    bool genouni=false;
    
    
    char* freqName=NULL;
    bool esdstd=false;
    char* vpFileName=NULL;
    
    bool metaflg=false;
    bool est_effe_spl_size_flg=false;
    // for SMR e2me
    int outcomePrbWind=2000; //kb
    char* eprobe=NULL;
    char* oprobe=NULL;
    char* eprobe2rm=NULL;
    char* oprobe2rm=NULL;
    bool prbmatch=false;
    bool e2e_cis=false;
    bool mediation=false;
    bool mediation_multiple_cis=false;
    bool mediation_multiple_trans=false;

    // for MR-IVW Steiger filtering
    double trev = -2;

    // for mediation analysis
    bool multi_corr_snp_flg = true;
    bool get_snp_effects_flg = false;
    int min_snp = 5;
    char* core_med_lst = NULL;
    char* qtl_med_lst = NULL;
    char* med_corr_snps_lst = NULL;
    bool incmp_expo = false;
    int num_med = -9;
    double med_R_thresh=-9;
    bool uncorr_med = false;
    double p_shrinkage=1;
    double p_expo_med = 0.01;
    bool calcul_cormat_flag = false;

    // for internal test
    char* smrRltFileName = NULL;    
    bool recodeflg=false;
    // for setBased SMR
    char* geneAnnoName=NULL;
    char* setlstName=NULL;
    int setWind=-9;
    bool ssmrflg=false;
    bool ssmrivwflg=false;
    bool ssmrrevflg=false;
    bool ssmrgwastogwasflg=false;
    bool ldmatrix=false;
    //
    char* queryFileName=NULL;
    bool queryfileflg=false;
    
    bool save_dense_flag=false;
    
    double pthres_me2esmr=1.0;
    
    double threshpsmrest=1.0;
    double threshphet=0.0;
    bool qcflag=false;
    int qcmtd=0;
    int z_thresh=3;
    bool new_het_mth=true;
    
    bool count_cis_flag = false;
    bool count_trans_flag = false;
    
    // for pick up instrument.// inflation, DO NOT USE
    bool opt_slct_flag = false;
    
    int meta_mtd = 0; //0 for meta, 1 for mecs
    double pmecs = 0.01;
    int nmecs = 100; // the number of common SNPs to calculation the correlation in MeCS
    
    bool extract_cis_only = false;
    bool rm_technical = false;
    char* prbseqregion = NULL;
    double ptech = 5.0e-8;
    
    int addn=-9;
    bool shownflag=false;
    
    char* refepiName = NULL;
    char* refesiName = NULL;
    
    //for molecular trait
    bool cis2all = false;
    
     bool sampleoverlap = false;
    int minsnpcor = 2;
    
    double diff_freq = 0.2;
    double diff_freq_ratio = 0.05;
    
    int ldWind = 4000;
    bool ldr = false;
    bool ldr2 = false;
    
    char* bldFileName = NULL;
    bool make_bld_flag = false;
    for(int i=0;i<option_num;i++)
    {
        // Plink files for LD from a reference sample
        if(0==strcmp(option_str[i],"--bfile")){
            bFileName=option_str[++i];
            FLAG_VALID_CK("--bfile", bFileName);
            printf("--bfile %s\n",bFileName);
        }
        else if(0==strcmp(option_str[i],"--bld")){
            bldFileName=option_str[++i];
            FLAG_VALID_CK("--bld", bldFileName);
            printf("--bld %s\n",bldFileName);
        }
        else if(strcmp(option_str[i],"--make-bld")==0){
            make_bld_flag=true;
            cout<<"--make-bld "<<endl;
        }
        // gwas data file as cojo format
        else if(0==strcmp(option_str[i],"--gwas-summary")){
            if(gwasFileName ==NULL){
                gwasFileName=option_str[++i];
                FLAG_VALID_CK("--gwas-summary", gwasFileName);
                printf("--gwas-summary %s\n",gwasFileName);
                CommFunc::FileExist(gwasFileName);
            }else
            {
                gwasFileName2=option_str[++i];
                FLAG_VALID_CK("--gwas-summary", gwasFileName2);
                printf("--gwas-summary %s\n",gwasFileName2);
                CommFunc::FileExist(gwasFileName2);
            }
        }
        // eQTL files
        else if(0==strcmp(option_str[i],"--eqtl-summary")){
            if(eqtlFileName==NULL)
            {
                eqtlFileName=option_str[++i];
                FLAG_VALID_CK("--eqtl-summary", eqtlFileName);
                printf("--eqtl-summary %s\n",eqtlFileName);
            }else{
                eqtlFileName2=option_str[++i];
                FLAG_VALID_CK("--eqtl-summary", eqtlFileName2);
                printf("--eqtl-summary %s\n",eqtlFileName2);
            }
        }
        else if(0==strcmp(option_str[i],"--beqtl-summary")){
             bFlag=true;
            if(eqtlFileName==NULL)
            {
                eqtlFileName=option_str[++i];
                FLAG_VALID_CK("--beqtl-summary", eqtlFileName);
                printf("--beqtl-summary %s\n",eqtlFileName);
            }else{
                eqtlFileName2=option_str[++i];
                FLAG_VALID_CK("--beqtl-summary", eqtlFileName2);
                printf("--beqtl-summary %s\n",eqtlFileName2);
            }
        }
        else if(strcmp(option_str[i],"--keep")==0){
            indilstName=option_str[++i];
            FLAG_VALID_CK("--keep", indilstName);
            cout<<"--keep "<<indilstName<<endl;
            CommFunc::FileExist(indilstName);
        }
        else if(strcmp(option_str[i],"--remove")==0){
            indilst2remove=option_str[++i];
            FLAG_VALID_CK("--remove", indilst2remove);
            cout<<"--remove "<<indilst2remove<<endl;
            CommFunc::FileExist(indilst2remove);
        }
        else if(strcmp(option_str[i],"--extract-snp")==0){
            snplstName=option_str[++i];
            FLAG_VALID_CK("--extract-snp", snplstName);
            cout<<"--extract-snp "<<snplstName<<endl;
            CommFunc::FileExist(snplstName);
        }
        else if(strcmp(option_str[i],"--extract-snp-p")==0){
            pinsnp=atof(option_str[++i]);
            cout<<"--extract-snp-p "<<pinsnp<<endl;
            if(pinsnp<0 || pinsnp>1)
            {
                fprintf (stderr, "Error: --extract-snp-p should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--exclude-snp-p")==0){
            pexsnp=atof(option_str[++i]);
            cout<<"--exclude-snp-p "<<pexsnp<<endl;
            if(pexsnp<0 || pexsnp>1)
            {
                fprintf (stderr, "Error: --exclude-snp-p should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }

        else if(strcmp(option_str[i],"--extract-probe")==0){
            problstName=option_str[++i];
            FLAG_VALID_CK("--extract-probe", problstName);
            cout<<"--extract-probe "<<problstName<<endl;
            CommFunc::FileExist(problstName);
        }
        else if(strcmp(option_str[i],"--extract-exposure-probe")==0){
            eproblstName=option_str[++i];
            FLAG_VALID_CK("--extract-exposure-probe", eproblstName);
            cout<<"--extract-exposure-probe "<<eproblstName<<endl;
            CommFunc::FileExist(eproblstName);
        }
        else if(strcmp(option_str[i],"--extract-outcome-probe")==0){
            oproblstName=option_str[++i];
            FLAG_VALID_CK("--extract-outcome-probe", oproblstName);
            cout<<"--extract-outcome-probe "<<oproblstName<<endl;
            CommFunc::FileExist(oproblstName);
        }
        else if(strcmp(option_str[i],"--exclude-snp")==0){
            snplst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-snp", snplst2exclde);
            cout<<"--exclude-snp "<<snplst2exclde<<endl;
            CommFunc::FileExist(snplst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-probe")==0){
            problst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-probe", problst2exclde);
            cout<<"--exclude-probe "<<problst2exclde<<endl;
            CommFunc::FileExist(problst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-exposure-probe")==0){
            eproblst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-exposure-probe", eproblst2exclde);
            cout<<"--exclude-exposure-probe "<<eproblst2exclde<<endl;
            CommFunc::FileExist(eproblst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-outcome-probe")==0){
            oproblst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-outcome-probe", oproblst2exclde);
            cout<<"--exclude-outcome-probe "<<oproblst2exclde<<endl;
            CommFunc::FileExist(oproblst2exclde);
        }

        else if(strcmp(option_str[i],"--maf")==0){
            maf=atof(option_str[++i]);
            cout<<"--maf "<<maf<<endl;
            if(maf<0 || maf>0.5)
            {
                fprintf (stderr, "Error: --maf should be within the range from 0 to 0.5.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--make-besd")==0){
            make_besd_flag=true;
            cout<<"--make-besd "<<endl;
        }
        else if (0 == strcmp(option_str[i], "--out")){
            outFileName = option_str[++i];
            if(outFileName !=NULL && has_prefix(outFileName, "--"))
            {
                outFileName=NULL;
                i--;
            }else   printf("--out %s\n", outFileName);
        }
        else if (0 == strcmp(option_str[i], "--peqtl-smr")){
            p_smr = atof(option_str[++i]);
            printf("--peqtl-smr %10.2e\n", p_smr);
            if(p_smr<0 || p_smr>1)
            {
                fprintf (stderr, "Error: --peqtl-smr should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--pmqtl-smr")){
            p_medsmr = atof(option_str[++i]);
            printf("--pmqtl-smr %10.2e\n", p_medsmr);
            if(p_medsmr<0 || p_medsmr>1)
            {
                fprintf (stderr, "Error: --pmqtl-smr should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--pgwas-smr")){
            p_gwas = atof(option_str[++i]);
            printf("--pgwas-smr %10.2e\n", p_gwas);
            if(p_gwas<0 || p_gwas>1)
            {
                fprintf (stderr, "Error: --pgwas-smr should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--peqtl-heidi")){
            p_hetero = atof(option_str[++i]);
            printf("--peqtl-heidi %10.2e\n", p_hetero);
            if(p_hetero<0 || p_hetero>1)
            {
                fprintf (stderr, "Error: --peqtl-heidi should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
            
        }
        else if (0 == strcmp(option_str[i], "--heidi-min-m")){
            m_hetero = atoi(option_str[++i]);
            printf("--heidi-min-m %d\n", m_hetero);
        }
        else if (0 == strcmp(option_str[i], "--heidi-max-m")){
            opt_hetero = atoi(option_str[++i]);
            printf("--heidi-max-m %d\n", opt_hetero);
        }
        else if (0 == strcmp(option_str[i], "--heidi-mtd")){
            int mtd=atoi(option_str[++i]);
            if(mtd==0) new_het_mth=false;
            else if(mtd==1) new_het_mth=true;
            else {
                printf("ERROR: only tow HEIDI methold avaliable.\n");
                exit(EXIT_FAILURE);
            }
            printf("--heidi-mtd %d\n", mtd);
        }
        else if (0 == strcmp(option_str[i], "--smr-multi")){
            ssmrflg=true;
            printf("--smr-multi \n");
        }
        else if (0 == strcmp(option_str[i], "--smr-ivw")){
            ssmrivwflg=true;
            printf("--smr-ivw \n");
        }
        else if (0 == strcmp(option_str[i], "--smr-gwas-to-gwas")){
            ssmrgwastogwasflg=true;
            printf("--smr-gwas-to-gwas \n");
        }
        else if (0 == strcmp(option_str[i], "--smr-multi-rev")){
            ssmrrevflg=true;
            printf("--smr-multi-rev \n");
        }
        else if (0 == strcmp(option_str[i], "--ld-matrix")){
            ldmatrix=true;
            printf("--ld-matrix \n");
        }
        else if (0 == strcmp(option_str[i], "--ld-upper-limit")){
            ld_prune = atof(option_str[++i]);
            printf("--ld-upper-limit %f\n", ld_prune);
            if(ld_prune<0 || ld_prune>1)
            {
                fprintf (stderr, "Error: --ld-upper-limit should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--ld-multi-snp")){
            ld_prune_multi = atof(option_str[++i]);
            printf("--ld-multi-snp %f\n", ld_prune_multi);
            if(ld_prune_multi<0 || ld_prune_multi>1)
            {
                fprintf (stderr, "Error: --ld-multi-snp should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--smr")){
            smr_flag=true;
            printf("--smr \n");
        }
        else if (0 == strcmp(option_str[i], "--make-besd-dense")){
            make_besd_flag=true;
            save_dense_flag=true;
            printf("--make-besd-dense \n");
        }
        else if(strcmp(option_str[i],"--cis-wind")==0){
            cis_flag=true;
            cis_itvl=atoi(option_str[++i]);
            printf("--cis-wind %d Kb\n", cis_itvl);
            if(cis_itvl<=0 )
            {
                fprintf (stderr, "Error: --cis-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--trans-wind")==0){
            trans_itvl=atoi(option_str[++i]);
            printf("--trans-wind %d Kb\n", trans_itvl);
            if(trans_itvl<0 )
            {
                fprintf (stderr, "Error: --trans-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--peqtl-trans")==0){
            transThres=atof(option_str[++i]);
            printf("--peqtl-trans %10.2e\n", transThres);
            if(transThres<0 || transThres>1)
            {
                fprintf (stderr, "Error: --peqtl-trans should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--peqtl-cis")==0){
            cisThres=atof(option_str[++i]);
            printf("--peqtl-cis %10.2e\n", cisThres);
            if(cisThres<0 || cisThres>1)
            {
                fprintf (stderr, "Error: --peqtl-cis should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--peqtl-other")==0){
            restThres=atof(option_str[++i]);
            printf("--peqtl-other %10.2e\n", restThres);
            if(restThres<0 || restThres>1)
            {
                fprintf (stderr, "Error: --peqtl-other should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--efile")){
            eremlFlag=true;
            eFileName = option_str[++i];
            FLAG_VALID_CK("--efile", eFileName);
            printf("--efile %s\n", eFileName);
        }
        else if (0 == strcmp(option_str[i], "--query")){
            lookup_flag=true;
            if(i+1==option_num || has_prefix(option_str[i+1],"--"))  plookup = 5.0e-8;
            else plookup = atof (option_str[++i]);
            if(bldFileName) printf("--query\n");
            else printf("--query %10.2e\n", plookup);
            if(plookup<0 || plookup>1)
            {
                fprintf (stderr, "Error: --query should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(0==strcmp(option_str[i],"--trans-dense")){
            trans_dense=true;
            printf("--trans-dense\n");
        }
        else if(strcmp(option_str[i],"--genes")==0){
            genelistName=option_str[++i];
            FLAG_VALID_CK("--genes", genelistName);
            cout<<"--genes "<<genelistName<<endl;
            CommFunc::FileExist(genelistName);
        }
        else if(strcmp(option_str[i],"--gene-list")==0){
            geneAnnoName=option_str[++i];
            FLAG_VALID_CK("--gene-list", geneAnnoName);
            cout<<"--gene-list "<<geneAnnoName<<endl;
            CommFunc::FileExist(geneAnnoName);
        }
        else if(strcmp(option_str[i],"--set-list")==0){
            setlstName=option_str[++i];
            FLAG_VALID_CK("--set-list", setlstName);
            cout<<"--set-list "<<setlstName<<endl;
            CommFunc::FileExist(setlstName);
        }
        else if (0 == strcmp(option_str[i], "--target-snp")){
            refSNP = option_str[++i];
            FLAG_VALID_CK("--target-snp", refSNP);
            printf("--target-snp %s\n", refSNP);
        }
        else if (0 == strcmp(option_str[i], "--heidi-off")){
            heidioffFlag = true;            
            printf("--heidi-off \n");
        }
        else if(0 == strcmp(option_str[i],"--thread-num")){
            thread_num=atoi(option_str[++i]);
            printf("--thread-num %d\n", thread_num);
        }       
        else if (0 == strcmp(option_str[i], "--besd-flist")){
            combineFlg=true;
            eqtlsmaslstName = option_str[++i];
            FLAG_VALID_CK("--besd-flist", eqtlsmaslstName);
            printf("--besd-flist %s\n", eqtlsmaslstName);
            CommFunc::FileExist(eqtlsmaslstName);
        }
        else if (0 == strcmp(option_str[i], "--plot")){
            plotflg = true;
            printf("--plot \n" );
        }
        else if (0 == strcmp(option_str[i], "--trans")){
            smr_trans_flag = true;
            if(!cis_flag) cis_itvl=5000;
            printf("--trans \n" );
        }
        else if (0 == strcmp(option_str[i], "--eqtl-flist")){
            syllabusName = option_str[++i];
             gctaflag = true;
            if(syllabusName !=NULL && has_prefix(syllabusName, "--"))
            {
                syllabusName=NULL;
                i--;
            }else   printf("--eqtl-flist %s\n", syllabusName);
        }
        else if (0 == strcmp(option_str[i], "--smr-format")){
            gctaflag = true;
            printf("--smr-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--plink-qassoc-format")){
            plinkflag = true;
             gctaflag = false;
            printf("--plink-qassoc-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--gemma-format")){
            gemmaflag = true;
            gctaflag = false;
            printf("--gemma-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--bolt-assoc-format")){
            boltflag = true;
            gctaflag = false;
            printf("--bolt-assoc-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--geno-uni")){
            genouni = true;
            printf("--geno-uni \n" );
        }
        else if (0 == strcmp(option_str[i], "--merlin-fastassoc-format")){
            merlinflag = true;
            gctaflag = false;
            printf("--merlin-fastassoc-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--matrix-eqtl-format")){
            mateqtlflag = true;
            gctaflag = false;
            printf("--matrix-eqtl-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--fastqtl-nominal-format")){
            fastqtlnflag = true;
            gctaflag = false;
            printf("--fastqtl-nominal-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--fastqtl-permu-format")){
            fastqtlpflag = true;
            gctaflag = false;
            printf("--fastqtl-permu-format \n" );
        }
        else if(strcmp(option_str[i],"--extract-trait")==0){
            traitlstName=option_str[++i];
            FLAG_VALID_CK("--extract-trait", traitlstName);
            cout<<"--extract-trait "<<traitlstName<<endl;
            CommFunc::FileExist(traitlstName);
        }
        else if(strcmp(option_str[i],"--esd-std")==0){
            esdstd=true;
            cout<<"--esd-std "<<endl;
        }
        else if(strcmp(option_str[i],"--update-freq")==0){
            freqName=option_str[++i];
            FLAG_VALID_CK("--update-freq", freqName);
            cout<<"--update-freq "<<freqName<<endl;
            CommFunc::FileExist(freqName);
        }
        else if(strcmp(option_str[i],"--probe-var")==0){
            vpFileName=option_str[++i];
            FLAG_VALID_CK("--probe-var", vpFileName);
            cout<<"--probe-var "<<vpFileName<<endl;
            CommFunc::FileExist(vpFileName);
        }
        else if(strcmp(option_str[i],"--meta")==0){
            metaflg=true;
            combineFlg=false;
            smr_flag=false;
            cout<<"--meta "<<endl;
        }
        else if(0==strcmp(option_str[i],"--mecs")){
            metaflg=true;
            meta_mtd=1;
            printf("--mecs\n");
        }
        else if(0==strcmp(option_str[i],"--pmecs")){
            pmecs=atof(option_str[++i]);
            if (pmecs < 0 || pmecs > 1) {
                printf("\nError: --pmecs should be between 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
            printf("--pmecs %0.2e\n", pmecs);
        }
        else if(0==strcmp(option_str[i],"--nmecs")){
            nmecs=atoi(option_str[++i]);
            if (nmecs <= 0 ) {
                printf("\nError: --nmecs should be over 0.\n");
                exit (EXIT_FAILURE);
            }
            printf("--nmecs %d\n", nmecs);
        }
        else if(strcmp(option_str[i],"--smr-wind")==0){
            if(smr_trans_flag)
            {
                trans_itvl=atoi(option_str[++i]);
                printf("--smr-wind %d Kb\n", trans_itvl);
                if(trans_itvl<0 )
                {
                    fprintf (stderr, "Error: --smr-wind should be over 0.\n");
                    exit (EXIT_FAILURE);
                }
            }
            else
            {
                cis_itvl=atoi(option_str[++i]);
                printf("--smr-wind %d Kb\n", cis_itvl);
                if(cis_itvl<0 )
                {
                    fprintf (stderr, "Error: --smr-wind should be over 0.\n");
                    exit (EXIT_FAILURE);
                }
            }
            
        }
        else if(strcmp(option_str[i],"--outcome-wind")==0){
            outcomePrbWind=atoi(option_str[++i]);
            printf("--outcome-wind %d Kb\n", outcomePrbWind);
            if(outcomePrbWind<0 )
            {
                fprintf (stderr, "Error: --outcome-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--set-wind")==0){
            setWind=atoi(option_str[++i]);
            printf("--set-wind %d kb\n", setWind);
            if(setWind<=0 )
            {
                fprintf (stderr, "Error: --set-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--est-n")==0){
            est_effe_spl_size_flg=true;
            combineFlg=false;
            cout<<"--est-n "<<endl;
        }
        else if(strcmp(option_str[i],"--smr-file")==0){
            smrRltFileName=option_str[++i];
            FLAG_VALID_CK("--smr-file", smrRltFileName);
            cout<<"--smr-file "<<smrRltFileName<<endl;
            CommFunc::FileExist(smrRltFileName);
        }
        else if(strcmp(option_str[i],"--recode")==0){
            recodeflg=true;
            cout<<"--recode "<<endl;
        }
        else if(strcmp(option_str[i],"--chr")==0){
            char* tmpstr=option_str[++i];
            if(strncmp(tmpstr,"X",1)==0) chr=23;
            else if(strncmp(tmpstr,"Y",1)==0) chr=24;
            else chr=atoi(tmpstr);
            FLAG_VALID_CK("--chr", tmpstr);
            cout<<"--chr "<<tmpstr<<endl;
        }
        else if(strcmp(option_str[i],"--probe-chr")==0){
            char* tmpstr=option_str[++i];
            if(strncmp(tmpstr,"X",1)==0) prbchr=23;
            else if(strncmp(tmpstr,"Y",1)==0) prbchr=24;
            else prbchr=atoi(tmpstr);
            FLAG_VALID_CK("--probe-chr", tmpstr);
            cout<<"--probe-chr "<<tmpstr<<endl;
        }
        else if(strcmp(option_str[i],"--snp-chr")==0){
            char* tmpstr=option_str[++i];
            if(strncmp(tmpstr,"X",1)==0) snpchr=23;
            else if(strncmp(tmpstr,"Y",1)==0) snpchr=24;
            else snpchr=atoi(tmpstr);
            FLAG_VALID_CK("--snp-chr", tmpstr);
            cout<<"--snp-chr "<<tmpstr<<endl;
        }
        else if(strcmp(option_str[i],"--snp-rm")==0){
            snp2rm=option_str[++i];
            FLAG_VALID_CK("--snp-rm", snp2rm);
            printf("--snp-rm %s\n", snp2rm);
        }
        else if (0 == strcmp(option_str[i], "--snp")){
            snprs = option_str[++i];
            FLAG_VALID_CK("--snp", snprs);
            printf("--snp %s\n", snprs);
        }
        else if (0 == strcmp(option_str[i], "--from-snp")){
            fromsnprs = option_str[++i];
            FLAG_VALID_CK("--from-snp", fromsnprs);
            printf("--from-snp %s\n", fromsnprs);
        }
        else if (0 == strcmp(option_str[i], "--to-snp")){
            tosnprs = option_str[++i];
            FLAG_VALID_CK("--to-snp", tosnprs);
            printf("--to-snp %s\n", tosnprs);
        }
        else if (0 == strcmp(option_str[i], "--probe")){
            prbname = option_str[++i];
            FLAG_VALID_CK("--probe", prbname);
            printf("--probe %s\n", prbname);
        }
        else if (0 == strcmp(option_str[i], "--from-probe")){
            fromprbname = option_str[++i];
            FLAG_VALID_CK("--from-probe", fromprbname);
            printf("--from-probe %s\n", fromprbname);
        }
        else if (0 == strcmp(option_str[i], "--to-probe")){
            toprbname = option_str[++i];
            FLAG_VALID_CK("--to-probe", toprbname);
            printf("--to-probe %s\n", toprbname);
        }
        else if(strcmp(option_str[i],"--snp-wind")==0){
            snpwindFlag=true;
            char* tmpstr=option_str[++i];
            if(tmpstr==NULL || SMRDATA::has_prefix(tmpstr, "--")) i--;
            else snpWind=atoi(tmpstr);
            printf("--snp-wind %d Kb\n", snpWind);
            if(snpWind<0 )
            {
                fprintf (stderr, "Error: --snp-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--probe-wind")==0){
            prbwindFlag=true;
            char* tmpstr=option_str[++i];
            if(tmpstr==NULL || SMRDATA::has_prefix(tmpstr, "--")) i--;
            else prbWind=atoi(tmpstr);
            printf("--probe-wind %d Kb\n", prbWind);
            if(prbWind<0 )
            {
                fprintf (stderr, "Error: --probe-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--gene")){
            genename = option_str[++i];
            FLAG_VALID_CK("--gene", genename);
            printf("--gene %s\n", genename);
        }
        else if(strcmp(option_str[i],"--from-snp-kb")==0){
            fromsnpkb=atoi(option_str[++i]);
            printf("--from-snp-kb %d Kb\n", fromsnpkb);
            if(fromsnpkb<0 )
            {
                fprintf (stderr, "Error: --from-snp-kb should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--to-snp-kb")==0){
            tosnpkb=atoi(option_str[++i]);
            printf("--to-snp-kb %d Kb\n", tosnpkb);
            if(tosnpkb<0 )
            {
                fprintf (stderr, "Error: --to-snp-kb should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--from-probe-kb")==0){
            fromprbkb=atoi(option_str[++i]);
            printf("--from-probe-kb %d Kb\n", fromprbkb);
            if(fromprbkb<0 )
            {
                fprintf (stderr, "Error: --from-probe-kb should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--to-probe-kb")==0){
            toprbkb=atoi(option_str[++i]);
            printf("--to-probe-kb %d Kb\n", toprbkb);
            if(toprbkb<0 )
            {
                fprintf (stderr, "Error: --to-probe-kb should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--extract-single-exposure-probe")==0){
            eprobe=option_str[++i];
            FLAG_VALID_CK("--extract-single-exposure-probe", eprobe);
            cout<<"--extract-single-exposure-probe "<<eprobe<<endl;
        }
        else if(strcmp(option_str[i],"--extract-single-outcome-probe")==0){
            oprobe=option_str[++i];
            FLAG_VALID_CK("--extract-single-outcome-probe", oprobe);
            cout<<"--extract-single-outcome-probe "<<oprobe<<endl;
        }
        else if(strcmp(option_str[i],"--exclude-single-exposure-probe")==0){
            eprobe2rm=option_str[++i];
            FLAG_VALID_CK("--exclude-single-exposure-probe", eprobe2rm);
            cout<<"--exclude-single-exposure-probe "<<eprobe2rm<<endl;
        }
        else if(strcmp(option_str[i],"--exclude-single-outcome-probe")==0){
            oprobe2rm=option_str[++i];
            FLAG_VALID_CK("--exclude-single-outcome-probe", oprobe2rm);
            cout<<"--exclude-single-outcome-probe "<<oprobe2rm<<endl;
        }
        else if(0==strcmp(option_str[i],"--qfile")){
            queryfileflg=true;
            queryFileName=option_str[++i];
            FLAG_VALID_CK("--qfile", queryFileName);
            printf("--qfile %s\n",queryFileName);
        }
        else if(strcmp(option_str[i],"--beqtl-qc")==0){
            qcflag=true;
            qcmtd = atoi(option_str[++i]);
            if(qcmtd<0 || qcmtd>2)
            {
                fprintf (stderr, "Error: --beqtl-qc should be 0,1 or 2.\n"); // 1 for SD, 2 for MAD
                exit (EXIT_FAILURE);
            }
            cout<<"--beqtl-qc "<<qcmtd<<endl;
        }
        else if(strcmp(option_str[i],"--z-thresh")==0){
           
            z_thresh = atoi(option_str[++i]);
            if(z_thresh<0 )
            {
                fprintf (stderr, "Error: --z-thresh should be positive.\n");
                exit (EXIT_FAILURE);
            }
            cout<<"--z-thresh "<<z_thresh<<endl;
        }
        else if (0 == strcmp(option_str[i], "--psmr")){
            threshpsmrest = atof(option_str[++i]);
            printf("--psmr %10.2e\n", threshpsmrest);
            if(threshpsmrest<0 || threshpsmrest>1)
            {
                fprintf (stderr, "Error: --psmr should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--phet")){
            threshphet = atof(option_str[++i]);
            printf("--phet %10.2e\n", threshphet);
            if(threshphet<0 || threshphet>1)
            {
                fprintf (stderr, "Error: --phet should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(0==strcmp(option_str[i],"--descriptive-cis")){
            count_cis_flag=true;
            printf("--descriptive-cis\n");
        }
        else if(0==strcmp(option_str[i],"--descriptive-trans")){
            count_trans_flag=true;
            if(!cis_flag) cis_itvl=5000;
            printf("--descriptive-trans\n");
        }
        else if(0==strcmp(option_str[i],"--opt-selection")){
            opt_slct_flag=true;
            printf("--opt-selection\n");
        }
        else if (0 == strcmp(option_str[i], "--extract-cis")){
            extract_cis_only=true;
            printf("--extract-cis \n");
        }
        else if (0 == strcmp(option_str[i], "--rm-technical")){
            rm_technical=true;
            prbseqregion=option_str[++i];
            printf("--rm-technical %s\n",prbseqregion);
        }
        else if (0 == strcmp(option_str[i], "--p-technical")){
            ptech = atof(option_str[++i]);
            printf("--p-technical %10.2e\n", ptech);
            if(ptech<0 || ptech>1)
            {
                fprintf (stderr, "Error: --p-technical should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--ld-lower-limit")){
            ld_min = atof(option_str[++i]);
            printf("--ld-lower-limit %f\n", ld_min);
            if(ld_min<0 || ld_min>1)
            {
                fprintf (stderr, "Error: --ld-lower-limit should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--add-n")==0){
            addn=atoi(option_str[++i]);
            printf("--add-n %d\n", addn);
            if(addn<0 )
            {
                fprintf (stderr, "Error: --add-n should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--show-n")==0){
            shownflag=true;
            printf("--show-n \n");
        }
        else if(strcmp(option_str[i],"--update-esi")==0){
            refesiName=option_str[++i];
            FLAG_VALID_CK("--update-esi", refesiName);
            cout<<"--update-esi "<<refesiName<<endl;
            FileExist(refesiName);
        }
        else if(strcmp(option_str[i],"--update-epi")==0){
            refepiName=option_str[++i];
            FLAG_VALID_CK("--update-epi", refepiName);
            cout<<"--update-epi "<<refepiName<<endl;
            FileExist(refepiName);
        }
        else if(strcmp(option_str[i],"--cis-to-all")==0){
            cis2all=true;
            printf("--cis-to-all \n");
        }
        else if(0==strcmp(option_str[i],"--sample-overlap")){
            sampleoverlap=true;
            printf("--sample-overlap \n");
        }
        else if(0==strcmp(option_str[i],"--mmecs")){
            minsnpcor=atoi(option_str[++i]);
            printf("--mmecs \n");
        }
        else if(0==strcmp(option_str[i],"--extract-target-snp-probe")){
            targetsnpproblstName=option_str[++i];
            FLAG_VALID_CK("--extract-target-snp-probe", targetsnpproblstName);
             printf("--extract-target-snp-probe %s\n",targetsnpproblstName);
            FileExist(targetsnpproblstName);
        }
        else if(0==strcmp(option_str[i],"--extract-snp-probe")){
            snpproblstName=option_str[++i];
            FLAG_VALID_CK("--extract-snp-probe", snpproblstName);
            printf("--extract-snp-probe %s\n",snpproblstName);
            FileExist(snpproblstName);
        }
        else if(0==strcmp(option_str[i],"--disable-freq-ck")){
            forcefrqck=false;
            printf("--disable-freq-ck\n");
        }
        else if(0==strcmp(option_str[i],"--diff-freq")){
            diff_freq=atof(option_str[++i]);
            if(diff_freq<0 || diff_freq>1)
            {
                printf ("Error: --diff-freq should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
            printf("--diff-freq %f\n",diff_freq);
        }
        else if(0==strcmp(option_str[i],"--diff-freq-prop")){
            diff_freq_ratio=atof(option_str[++i]);
            if(diff_freq_ratio<0 || diff_freq_ratio>1)
            {
                printf ("Error: --diff-freq-prop should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
            printf("--diff-freq-prop %f\n",diff_freq_ratio);
        }
        else if(0==strcmp(option_str[i],"--r")){
            ldr=true;
            printf("--r\n");
        }
        else if(0==strcmp(option_str[i],"--r2")){
            ldr2=true;
            printf("--r2\n");
        }
        else if(0==strcmp(option_str[i],"--ld-wind")){
            ldWind=atoi(option_str[++i]);
            printf("--ld-wind %d kb\n", ldWind);
            if(ldWind<=0 )
            {
                fprintf (stderr, "Error: --ld-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--trev")){
            trev = atof(option_str[++i]);
            printf("--trev %f\n", trev);
            if(trev>0)
            {
                fprintf (stderr, "Error: --trev should be negative or 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--qtltools-nominal-format")){
            qtltoolsnflag = true;
            gctaflag = false;
            printf("--qtltools-nominal-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--qtltools-permu-format")){
            qtltoolspflag = true;
            gctaflag = false;
            printf("--qtltools-permu-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--beqtl-probe-match")){
            prbmatch=true;
            printf("--beqtl-probe-match \n");
        }
        else if (0 == strcmp(option_str[i], "--beqtl-cis")){
            e2e_cis=true;
            printf("--beqtl-cis \n");
        }
        else if (0 == strcmp(option_str[i], "--mediation-probe-match")){
            mediation=true;
            printf("--mediation-probe-match \n");
        }
        else if (0 == strcmp(option_str[i], "--mediation-cis")){
            mediation_multiple_cis=true;
            printf("--mediation-cis \n");
        }
        else if (0 == strcmp(option_str[i], "--mediation-trans")){
            mediation_multiple_trans=true;
            printf("--mediation-trans \n");
        }
        else if (0 == strcmp(option_str[i], "--calcul-cormat")){
            calcul_cormat_flag=true;
            printf("--calcul-cormat \n");
        }
        else if (0 == strcmp(option_str[i], "--core-mediators")){
            core_med_lst=option_str[++i];
            printf("--core-mediators %s\n", core_med_lst);
            FileExist(core_med_lst);
        }  
        else if(strcmp(option_str[i],"--num-med")==0){
            num_med=atoi(option_str[++i]);
            printf("--num-med %d\n", num_med);
            if(num_med<1 )
            {
                fprintf (stderr, "Error: --num-med should be at least 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--uncorr-med")==0){
            med_R_thresh=atof(option_str[++i]);
            uncorr_med = true;
            cout<<"--uncorr-med "<<med_R_thresh<<endl;
            if(med_R_thresh<0 || med_R_thresh>1)
            {
                fprintf (stderr, "Error: --uncorr-med should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--mediator-qtls")){
            qtl_med_lst=option_str[++i];
            printf("--mediator-qtls %s\n", qtl_med_lst);
            FileExist(qtl_med_lst);
        }  
        else if (0 == strcmp(option_str[i], "--mediator-corr-snps")){
            med_corr_snps_lst=option_str[++i];
            printf("--mediator-corr-snps %s\n", med_corr_snps_lst);
            FileExist(qtl_med_lst);
        }  
        else if(strcmp(option_str[i],"--p-shrinkage")==0){
            p_shrinkage=atof(option_str[++i]);
            printf("--p-shrinkage %10.2e\n", p_shrinkage);
            if(p_shrinkage<0 || p_shrinkage>1)
            {
                fprintf (stderr, "Error: --p-shrinkage should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--p-expo-med")==0){
            p_expo_med=atof(option_str[++i]);
            printf("--p-expo-med %10.2e\n", p_expo_med);
            if(p_expo_med<0 || p_expo_med>1)
            {
                fprintf (stderr, "Error: --p-expo-med should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }

        else if (0 == strcmp(option_str[i], "--incomplete-expo")){
            incmp_expo=true;
            printf("--incomplete-expo \n");
        }
        else if (0 == strcmp(option_str[i], "--multi-uncorr-snp")){
            multi_corr_snp_flg=false;
            printf("--multi-uncorr-inst\n");
        }
        else if (0 == strcmp(option_str[i], "--get-snp-effects")){
            get_snp_effects_flg=true;
            printf("--get-snp-effects\n");
        }
        else if(strcmp(option_str[i],"--min-snp")==0){
            min_snp=atoi(option_str[++i]);
            printf("--min-snp %d\n", min_snp);
            if(min_snp<1 )
            {
                fprintf (stderr, "Error: --min-snp should be at least 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        
    }
    
#ifndef __APPLE__
#if defined _WIN64 || defined _WIN32
    omp_set_num_threads(thread_num);
#else
    stringstream ss;
    ss << thread_num;
    setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
    omp_set_num_threads(thread_num);
    
#endif
#endif
    
    cout<<endl;
    
    char tmpch[4]="smr";
    if(outFileName == NULL) outFileName=tmpch;
    if(shownflag) shown(eqtlFileName);
    else if(make_bld_flag) ld_report(outFileName, bFileName,indilstName, indilst2remove, snplstName,  snplst2exclde,chr, snprs, maf,ldr,ldr2,ldWind);
    else if(make_besd_flag && (mateqtlflag || fastqtlnflag || fastqtlpflag || qtltoolsnflag || qtltoolspflag) )  make_besd_fmat(eqtlFileName,outFileName,mateqtlflag, fastqtlnflag, fastqtlpflag,qtltoolsnflag,qtltoolspflag,addn);
    else if(make_besd_flag && (gctaflag || plinkflag || gemmaflag || boltflag || merlinflag) ) make_besd(outFileName, syllabusName, gctaflag, plinkflag, gemmaflag,merlinflag,boltflag,save_dense_flag, cis_itvl,  trans_itvl,  transThres, restThres,genouni,addn);
    else if (make_besd_flag && queryfileflg) make_besd_byQfile(queryFileName,outFileName, save_dense_flag, cis_itvl,  trans_itvl,  transThres,  restThres,addn);
    else if (metaflg) meta( eqtlsmaslstName ,outFileName,meta_mtd,pmecs, cis_flag,  cis_itvl, nmecs,problstName,  problst2exclde,  genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,snplstName,  snplst2exclde, snpchr, snprs,  fromsnprs, tosnprs, snpWind, fromsnpkb,  tosnpkb,  snpwindFlag);
    else if (combineFlg && !metaflg) combineBesd(eqtlsmaslstName, outFileName, save_dense_flag, cis_itvl,trans_itvl, transThres, restThres, genouni,addn);
    else if(save_dense_flag) make_full_besd(outFileName, eqtlFileName, snplstName,problstName,bFlag,make_besd_flag, snplst2exclde, problst2exclde,  cis_itvl,genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,addn);
    else if(make_besd_flag) make_sparse_besd(eqtlFileName, outFileName, cis_itvl,trans_itvl, transThres, restThres,genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag, snplstName, problstName,  snplst2exclde,  problst2exclde, qcflag ,qcmtd, z_thresh,extract_cis_only,prbseqregion,ptech,pinsnp,pexsnp,addn);
    else if(plotflg && eqtlFileName2 != NULL) plot_triple( outFileName,  bFileName, gwasFileName,  eqtlFileName2,eqtlFileName, maf, indilstName,  snplstName, p_hetero, ld_prune, m_hetero ,opt_hetero, indilst2remove,  snplst2exclde,  p_smr,  refSNP,  cis_itvl,  prbname,  prbWind, prbwindFlag,  snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  geneAnnoName,  pthres_me2esmr,threshpsmrest,new_het_mth,threshphet,opt_slct_flag,ld_min,sampleoverlap,pmecs,minsnpcor,targetsnpproblstName, diff_freq,diff_freq_ratio);
    else if(plotflg) plot_newheidi(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,p_hetero,ld_prune,m_hetero, opt_hetero,indilst2remove, snplst2exclde, problst2exclde,p_smr,refSNP, heidioffFlag,cis_itvl, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,geneAnnoName,threshpsmrest,new_het_mth, threshphet,ld_min,sampleoverlap,pmecs,minsnpcor, diff_freq,diff_freq_ratio);
    else if (calcul_cormat_flag){
        if (eqtlFileName2 != NULL) calcul_cormat(outFileName, eqtlFileName, eqtlFileName2, snplstName, problstName, genelistName, chr, prbchr, snpchr, snprs, fromsnprs, tosnprs, prbname, fromprbname, toprbname, snpWind, prbWind, genename, fromsnpkb, tosnpkb, fromprbkb, toprbkb, snpwindFlag, prbwindFlag, cis_flag, cis_itvl);
        else calcul_cormat(outFileName, eqtlFileName, snplstName, problstName, genelistName, chr, prbchr, snpchr, snprs, fromsnprs, tosnprs, prbname, fromprbname, toprbname, snpWind, prbWind, genename, fromsnpkb, tosnpkb, fromprbkb, toprbkb, snpwindFlag, prbwindFlag, cis_flag, cis_itvl);
    }
    else if((eqtlFileName2 != NULL) & !prbmatch & !e2e_cis & !mediation & !mediation_multiple_cis & !mediation_multiple_trans) smr_e2e(outFileName, bFileName,eqtlFileName2, eqtlFileName, maf,indilstName, snplstName,problstName,oproblstName,eproblstName,bFlag,p_hetero,ld_prune,m_hetero, opt_hetero,indilst2remove, snplst2exclde, problst2exclde,oproblst2exclde,eproblst2exclde,p_smr,refSNP, heidioffFlag,cis_itvl,traitlstName,outcomePrbWind,oprobe, eprobe, oprobe2rm, eprobe2rm,threshpsmrest,new_het_mth,opt_slct_flag,ld_min,cis2all,sampleoverlap,pmecs,minsnpcor,ssmrflg, setWind,ld_prune_multi,targetsnpproblstName,snpproblstName, diff_freq,diff_freq_ratio);
    else if((eqtlFileName2 != NULL) & prbmatch) smr_e2e_prbmatch(outFileName, bFileName,eqtlFileName2, eqtlFileName, maf,indilstName, snplstName,problstName,oproblstName,eproblstName,indilst2remove, snplst2exclde, problst2exclde,oproblst2exclde,eproblst2exclde,p_smr, cis_itvl,traitlstName,oprobe, eprobe, oprobe2rm, eprobe2rm,cis2all,ld_prune_multi, ldmatrix, trev, snpproblstName, diff_freq,diff_freq_ratio, get_snp_effects_flg, min_snp);
    else if((eqtlFileName2 != NULL) & e2e_cis) smr_e2e_cis(outFileName, bFileName,eqtlFileName2, eqtlFileName, maf,indilstName, snplstName,problstName,oproblstName,eproblstName,indilst2remove, snplst2exclde, problst2exclde,oproblst2exclde,eproblst2exclde,p_smr, cis_itvl,traitlstName,oprobe, eprobe, oprobe2rm, eprobe2rm,cis2all,ld_prune_multi,ldmatrix, trev, snpproblstName, diff_freq,diff_freq_ratio,multi_corr_snp_flg, get_snp_effects_flg, min_snp);
    else if((eqtlFileName2 != NULL) & mediation) smr_mediation_prbmatch(outFileName, bFileName,eqtlFileName2, eqtlFileName, gwasFileName, maf,indilstName, snplstName,problstName,oproblstName,eproblstName,indilst2remove, snplst2exclde, problst2exclde,oproblst2exclde,eproblst2exclde,p_smr, p_medsmr, cis_itvl,traitlstName,oprobe, eprobe, oprobe2rm, eprobe2rm,cis2all,ld_prune_multi, ldmatrix, trev, snpproblstName, diff_freq,diff_freq_ratio,multi_corr_snp_flg, get_snp_effects_flg, min_snp);
    else if((eqtlFileName2 != NULL) & mediation_multiple_cis) smr_multi_mediation_cis(outFileName, bFileName,eqtlFileName2, eqtlFileName, gwasFileName, maf,indilstName, snplstName,problstName,oproblstName,eproblstName,indilst2remove, snplst2exclde, problst2exclde,oproblst2exclde,eproblst2exclde,p_smr, p_medsmr, cis_itvl,traitlstName,oprobe, eprobe, oprobe2rm, eprobe2rm,cis2all,ld_prune_multi, ldmatrix, trev, diff_freq,diff_freq_ratio,multi_corr_snp_flg, get_snp_effects_flg, min_snp, incmp_expo, med_R_thresh, uncorr_med, p_shrinkage, p_expo_med);
    else if((eqtlFileName2 != NULL) & mediation_multiple_trans) smr_multi_mediation_trans(outFileName, bFileName, eqtlFileName2, eqtlFileName, gwasFileName, maf,indilstName, snplstName,problstName,oproblstName,eproblstName,indilst2remove, snplst2exclde, problst2exclde,oproblst2exclde,eproblst2exclde,p_smr, p_medsmr, cis_itvl,oprobe, eprobe, oprobe2rm, eprobe2rm, ld_prune_multi, ldmatrix, trev, snpproblstName, diff_freq,diff_freq_ratio, get_snp_effects_flg, min_snp, incmp_expo, core_med_lst, qtl_med_lst, med_corr_snps_lst, num_med, med_R_thresh, uncorr_med, p_shrinkage, p_expo_med);
    else if(lookup_flag) {
        if(bldFileName)  lookup(outFileName, bldFileName,  snplstName, snplst2exclde, chr, snprs, snp2rm,  fromsnprs,  tosnprs, snpWind,  snpwindFlag,  fromsnpkb,  tosnpkb, ldWind);
        else if(trans_dense) lookup_dense(outFileName,eqtlFileName, snplstName, problstName, genelistName, plookup, chr, prbchr,snpchr, snprs, fromsnprs, tosnprs, prbname, fromprbname, toprbname,snpWind,prbWind,genename,fromsnpkb,tosnpkb,fromprbkb, toprbkb, snpwindFlag, prbwindFlag,cis_flag, cis_itvl,snpproblstName);
        else lookup(outFileName,eqtlFileName, snplstName, problstName, genelistName, plookup, bFlag, chr, prbchr,snpchr, snprs, fromsnprs, tosnprs, prbname, fromprbname, toprbname,snpWind,prbWind,genename,fromsnpkb,tosnpkb,fromprbkb, toprbkb, snpwindFlag, prbwindFlag,cis_flag, cis_itvl,snpproblstName);
    }
    else if(count_cis_flag) count_cis(outFileName,eqtlFileName, cisThres, cis_itvl); //output top-SNP in the cis-regions
    else if(count_trans_flag) count_trans(outFileName,eqtlFileName, transThres, cis_itvl, trans_itvl); // output top-SNP in the trans-regions. trans-region is defined as the region aside of the cis_itvl
    else if (freqName) update_freq(eqtlFileName, freqName);
    else if(refepiName || refesiName) {
        if(refepiName) update_epifile(eqtlFileName, refepiName);
        if(refesiName) update_esifile(eqtlFileName, refesiName);
    }
    else if(recodeflg) make_cojo(outFileName, eqtlFileName, snplstName, snplst2exclde,  problstName,  problst2exclde,  genelistName, bFlag);
    else if(gwasFileName2 != NULL && ssmrgwastogwasflg)  smr_ivw_gwas_gwas(outFileName, bFileName, gwasFileName, gwasFileName2, maf,indilstName, snplstName,problstName, indilst2remove, snplst2exclde, problst2exclde, p_smr, cis_itvl, genelistName, chr, prbchr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, snpchr, snprs, fromsnprs, tosnprs, snpWind, fromsnpkb, tosnpkb, snpwindFlag, cis_flag, setlstName, geneAnnoName, ld_prune_multi, ldmatrix, trev, diff_freq, diff_freq_ratio, get_snp_effects_flg, min_snp);
    else if(ssmrflg)  smr_multipleSNP(outFileName, bFileName, gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName, bFlag, p_hetero, ld_prune, m_hetero ,opt_hetero, indilst2remove, snplst2exclde, problst2exclde, p_smr, refSNP, heidioffFlag, cis_itvl, genelistName, chr, prbchr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, snpchr, snprs, fromsnprs, tosnprs, snpWind, fromsnpkb, tosnpkb, snpwindFlag, cis_flag, setlstName, geneAnnoName, setWind, ld_min,threshpsmrest,  sampleoverlap,  pmecs,  minsnpcor, ld_prune_multi, diff_freq, diff_freq_ratio);
    else if(ssmrivwflg)  smr_ivw_analysis(outFileName, bFileName, gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName, indilst2remove, snplst2exclde, problst2exclde, p_smr, cis_itvl, genelistName, chr, prbchr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, snpchr, snprs, fromsnprs, tosnprs, snpWind, fromsnpkb, tosnpkb, snpwindFlag, cis_flag, setlstName, geneAnnoName, ld_prune_multi, ldmatrix, trev, diff_freq, diff_freq_ratio, get_snp_effects_flg, min_snp);
    else if(ssmrrevflg)  smr_rev_ivw_analysis(outFileName, bFileName, gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName, indilst2remove, snplst2exclde, problst2exclde, p_gwas, cis_itvl, genelistName, chr, prbchr, prbname, fromprbname, toprbname, prbWind, fromprbkb, toprbkb, prbwindFlag, genename, snpchr, snprs, fromsnprs, tosnprs, snpWind, fromsnpkb, tosnpkb, snpwindFlag, cis_flag, setlstName, geneAnnoName, ld_prune_multi, ldmatrix, trev, diff_freq, diff_freq_ratio, get_snp_effects_flg, min_snp);
    else if(smr_flag && !smr_trans_flag) smr(outFileName, bFileName,bldFileName, gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,p_hetero,ld_prune,m_hetero,opt_hetero, indilst2remove, snplst2exclde, problst2exclde,p_smr,refSNP, heidioffFlag,cis_itvl, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,threshpsmrest,new_het_mth,opt_slct_flag,prbseqregion,ld_min,sampleoverlap,pmecs,minsnpcor,targetsnpproblstName,snpproblstName, diff_freq,diff_freq_ratio);
    else if (est_effe_spl_size_flg) est_effect_splsize(eqtlsmaslstName,eqtlFileName, snplstName,problstName,snplst2exclde, problst2exclde,p_smr);
    else if(smr_flag && smr_trans_flag) smr_trans_region(outFileName, bFileName,gwasFileName, eqtlFileName,  maf, indilstName, snplstName,problstName, bFlag, p_hetero, ld_prune, m_hetero ,opt_hetero,indilst2remove,snplst2exclde, problst2exclde,  transThres, refSNP,  heidioffFlag, cis_itvl, trans_itvl,genelistName,  chr, prbchr, prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag, genename,snpchr,snprs,fromsnprs, tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag, threshpsmrest,  new_het_mth, p_smr,opt_slct_flag,ld_min, diff_freq,diff_freq_ratio);
    
   }
