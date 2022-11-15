/*
 * SMR_data_p4.cpp
 * Implementations of univariable and multivariable MR-IVW implementations
 *
 * 2020 by Marie Sadler
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  
 */


#include "SMR_data_p4.h"
namespace SMRDATA
{
    void lookup_dense(char* outFileName,char* eqtlFileName, char* snplstName, char* problstName,char* genelistName, double plookup, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl,char* snpproblstName)
    {
        string logstr;     
        map<string, string> prb_snp;
        map<string, string>::iterator iter;
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading QTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
        }

        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
        
        uint64_t epinum = eqtlinfo._probNum;
        uint64_t esinum = eqtlinfo._snpNum;

        vector<uint64_t> e_snp_idx;
        vector<uint64_t> e_prb_idx;

        for (int i = 0; i < eqtlinfo._esi_include.size(); i++) e_snp_idx.push_back(eqtlinfo._snp_name_map[eqtlinfo._esi_rs[eqtlinfo._esi_include[i]]]);
        for (int i = 0; i < eqtlinfo._include.size(); i++) e_prb_idx.push_back(eqtlinfo._probe_name_map[eqtlinfo._epi_prbID[eqtlinfo._include[i]]]);

        // get the effect sizes

        vector<int> out_esi_id;
        vector<int> out_epi_id;
        vector<float> out_beta;
        vector<float> out_se;
        vector<double> out_pval;

        
        FILE * fptr;
        string besd_file = string(eqtlFileName)+".besd";
        fptr = fopen(besd_file.c_str(),"rb");
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        int infoLen=sizeof(uint32_t);
        if(indicator==DENSE_FILE_TYPE_3) infoLen=RESERVEDUNITS*sizeof(int);

        if(indicator!=DENSE_FILE_TYPE_1 && indicator!=DENSE_FILE_TYPE_3) {
            cout << "ERROR: the BESD file is not in dense BESD format." << endl;
            exit(EXIT_FAILURE);
        }
        

        float* tmp=(float*)malloc(sizeof(float)*esinum<<1);

        for (int i = 0; i < eqtlinfo._include.size(); i++){

            fseek(fptr,((e_prb_idx[i]*esinum)<<3)+infoLen, SEEK_SET);
            fread(tmp, sizeof(float), esinum*2,fptr);

            for (int j = 0; j< eqtlinfo._esi_include.size(); j++){

                double beta=tmp[e_snp_idx[j]];
                double se=tmp[e_snp_idx[j] + esinum];
                if(fabs(se+9)<1e-6) continue;
                double zsxz=beta/se;
                double pxz=pchisq(zsxz*zsxz, 1);
                if(pxz<=plookup){

                    out_esi_id.push_back(eqtlinfo._esi_include[j]);
                    out_epi_id.push_back(eqtlinfo._include[i]);
                    out_beta.push_back(beta);
                    out_se.push_back(se);
                    out_pval.push_back(pxz);
                }
            }
        }

        fclose(fptr);

        string smrfile = string(outFileName)+".txt";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
        
        smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' <<"Freq"<<'\t'<< "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"Orientation"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_bp[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<((eqtlinfo._esi_freq[out_esi_id[i]]+9>1e-6)?atos(eqtlinfo._esi_freq[out_esi_id[i]]):"NA")<<'\t'<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_chr[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_bp[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_gene[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_orien[out_epi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
        }
        
        smr.close();
        cout<<"Extracted results of "<<out_esi_id.size()<<" SNPs have been saved in the file [" + smrfile + "]."<<endl;
      
    }

    void allele_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        vector<int> mdId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the exposure QTL summary data, the mediation QTL summary data, the GWAS summary dataset and the LD reference data).\n";
        cout<<logstr<<endl;
        
        vector<string> bsnp(bdata->_include.size());
        vector<string> gsnp(gdata->_include.size());
        vector<string> msnp(mdata->_esi_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum || mdata->_esi_include.size()<mdata->_snpNum )
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<gdata->_include.size();i++)
                gsnp[i]=gdata->snpName[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<mdata->_esi_include.size();i++)
                msnp[i]=mdata->_esi_rs[mdata->_esi_include[i]];
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(bsnp, msnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and mediation data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=msnp[edId[i]];
            edId.clear();
            match_only(cmmnSNPs, gsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gsnp[edId[i]];
            edId.clear();
            match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=essnp[edId[i]];
        }else
        {
            match_only(bdata->_snp_name, mdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and mediation data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=mdata->_esi_rs[edId[i]];
            edId.clear();
            match_only(cmmnSNPs, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gdata->snpName[edId[i]];
            edId.clear();
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        match(slctSNPs, bdata->_snp_name, bdId);
        match(slctSNPs, mdata->_esi_rs, mdId);
        match(slctSNPs, gdata->snpName, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        mdata->_esi_include.clear();
        esdata->_esi_include.clear();
        double disp=0;
        for (int i = 0; i<edId.size(); i++)
        {
            progress(i, disp, (int)edId.size());
           
            string a1, a2, ga1, ga2, ea1, ea2, ma1, ma2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ma1 = mdata->_esi_allele1[mdId[i]];
            ma2 = mdata->_esi_allele2[mdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
         
            // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memory
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ma1 && ea2 == ma2 && ea1 == ga1 && ea2 == ga2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ma2 && ea2 == ma1 && ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=mdata->_esi_allele1[mdId[i]];
                    mdata->_esi_allele1[mdId[i]] = mdata->_esi_allele2[mdId[i]];
                    mdata->_esi_allele2[mdId[i]] = tmpch;
                    if(mdata->_esi_freq[mdId[i]] >0) mdata->_esi_freq[mdId[i]] = 1 - mdata->_esi_freq[mdId[i]];
                    
                    if(mdata->_val.size()>0)
                    {
                        
                        int count=0;
                        for(int j=0;j<mdata->_rowid.size();j++)
                        {
                            if(mdata->_rowid[j]==mdId[i])
                            {
                                count++;
                                if(count & 1)
                                    mdata->_val[j]=-1*mdata->_val[j];
                            }
                        }
                    }
                    else
                    {
                        int mid=mdId[i];
                        #pragma omp parallel for
                        for(int j=0;j<mdata->_include.size();j++)
                            if( fabs(mdata->_sexz[j][mid]+9) > 1e-6 )
                                mdata->_bxz[j][mid]=-1*mdata->_bxz[j][mid];
                    }
                }
                else if(ea1 == ma1 && ea2 == ma2 && ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                }
                else if(ea1 == ma2 && ea2 == ma1 && ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);

                    string tmpch=mdata->_esi_allele1[mdId[i]];
                    mdata->_esi_allele1[mdId[i]] = mdata->_esi_allele2[mdId[i]];
                    mdata->_esi_allele2[mdId[i]] = tmpch;
                    if(mdata->_esi_freq[mdId[i]] >0) mdata->_esi_freq[mdId[i]] = 1 - mdata->_esi_freq[mdId[i]];
                    
                    if(mdata->_val.size()>0)
                    {                       
                        int count=0;
                        for(int j=0;j<mdata->_rowid.size();j++)
                        {
                            if(mdata->_rowid[j]==mdId[i])
                            {
                                count++;
                                if(count & 1)
                                    mdata->_val[j]=-1*mdata->_val[j];
                            }
                        }
                    }
                    else
                    {
                        int mid=mdId[i];
                        #pragma omp parallel for
                        for(int j=0;j<mdata->_include.size();j++)
                            if( fabs(mdata->_sexz[j][mid]+9) > 1e-6 )
                                mdata->_bxz[j][mid]=-1*mdata->_bxz[j][mid];
                    }
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ma1 && ea2 == ma2 && ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ma2 && ea2 == ma1 && ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=mdata->_esi_allele1[mdId[i]];
                    mdata->_esi_allele1[mdId[i]] = mdata->_esi_allele2[mdId[i]];
                    mdata->_esi_allele2[mdId[i]] = tmpch;
                    if(mdata->_esi_freq[mdId[i]] >0) mdata->_esi_freq[mdId[i]] = 1 - mdata->_esi_freq[mdId[i]];
                    
                    if(mdata->_val.size()>0)
                    {
                        
                        int count=0;
                        for(int j=0;j<mdata->_rowid.size();j++)
                        {
                            if(mdata->_rowid[j]==mdId[i])
                            {
                                count++;
                                if(count & 1)
                                    mdata->_val[j]=-1*mdata->_val[j];
                            }
                        }
                    }
                    else
                    {
                        int mid=mdId[i];
                        #pragma omp parallel for
                        for(int j=0;j<mdata->_include.size();j++)
                            if( fabs(mdata->_sexz[j][mid]+9) > 1e-6 )
                                mdata->_bxz[j][mid]=-1*mdata->_bxz[j][mid];
                    }

                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ma1 && ea2 == ma2 && ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];

                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ma2 && ea2 == ma1 && ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);

                    string tmpch=mdata->_esi_allele1[mdId[i]];
                    mdata->_esi_allele1[mdId[i]] = mdata->_esi_allele2[mdId[i]];
                    mdata->_esi_allele2[mdId[i]] = tmpch;
                    if(mdata->_esi_freq[mdId[i]] >0) mdata->_esi_freq[mdId[i]] = 1 - mdata->_esi_freq[mdId[i]];
                    
                    if(mdata->_val.size()>0)
                    {
                        
                        int count=0;
                        for(int j=0;j<mdata->_rowid.size();j++)
                        {
                            if(mdata->_rowid[j]==mdId[i])
                            {
                                count++;
                                if(count & 1)
                                    mdata->_val[j]=-1*mdata->_val[j];
                            }
                        }
                    }
                    else
                    {
                        int mid=mdId[i];
                        #pragma omp parallel for
                        for(int j=0;j<mdata->_include.size();j++)
                            if( fabs(mdata->_sexz[j][mid]+9) > 1e-6 )
                                mdata->_bxz[j][mid]=-1*mdata->_bxz[j][mid];
                    }
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];

                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }          
        }
        
        logstr= "\n" + atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata, bool incmp_expo)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        vector<int> mdId;
        cmmnSNPs.clear();
        mdId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the exposure QTL summary data, the mediation QTL summary data, the GWAS summary dataset and the LD reference data).\n";
        cout<<logstr<<endl;
        
        vector<string> bsnp(bdata->_include.size());
        vector<string> gsnp(gdata->_include.size());
        vector<string> msnp(mdata->_esi_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum || mdata->_esi_include.size()<mdata->_snpNum )
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<gdata->_include.size();i++)
                gsnp[i]=gdata->snpName[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<mdata->_esi_include.size();i++)
                msnp[i]=mdata->_esi_rs[mdata->_esi_include[i]];
            match_only(bsnp, gsnp, mdId);
            if(mdId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            slctSNPs.resize(mdId.size());
            #pragma omp parallel for
            for(int i=0;i<mdId.size();i++)
                slctSNPs[i]=gsnp[mdId[i]];
            mdId.clear();
            match_only(slctSNPs, msnp, mdId);
            if(mdId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(mdId.size());
            #pragma omp parallel for
            for(int i=0;i<mdId.size();i++)
                slctSNPs[i]=msnp[mdId[i]];
            if (!incmp_expo){
                #pragma omp parallel for
                for(int i=0;i<esdata->_esi_include.size();i++)
                    essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
                mdId.clear();
                match_only(slctSNPs, essnp, mdId);
                if(mdId.empty()) throw("Error: no common SNPs found.");
                slctSNPs.resize(mdId.size());
                #pragma omp parallel for
                for(int i=0;i<mdId.size();i++)
                    slctSNPs[i]=essnp[mdId[i]];
                match(slctSNPs, mdata->_esi_rs, mdId);
            }

        }else
        {
            match_only(bdata->_snp_name, gdata->snpName, mdId);
            if(mdId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            slctSNPs.resize(mdId.size());
            #pragma omp parallel for
            for(int i=0;i<mdId.size();i++)
                slctSNPs[i]=gdata->snpName[mdId[i]];
            mdId.clear();
            match_only(slctSNPs, mdata->_esi_rs, mdId);
            if(mdId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(mdId.size());
            #pragma omp parallel for
            for(int i=0;i<mdId.size();i++)
                slctSNPs[i]=mdata->_esi_rs[mdId[i]];         
            if (!incmp_expo){
                mdId.clear();
                match_only(slctSNPs, esdata->_esi_rs, mdId);
                if(mdId.empty()) throw("Error: no common SNPs found.");
                slctSNPs.resize(mdId.size());
                #pragma omp parallel for
                for(int i=0;i<mdId.size();i++)
                    slctSNPs[i]=esdata->_esi_rs[mdId[i]];
                match(slctSNPs, mdata->_esi_rs, mdId);
            }
        }
        logstr= atos(slctSNPs.size())+" SNPs are included after allele check based on SNP ID. ";
        cout<<logstr<<endl;
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        match(slctSNPs, bdata->_snp_name, bdId);
        match(slctSNPs, gdata->snpName, gdId);
        //match(slctSNPs, mdata->_esi_rs, mdId);
        match(slctSNPs, esdata->_esi_rs, edId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        mdata->_esi_include.clear();
        esdata->_esi_include.clear();
        double disp=0;
        for (int i = 0; i<mdId.size(); i++)
        {
            progress(i, disp, (int)mdId.size());
           
            string a1, a2, ga1, ga2, ea1, ea2, ma1, ma2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ma1 = mdata->_esi_allele1[mdId[i]];
            ma2 = mdata->_esi_allele2[mdId[i]];
            if (edId[i] > -9){
                ea1 = esdata->_esi_allele1[edId[i]];
                ea2 = esdata->_esi_allele2[edId[i]];
            } else {
                ea1 = a1;
                ea2 = a2;
            }
         
            // use the allele in mediator summary data "mdata" as the reference allele. so we won't get the whole besd into memory
            if(ma1 == a1 &&  ma2 == a2)
            {
                if( ma1 == ea1 && ma2 == ea2 && ma1 == ga1 && ma2 == ga2)
                {                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ma1 == ea2 && ma2 == ea1 && ma1 == ga1 && ma2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    if (edId[i] > -9){        
                        string tmpch=esdata->_esi_allele1[edId[i]];
                        esdata->_esi_allele1[edId[i]] = esdata->_esi_allele2[edId[i]];
                        esdata->_esi_allele2[edId[i]] = tmpch;
                        if(esdata->_esi_freq[edId[i]] >0) esdata->_esi_freq[edId[i]] = 1 - esdata->_esi_freq[edId[i]];
                        
                        if(esdata->_val.size()>0)
                        {
                            int count=0;
                            for(int j=0;j<esdata->_rowid.size();j++)
                            {
                                if(esdata->_rowid[j]== edId[i])
                                {
                                    count++;
                                    if(count & 1)
                                        esdata->_val[j]=-1*esdata->_val[j];
                                }
                            }
                        }
                        else
                        {
                            int eid=edId[i];
                            #pragma omp parallel for
                            for(int j=0;j<esdata->_include.size();j++)
                                if(fabs(esdata->_sexz[j][eid]+9) > 1e-6 )
                                    esdata->_bxz[j][eid]=-1*esdata->_bxz[j][eid];
                        }
                    }
                }
                else if(ma1 == ea1 && ma2 == ea2 && ma1 == ga2 && ma2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                }
                else if(ma1 == ea2 && ma2 == ea1 && ma1 == ga2 && ma2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);

                    if (edId[i] > -9){        
                        string tmpch=esdata->_esi_allele1[edId[i]];
                        esdata->_esi_allele1[edId[i]] = esdata->_esi_allele2[edId[i]];
                        esdata->_esi_allele2[edId[i]] = tmpch;
                        if(esdata->_esi_freq[edId[i]] >0) esdata->_esi_freq[edId[i]] = 1 - esdata->_esi_freq[edId[i]];
                        
                        if(esdata->_val.size()>0)
                        {
                            int count=0;
                            for(int j=0;j<esdata->_rowid.size();j++)
                            {
                                if(esdata->_rowid[j]==edId[i])
                                {
                                    count++;
                                    if(count & 1)
                                        esdata->_val[j]=-1*esdata->_val[j];
                                }
                            }
                        }
                        else
                        {
                            int eid=edId[i];
                            #pragma omp parallel for
                            for(int j=0;j<esdata->_include.size();j++)
                                if(fabs(esdata->_sexz[j][eid]+9) > 1e-6 )
                                    esdata->_bxz[j][eid]=-1*esdata->_bxz[j][eid];
                        }
                    }
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                }
            }
            else if(ma1 == a2 &&  ma2 == a1)
            {
                
                if(ma1 == ea1 && ma2 == ea2 && ma1 == ga1 && ma2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ma1 == ea2 && ma2 == ea1 && ma1 == ga1 && ma2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    if (edId[i] > -9){
                        string tmpch=esdata->_esi_allele1[edId[i]];
                        esdata->_esi_allele1[edId[i]] = esdata->_esi_allele2[edId[i]];
                        esdata->_esi_allele2[edId[i]] = tmpch;
                        if(esdata->_esi_freq[edId[i]] >0) esdata->_esi_freq[edId[i]] = 1 - esdata->_esi_freq[edId[i]];
                        
                        if(esdata->_val.size()>0)
                        {
                            int count=0;
                            for(int j=0;j<esdata->_rowid.size();j++)
                            {
                                if(esdata->_rowid[j]==edId[i])
                                {
                                    count++;
                                    if(count & 1)
                                        esdata->_val[j]=-1*esdata->_val[j];
                                }
                            }
                        }
                        else
                        {
                            int eid=edId[i];
                            #pragma omp parallel for
                            for(int j=0;j<esdata->_include.size();j++)
                                if( fabs(esdata->_sexz[j][eid]+9) > 1e-6 )
                                    esdata->_bxz[j][eid]=-1*esdata->_bxz[j][eid];
                        }
                    }

                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ma1 == ea1 && ma2 == ea2 && ma1 == ga2 && ma2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];

                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ma1 == ea2 && ma2 == ea1 && ma1 == ga2 && ma2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    mdata->_esi_include.push_back(mdId[i]);
                    esdata->_esi_include.push_back(edId[i]);

                    if (edId[i] > -9){
                        string tmpch=esdata->_esi_allele1[edId[i]];
                        esdata->_esi_allele1[edId[i]] = esdata->_esi_allele2[edId[i]];
                        esdata->_esi_allele2[edId[i]] = tmpch;
                        if(esdata->_esi_freq[edId[i]] >0) esdata->_esi_freq[edId[i]] = 1 - esdata->_esi_freq[edId[i]];
                        
                        if(esdata->_val.size()>0)
                        {
                            int count=0;
                            for(int j=0;j<esdata->_rowid.size();j++)
                            {
                                if(esdata->_rowid[j]==edId[i])
                                {
                                    count++;
                                    if(count & 1)
                                        esdata->_val[j]=-1*esdata->_val[j];
                                }
                            }
                        }
                        else
                        {
                            int eid=edId[i];
                            #pragma omp parallel for
                            for(int j=0;j<esdata->_include.size();j++)
                                if( fabs(esdata->_sexz[j][eid]+9) > 1e-6 )
                                    esdata->_bxz[j][eid]=-1*esdata->_bxz[j][eid];
                        }
                    }
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];

                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }          
        }
        
        logstr= "\n" + atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // update _snp_name_map which will be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check(bInfo* bdata, gwasData* gdata, gwasData* gdata2)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the exposure GWAS summary data, the outcome GWAS summary data and the LD reference data).";
        cout<<logstr<<endl;
        
        vector<string> bsnp;
        if(bdata->_include.size()< bdata->_snp_num)
        {
            for(int i=0;i<bdata->_include.size();i++) bsnp.push_back(bdata->_snp_name[bdata->_include[i]]);
            StrFunc::match_only(bsnp, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no SNPs in common between reference data and exposure GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            StrFunc::match_only(cmmnSNPs, gdata2->snpName, edId);
            if(edId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(gdata2->snpName[edId[i]]);
        }else
        {
            StrFunc::match_only(bdata->_snp_name, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no SNPs in common  between reference data and GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            edId.clear();
            StrFunc::match_only(cmmnSNPs, gdata2->snpName, edId);
            if(edId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(gdata2->snpName[edId[i]]);
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        //logstr=itos(edId.size())+" SNPs are included after rsID check. ";
        //cout<<logstr<<endl;
        
        //alleles check
        StrFunc::match(slctSNPs, bdata->_snp_name, bdId);
        StrFunc::match(slctSNPs, gdata->snpName, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        gdata2->_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string a1, a2, ga1, ga2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = gdata2->allele_1[edId[i]];
            ea2 = gdata2->allele_2[edId[i]];
            // use the allele in outcome GWAS summary data as the reference allele

            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ga1 && ea2 == ga2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    gdata2->_include.push_back(edId[i]);
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    gdata2->_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    gdata2->_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                   
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    gdata2->_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }

            
        }
        
        logstr=itos(bdata->_include.size())+" SNPs are included after allele checking. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    double freq_check(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata, double &freqthresh, double &percenthresh)
    {
        cout << "Checking the consistency of allele frequency of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data, the mediation data and the LD reference data).\n" << endl;
        bool success=true;
        long failcount=0, snpnum=bdata->_include.size();
        vector<int> pasbid, pasgid, paseid, pasmid;
        string fname=string(outFileName) + ".snp_failed_freq_ck.list";
        FILE* tmpfile=fopen(fname.c_str(), "w");
        if (!(tmpfile)) {
            printf("Error: Failed to open file %s.\n",fname.c_str());
            exit(EXIT_FAILURE);
        }
        string tmpstr = "SNP\tAllele\tReferencePanel\tGWAS\teQTL\tmdata\n";
        fputs(tmpstr.c_str(),tmpfile);
        for(int i=0;i<bdata->_include.size();i++)
        {
            double bfreq=bdata->_mu[bdata->_include[i]]/2;
            double gfreq=gdata->freq[gdata->_include[i]];
            double efreq;
            if (esdata->_esi_include[i] > -9){
                efreq=esdata->_esi_freq[esdata->_esi_include[i]];
            } else {
                efreq = bdata->_mu[bdata->_include[i]]/2;
            }
            double mfreq;
            if (mdata->_esi_include[i] > -9){
                mfreq=mdata->_esi_freq[mdata->_esi_include[i]];
            } else {
                mfreq = bdata->_mu[bdata->_include[i]]/2;
            }
            if(gfreq<0 && efreq>0 && mfreq>0)
            {
                if(fabs(bfreq-efreq)>freqthresh || fabs(bfreq-mfreq)>freqthresh || fabs(efreq-mfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else if(gfreq>0 && efreq<0 && mfreq>0)
            {
                
                if(fabs(bfreq-gfreq)>freqthresh || fabs(bfreq-mfreq)>freqthresh || fabs(gfreq-mfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else if(gfreq>0 && efreq>0 && mfreq<0)
            {
                
                if(fabs(bfreq-gfreq)>freqthresh || fabs(bfreq-efreq)>freqthresh || fabs(gfreq-efreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else if(gfreq>0 && efreq<0 && mfreq<0)
            {
                
                if(fabs(bfreq-gfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else if(mfreq>0 && efreq<0 && gfreq<0)
            {
                
                if(fabs(bfreq-mfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else if(gfreq<0 && efreq>0 && mfreq<0)
            {
                
                if(fabs(bfreq-efreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else if (gfreq>0 && efreq>0 && mfreq>0)
            {
                if( fabs(bfreq-gfreq)>freqthresh || fabs(bfreq-efreq)>freqthresh || fabs(bfreq-mfreq)>freqthresh || fabs(gfreq-efreq)>freqthresh || fabs(gfreq-mfreq)>freqthresh || fabs(efreq-mfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=mdata->_esi_rs[mdata->_esi_include[i]] + '\t' + mdata->_esi_allele1[mdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\t' + atos(mfreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                    pasmid.push_back(mdata->_esi_include[i]);
                }
            }
            else
            {
                pasbid.push_back(bdata->_include[i]);
                pasgid.push_back(gdata->_include[i]);
                paseid.push_back(esdata->_esi_include[i]);
                pasmid.push_back(mdata->_esi_include[i]);
            }
        }
        if(paseid.size()==0)
        {
            printf("Error: No SNP passed Allele Frequency check. please check the frequecny colomn in your GWAS summmay data and your .esi file.\n");
            exit(EXIT_FAILURE);
        }
        bdata->_include = pasbid;
        gdata->_include = pasgid;
        esdata->_esi_include = paseid;
        mdata->_esi_include = pasmid;

        double prop=1.0*failcount/snpnum;
        if(prop > percenthresh) {
            success = false;
            printf("%ld SNPs (%0.2f%% > %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        } else  printf("%ld SNPs (%0.2f%% <= %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
       
        fclose(tmpfile);
        return prop;
    }

    double freq_check(bInfo* bdata, gwasData* gdata, gwasData* gdata2, double &freqthresh, double &percenthresh)
    {
        printf("Checking the consistency of allele frequency of each SNP between pairwise data sets (including the exposure GWAS summary data, the outcome GWAS summary data and the LD reference data).\n");
        bool success=true;
        long failcount=0, snpnum=bdata->_include.size();
        vector<int> pasbid, pasgid, paseid;
        string fname=string(outFileName) + ".snp_failed_freq_ck.list";
        FILE* tmpfile=fopen(fname.c_str(), "w");
        if (!(tmpfile)) {
            printf("Error: Failed to open file %s.\n",fname.c_str());
            exit(EXIT_FAILURE);
        }
        string tmpstr = "SNP\tAllele\tReferencePanel\texpGWAS\toutGWAS\n";
        fputs(tmpstr.c_str(),tmpfile);
        for(int i=0;i<bdata->_include.size();i++)
        {
            double bfreq=bdata->_mu[bdata->_include[i]]/2;
            double gfreq=gdata->freq[gdata->_include[i]];
            double efreq=gdata2->freq[gdata2->_include[i]];
            if(gfreq<0 && efreq>0)
            {
                if(fabs(bfreq-efreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=gdata2->snpName[gdata2->_include[i]]+ '\t' + gdata2->allele_1[gdata2->_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(gdata2->_include[i]);
                }
            }
            else if(gfreq>0 && efreq<0)
            {
                
                if(fabs(bfreq-gfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=gdata2->snpName[gdata2->_include[i]]+ '\t' + gdata2->allele_1[gdata2->_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(gdata2->_include[i]);
                }
            }
            else if (gfreq>0 && efreq>0)
            {
                if( fabs(bfreq-gfreq)>freqthresh || fabs(bfreq-efreq)>freqthresh || fabs(gfreq-efreq)>freqthresh )
                {
                    failcount++;
                    tmpstr=gdata2->snpName[gdata2->_include[i]]+ '\t' + gdata2->allele_1[gdata2->_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(gdata2->_include[i]);
                }
            }
            else
            {
                pasbid.push_back(bdata->_include[i]);
                pasgid.push_back(gdata->_include[i]);
                paseid.push_back(gdata2->_include[i]);
            }
        }
        if(paseid.size()==0)
        {
            printf("Error: No SNP passed Allele Frequency check. please check the frequecny colomn in your exposure and outcome GWAS summmay data.\n");
            exit(EXIT_FAILURE);
        }
        bdata->_include = pasbid;
        gdata->_include = pasgid;
        gdata2->_include = paseid;

        double prop=1.0*failcount/snpnum;
        if(prop > percenthresh) {
            success = false;
            printf("%ld SNPs (%0.2f%% > %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        } else  printf("%ld SNPs (%0.2f%% <= %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
       
        fclose(tmpfile);
        return prop;
    }

    void update_geIndx(bInfo* bdata, gwasData* gdata, gwasData* gdata2)
    {
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx1.push_back(gdata->_include[bdata->_include[i]]);
            tmpIdx2.push_back(gdata2->_include[bdata->_include[i]]);
        }
        gdata->_include.clear();
        gdata2->_include.clear();
        gdata->_include = tmpIdx1;
        gdata2->_include = tmpIdx2;
    }

    void update_geIndx(bInfo* bdata, eqtlInfo* mdata, eqtlInfo* esdata, gwasData* gdata)
    {
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        vector<int> tmpIdx3;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx1.push_back(mdata->_esi_include[bdata->_include[i]]);
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
            tmpIdx3.push_back(gdata->_include[bdata->_include[i]]);
        }
        mdata->_esi_include.clear();
        esdata->_esi_include.clear();
        gdata->_include.clear();
        mdata->_esi_include = tmpIdx1;
        esdata->_esi_include = tmpIdx2;
        gdata->_include = tmpIdx3;
    }

    void init_smr_wk(SMRWKM* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear(),smrwk->bmz.clear(),smrwk->semz.clear(),smrwk->zmz.clear();
        smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear(),smrwk->freq.clear(),smrwk->zxz.clear();
    }
    void init_med_info(MED_info* med_info)
    {
        med_info->Med_ID.clear(), med_info->Gene_ID.clear(), med_info->Med_Chr.clear(), med_info->Med_bp.clear(), med_info->b_xm.clear(), med_info->se_xm.clear(), med_info->p_xm.clear();
        med_info->n_xm.clear(), med_info->b_my.clear(), med_info->se_my.clear(), med_info->p_my.clear(), med_info->n_my.clear();
    }
    void init_smr_wk(SMRWKMULT* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear();
        smrwk->Xm_b.resize(0,0), smrwk->Xm_se.resize(0,0);
        smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear(),smrwk->freq.clear(),smrwk->zxz.clear();
    }
    void init_smr_wk(SMRWKMULTTRANS* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->bmz.clear(),smrwk->semz.clear(),smrwk->seyz.clear(),smrwk->pyz.clear();
        smrwk->Xm_b.resize(0,0), smrwk->Xm_se.resize(0,0), smrwk->Xm_z.resize(0,0);
        smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear(),smrwk->freq.clear(),smrwk->zxz.clear(),smrwk->zmz.clear();
    }
    void fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, gwasData* meddata, SMRWKM* smrwk,int cis_itvl)
    {
        int i=smrwk->cur_prbidx;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_sexz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==esdata->_epi_chr[i] && fabs(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[j]+9>1e-6 && meddata->seyz[j]+9>1e-6)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->bmz.push_back(meddata->byz[j]);
                        smrwk->semz.push_back(meddata->seyz[j]);
                        smrwk->zmz.push_back(meddata->byz[j]/meddata->seyz[j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                    }
                }
            }   
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=esdata->_rowid[beta_start+j];
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                if(snpchr==esdata->_epi_chr[i] && abs(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[ge_rowid]+9>1e-6 && meddata->seyz[ge_rowid]+9>1e-6)
                {
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    smrwk->bmz.push_back(meddata->byz[ge_rowid]);
                    smrwk->semz.push_back(meddata->seyz[ge_rowid]);
                    smrwk->zmz.push_back(meddata->byz[ge_rowid]/meddata->seyz[ge_rowid]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                }
            }
        }
    }

    void fill_smr_wk_trans(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, gwasData* meddata, SMRWKM* smrwk)
    {
        int i=smrwk->cur_prbidx;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_sexz[i][j] + 9) > 1e-6)
                {
                    if(gdata->seyz[j]+9>1e-6 && meddata->seyz[j]+9>1e-6)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->bmz.push_back(meddata->byz[j]);
                        smrwk->semz.push_back(meddata->seyz[j]);
                        smrwk->zmz.push_back(meddata->byz[j]/meddata->seyz[j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                    }
                }
            }   
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=esdata->_rowid[beta_start+j];
                if(gdata->seyz[ge_rowid]+9>1e-6 && meddata->seyz[ge_rowid]+9>1e-6)
                {
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    smrwk->bmz.push_back(meddata->byz[ge_rowid]);
                    smrwk->semz.push_back(meddata->seyz[ge_rowid]);
                    smrwk->zmz.push_back(meddata->byz[ge_rowid]/meddata->seyz[ge_rowid]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                }
            }
        }
    }

    long fill_smr_wk_trans(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_sexz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(gdata->seyz[j]+9>1e-6)
                    {
                        if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                            {
                                smrwk->bxz.push_back(esdata->_bxz[i][j]);
                                smrwk->sexz.push_back(esdata->_sexz[i][j]);
                                smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                smrwk->byz.push_back(gdata->byz[j]);
                                smrwk->seyz.push_back(gdata->seyz[j]);
                                smrwk->pyz.push_back(gdata->pvalue[j]);
                                smrwk->curId.push_back(j);
                                smrwk->rs.push_back(esdata->_esi_rs[j]);
                                smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                                smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                                smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                                smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                                smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata->_esi_rs[j])+"\t"+ atos(esdata->_esi_chr[j])+"\t"+ atos(esdata->_esi_bp[j])+"\t"+ atos(esdata->_esi_allele1[j])+"\t"+ atos(esdata->_esi_allele2[j])+"\t"+ atos(esdata->_esi_freq[j])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][j])+"\t"+ atos(esdata->_sexz[i][j])+"\t"+ dtos(p)+"\n";
                        
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            smrwk->bxz.push_back(esdata->_bxz[i][j]);
                            smrwk->sexz.push_back(esdata->_sexz[i][j]);
                            smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                            smrwk->byz.push_back(gdata->byz[j]);
                            smrwk->seyz.push_back(gdata->seyz[j]);
                            smrwk->pyz.push_back(gdata->pvalue[j]);
                            smrwk->curId.push_back(j);
                            smrwk->rs.push_back(esdata->_esi_rs[j]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                            smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.

                        }
                    }
                }
            }
            
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=esdata->_rowid[beta_start+j];
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                if(gdata->seyz[ge_rowid]+9>1e-6)
                {
                    if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0)
                    {
                        if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                        {
                            smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                            smrwk->sexz.push_back(esdata->_val[se_start+j]);
                            smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                            smrwk->byz.push_back(gdata->byz[ge_rowid]);
                            smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                            smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                            smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                            smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                            smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                        } else {
                            printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                            double z=(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                            double p=pchisq(z*z, 1);
                            string tmp=atos(esdata->_esi_rs[j])+"\t"+ atos(esdata->_esi_chr[j])+"\t"+ atos(esdata->_esi_bp[j])+"\t"+ atos(esdata->_esi_allele1[j])+"\t"+ atos(esdata->_esi_allele2[j])+"\t"+ atos(esdata->_esi_freq[j])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][j])+"\t"+ atos(esdata->_sexz[i][j])+"\t"+ dtos(p)+"\n";
                            
                            printf("%s\n",tmp.c_str());
                        }
                        
                    } else {
                        smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                        smrwk->sexz.push_back(esdata->_val[se_start+j]);
                        smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                        smrwk->byz.push_back(gdata->byz[ge_rowid]);
                        smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                        smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                        smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                        smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                    }
                  
                }
            }
        }
        
        return maxid;
    }

    void sbat_calcu_lambda(MatrixXd &X, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx,vector<double> &zxz4smr)
    {
        int m = snp_count;
       
        vector<int> rm_ID1;
        double R_cutoff = sbat_ld_cutoff;
        int qi = 0; //alternate index
        
        MatrixXd C;
        cor_calc(C, X);
        if (sbat_ld_cutoff < 1) rm_cor_sbat(C, R_cutoff, m, rm_ID1,zxz4smr);
        //Create new index
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sub_indx.push_back(i);
            else {
                if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                else sub_indx.push_back(i);
            }
        }
        snp_count = (int)sub_indx.size();
        if (sub_indx.size() < C.size()) { //Build new matrix
            MatrixXd D(sub_indx.size(),sub_indx.size());
            for (int i = 0 ; i < sub_indx.size() ; i++) {
                for (int j = 0 ; j < sub_indx.size() ; j++) {
                    D(i,j) = C(sub_indx[i],sub_indx[j]);
                }
            }
            C = D;
        }    
    }

    vector<int> sort_indexes(vector<double> &v) {

    // adapted from Łukasz Wiklendt, Stack Overflow
    // initialize original index locations
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
    }

    double qnorm5(double p, int lower_tail)
    {
        double p_, q, r, val;
        int log_p = 0;

        p_ = p;/* real lower_tail prob. p */
        q = p_ - 0.5;

    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

            Produces the normal deviate Z corresponding to a given lower
            tail area of P; Z is accurate to about 1 part in 10**16.

            (original fortran code used PARAMETER(..) for the coefficients
            and provided hash codes for checking them...)
    */
        if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
            r = .180625 - q * q;
        val =
                q * (((((((r * 2509.0809287301226727 +
                        33430.575583588128105) * r + 67265.770927008700853) * r +
                        45921.953931549871457) * r + 13731.693765509461125) * r +
                    1971.5909503065514427) * r + 133.14166789178437745) * r +
                    3.387132872796366608)
                / (((((((r * 5226.495278852854561 +
                        28729.085735721942674) * r + 39307.89580009271061) * r +
                    21213.794301586595867) * r + 5394.1960214247511077) * r +
                    687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
        }
        else { /* closer than 0.075 from {0,1} boundary */

        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = 1-p;/* 1-p */
        else
            r = p_;/* = R_DT_Iv(p) ^=  p */

        r = sqrt(- ((log_p &&
                ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
                p : /* else */ log(r)));
            /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

            if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
                r += -1.6;
                val = (((((((r * 7.7454501427834140764e-4 +
                        .0227238449892691845833) * r + .24178072517745061177) *
                        r + 1.27045825245236838258) * r +
                        3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                    1.42343711074968357734)
                    / (((((((r *
                            1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                            r + .0151986665636164571966) * r +
                        .14810397642748007459) * r + .68976733498510000455) *
                        r + 1.6763848301838038494) * r +
                        2.05319162663775882187) * r + 1.);
            }
            else { /* very close to  0 or 1 */
                r += -5.;
                val = (((((((r * 2.01033439929228813265e-7 +
                        2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                        r + .29656057182850489123) * r +
                    1.7848265399172913358) * r + 5.4637849111641143699) *
                    r + 6.6579046435011037772)
                    / (((((((r *
                            2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                            r + 1.8463183175100546818e-5) * r +
                        7.868691311456132591e-4) * r + .0148753612908506148525)
                        * r + .13692988092273580531) * r +
                        .59983220655588793769) * r + 1.);
            }

        if(q < 0.0)
            val = -val;
            /* return (q >= 0.)? r : -r ;*/
        }
        return val;
    }

    int min_id(vector<double> &pyz)
    {
        int id=0;
        double tmpVal, cmpVal=abs(pyz[0]);
        for( int i=1;i<pyz.size();i++)
        {
            tmpVal=abs(pyz[i]);
            if( cmpVal-tmpVal > 1e-10)
            {
                cmpVal=tmpVal;
                id=i;
            }
        }
        return(id);
    }

    int min_id(VectorXd &pyz)
    {
        int id=0;
        double tmpVal, cmpVal=pyz[0];
        for( int i=1;i<pyz.size();i++)
        {        
            tmpVal=pyz[i];
            if( cmpVal-tmpVal > 1e-50)
            {
                cmpVal=tmpVal;
                id=i;
            }
        }
        return(id);
    }

    void mr_ivw(VectorXd &bzx, VectorXd &bzy, VectorXd &sezx, VectorXd &sezy, double &beta_ivw, double &se_ivw, double &p_ivw)
    {
        int m = bzx.size();
        
        VectorXd varxy(m), one(m), res(m), epsilon(m);
        double sigma_2, z_ivw;

        beta_ivw = (bzy.dot(bzx))/(bzx.dot(bzx));

        if (m < 2){
            //se_ivw = sqrt(pow(sezy[0], 2)/pow(bzx[0],2) + pow(bzy[0]*sezx[0],2)/pow(bzx[0],4));
            se_ivw = sqrt(pow(sezy[0], 2)/pow(bzx[0],2));
        } else {
            res = (bzy - bzx*beta_ivw);
            sigma_2 = (res.dot(res))/(m-1);
            for (int j=0; j<m; j++){
                epsilon[j] = max(sigma_2, pow(sezy[j],2));
                //varxy[j] = epsilon[j]/pow(bzx[j],2) + pow(bzy[j]*sezx[j],2)/pow(bzx[j],4);
                varxy[j] = epsilon[j]/pow(bzx[j],2);
                one[j] = 1;
            }
            se_ivw = sqrt(1/(one.array()/varxy.array()).sum());
        }
        z_ivw = beta_ivw/se_ivw;
        p_ivw = 2*pnorm(fabs(z_ivw));
    }

    void mr_ivw_LD(VectorXd &bzx, VectorXd &bzy, VectorXd &sezx, VectorXd &sezy, MatrixXd &C, double &beta_ivw, double &se_ivw, double &p_ivw, int N_ref)
    {
        int m = bzx.size();
        MatrixXd W(m,m), D(m,m), Sigma(m,m), I_mat, D_sqrt_inv(m,m);
        VectorXd epsilon(m), res(m);
        double sigma_2, z_ivw;
        I_mat.setIdentity(m, m);

        D =  (1-2/sqrt(N_ref))*C+(2/sqrt(N_ref))*I_mat;
        W = D.inverse();

        SelfAdjointEigenSolver<MatrixXd> es(D);
        D_sqrt_inv = es.operatorInverseSqrt();

        beta_ivw = (bzx.transpose()*W*bzx).inverse()*bzx.transpose()*W*bzy;

        if (m < 2){
            se_ivw = sqrt(pow(sezy[0], 2)/pow(bzx[0],2));
        } else {
            res = (bzy - bzx*beta_ivw);
            sigma_2 = res.dot(W*res)/(m-1);
            for (int j=0; j<m; j++){
                epsilon[j] = max(sigma_2, pow(sezy[j],2));
            }
            Sigma = epsilon.asDiagonal();
            se_ivw = sqrt((bzx.transpose()*W*bzx).inverse()*bzx.dot(D_sqrt_inv*Sigma*D_sqrt_inv*bzx)*(bzx.transpose()*W*bzx).inverse());
        }
        z_ivw = beta_ivw/se_ivw;
        p_ivw = 2*pnorm(fabs(z_ivw));
    }

    void mvmr_ivw_LD(MatrixXd &X, VectorXd &y, VectorXd &SEs, MatrixXd &C, VectorXd &betas, VectorXd &vars, int N_ref)
    {
        int m = X.rows();
        int k = X.cols();
        
        MatrixXd D(m,m), W(m,m), M(k,k), J(k, (k+1)*m), df_dg(k,m), df_dG(k, k*m), Sigma((k+1)*m, (k+1)*m), I_mat, R;
        VectorXd res(m);
        I_mat.setIdentity(m, m);
        R.setIdentity(k+1, k+1);
        double sigma_2;

        betas.resize(k), vars.resize(k);

        D =  (1-2/sqrt(N_ref))*C+(2/sqrt(N_ref))*I_mat;

        // Calculate the MR estimates
        W = D.inverse();
        M = (X.transpose()*W*X).inverse();
        betas = M*(X.transpose()*W*y);

        // Calculate the standard errors
        df_dg = M*X.transpose()*W;
        df_dG = kroneckerProduct(M, y.transpose()*W*(-(X*M*X.transpose())*W + I_mat)) + kroneckerProduct(-y.transpose()*W*X*M, M*X.transpose()*W);

        J << df_dG, df_dg; 

        // outcome error is the maximum between se^2_zy and residual
        res = (y - X*betas);
        sigma_2 = (res.dot(res))/(m-k-1);
        double epsilon;

        for (int j=0; j<m; j++){
            epsilon = max(sigma_2, pow(SEs[m*k+j],2));
            SEs[m*k+j] = sqrt(epsilon);
        }

        // standard errors              
        Sigma = (SEs*SEs.transpose()).array() * kroneckerProduct(R, D).array();
        vars = (J*Sigma*J.transpose()).diagonal();        
    }

    int smr_ivw_test(bInfo* bdata, vector<uint32_t> &slctId, vector<string> &slct_snpName, vector<string> &slct_a1, vector<string> &slct_a2, vector<double> &slct_bxz,vector<double> &slct_sexz, vector<double> &slct_pxz, vector<double> &slct_byz,vector<double> &slct_seyz, double &beta_ivw, double &se_ivw, double &p_ivw, double &h2cis, double p_smr, double ld_top_multi, bool ldmatrix, int N_ref, double trev, vector<string> &snp4msmr, bool p_expo_provided, bool get_snp_effects_flg, int min_snp, FILE* snpfile, string snpfilename, string exponame)
    {
        /* step4: Filter out the SNPs with p-smr threshold and ld-pruning */
        printf("Conducting IVW multi-SNP SMR test...\n");
        
        vector<uint32_t> Id4smr;
        vector<double> bxz4smr, sexz4smr, byz4smr, seyz4smr, zxz4smr, zyz4smr;
        vector<string> snpname4mr, allele14mr, allele24mr;
        double z_smr= fabs(qnorm5(p_smr/2));
        
        for(int j=0;j<slctId.size();j++)
        {
            double ztmp;

            if (p_expo_provided) {
                ztmp = fabs(qnorm5(slct_pxz[j]/2));
            } else {
                ztmp=fabs(slct_bxz[j]/slct_sexz[j]);
            }
            
            double zrev = (fabs(slct_bxz[j]) - fabs(slct_byz[j]))/sqrt(pow(slct_sexz[j],2) + pow(slct_seyz[j],2)); 
            if((ztmp>=z_smr) && (zrev > trev))
            {
                Id4smr.push_back(slctId[j]);
                snpname4mr.push_back(slct_snpName[j]);
                allele14mr.push_back(slct_a1[j]);
                allele24mr.push_back(slct_a2[j]);
                bxz4smr.push_back(slct_bxz[j]);
                sexz4smr.push_back(slct_sexz[j]);
                zxz4smr.push_back(ztmp);
                byz4smr.push_back(slct_byz[j]);
                seyz4smr.push_back(slct_seyz[j]);
                zyz4smr.push_back(slct_byz[j]/slct_seyz[j]);
            }
        }
        int snp_count=(int)Id4smr.size();
        printf("%ld SNPs passed the p-value threshold %6.2e and %ld SNPs are excluded.\n",Id4smr.size(), p_smr,slctId.size()-Id4smr.size());
        if(snp_count==0) return -9;
        /* step5: multiple-SNP SMR test */
    
        vector<int> sub_indx;
        MatrixXd _X;
        make_XMat(bdata, Id4smr, _X); //_X: one row one individual, one column one SNP

        double sbat_ld_cutoff=sqrt(ld_top_multi);
        sbat_calcu_lambda(_X, snp_count,  sbat_ld_cutoff, sub_indx, zxz4smr); //the index of slectId, snp_count can change here
        printf("%ld SNPs passed LD-square threshold of %6.2f and %ld SNPs are excluded.\n",sub_indx.size(), ld_top_multi, Id4smr.size()-sub_indx.size());
        // IVW calculation
        printf("%ld SNPs are included in the IVW multi-SNP SMR test.\n",sub_indx.size());
    
        VectorXd bzy(sub_indx.size()), bzx(sub_indx.size()), sezx(sub_indx.size()), sezy(sub_indx.size());
        vector<uint32_t> Id4mrdelta;
        h2cis = 0;
        string snpstr;

        for(int j=0; j<sub_indx.size();j++){
            bzy[j] = byz4smr[sub_indx[j]];
            bzx[j] = bxz4smr[sub_indx[j]];
            h2cis+= bxz4smr[sub_indx[j]]*bxz4smr[sub_indx[j]] - sexz4smr[sub_indx[j]]*sexz4smr[sub_indx[j]];
            sezx[j] = sexz4smr[sub_indx[j]];
            sezy[j] = seyz4smr[sub_indx[j]];
            Id4mrdelta.push_back(Id4smr[sub_indx[j]]);

            if (get_snp_effects_flg && (snp_count >= min_snp)){
                snpstr = exponame + '\t' + snpname4mr[sub_indx[j]] + '\t' + allele14mr[sub_indx[j]] + '\t' + allele24mr[sub_indx[j]] + '\t';
                snpstr += atos(bxz4smr[sub_indx[j]]) + '\t' + atos(sexz4smr[sub_indx[j]]) + '\t';
                snpstr += atos(byz4smr[sub_indx[j]]) + '\t' + atos(seyz4smr[sub_indx[j]]) + '\n'; 

                if(fputs_checked(snpstr.c_str(),snpfile))
                {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
                }
            }     
        }

        if (!ldmatrix) mr_ivw(bzx, bzy, sezx, sezy, beta_ivw, se_ivw, p_ivw);
        else {
            MatrixXd C, _X;
            make_XMat(bdata, Id4mrdelta, _X);
            cor_calc(C, _X);     
            mr_ivw_LD(bzx, bzy, sezx, sezy, C, beta_ivw, se_ivw, p_ivw, N_ref);            
        }

        cout << "MR effects: b = " << beta_ivw << "; p = " << p_ivw << endl;

        for (int j = 0; j < sub_indx.size(); j++) snp4msmr.push_back(bdata->_snp_name[bdata->_include[Id4mrdelta[j]]]);
        return snp_count;
    }

    void smr_ivw_gwas_gwas(char* outFileName, char* bFileName,char* gwasFileName, char* gwasFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp)
    {
        double theta=0;
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        gwasData gdata2;
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(gwasFileName2==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        read_gwas_data( &gdata, gwasFileName);
        read_gwas_data( &gdata2, gwasFileName2);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &gdata, &gdata2);
        
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &gdata, &gdata2);
        }
        if(forcefrqck)
        {
            double prop= freq_check(&bdata, &gdata, &gdata2,afthresh,percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
        update_gwas(&gdata);
        update_gwas(&gdata2);
        
        
        cout<<endl<<"Performing IVW SMR analysis..... "<<endl;
           
        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }
        
        string smrfile = string(outFileName)+".msmr";
        FILE* smr=NULL;
        smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("Open error %s\n", smrfile.c_str());
            exit(1);
        }
        
        string outstr="ExposureGWAS\tOutcomeGWAS\tb_ivw\tse_ivw\tp_ivw\tn_iv\th2cis\n";
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }

        FILE* snpfile;
        string snpfilename = string(outFileName) + ".snps";
        string snpstr;
        if (get_snp_effects_flg){
            snpfile = fopen(snpfilename.c_str(), "w");
            if (!(snpfile)) {
                printf("ERROR: open error %s\n", snpfilename.c_str());
                exit(1);
            }
            snpstr="ExposureGWAS\tSNP\tA1\tA2\tbeta.exp\tse.exp\tbeta.out\tse.out\n";
            if(fputs_checked(snpstr.c_str(),snpfile))
            {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }

        vector<uint32_t> slctId;
        vector<int> slct_bpsnp,slct_snpchr;
        vector<double> slct_bxz, slct_sexz, slct_byz, slct_seyz, slct_pxz; //slct_zsxz,slct_zxz would be removed one of them
        vector<string> slct_snpName, slct_a1, slct_a2;

        for (int j = 0; j<bdata._include.size(); j++){
            
            slct_snpName.push_back(gdata.snpName[j]);
            slct_a1.push_back(gdata.allele_1[j]);
            slct_a2.push_back(gdata.allele_2[j]);
            slctId.push_back(j);
            slct_bxz.push_back(gdata.byz[j]);
            slct_sexz.push_back(gdata.seyz[j]);
            slct_byz.push_back(gdata2.byz[j]);
            slct_seyz.push_back(gdata2.seyz[j]);
            slct_pxz.push_back(gdata.pvalue[j]);

        } 
            
        vector<string> snp4msmr;
        double beta_ivw, se_ivw, p_ivw, h2cis;
        beta_ivw = -9;
        se_ivw = -9;
        h2cis = -9;
        bool p_expo_provided = true;
        int snp_count=smr_ivw_test(&bdata, slctId, slct_snpName, slct_a1, slct_a2, slct_bxz,slct_sexz,slct_pxz, slct_byz,slct_seyz, beta_ivw, se_ivw, p_ivw, h2cis, p_smr, ld_top_multi, ldmatrix, N_ref, trev, snp4msmr, p_expo_provided, get_snp_effects_flg, min_snp, snpfile, snpfilename, gwasFileName);

        // output snp set list
        string setstr="SNP_instruments\n";
        if(fputs_checked(setstr.c_str(),setlst))
        {
            printf("ERROR: in writing file %s .\n", setlstfile.c_str());
            exit(EXIT_FAILURE);
        }
        for(int j=0;j<snp4msmr.size();j++)
        {
            setstr=snp4msmr[j]+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
        }
        // end of output

        outstr = std::string(gwasFileName) + '\t' + gwasFileName2 + atos(beta_ivw) + '\t' + atos(se_ivw) + '\t' + dtos(p_ivw) + '\t' + atos(snp_count) + '\t' + atos(h2cis) + '\n';
        
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }

        cout<<"\nGWAS to GWAS MR analysis has been saved in the file [" + smrfile + "]."<<endl;
        fclose(smr);
        fclose(setlst);
        if (get_snp_effects_flg){
            fclose(snpfile);
        }
    }

    void smr_ivw_analysis(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, int cis_itvl, char* genelistName, int chr, int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp)
    {
        double theta=0;
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, prbname);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &gdata, &esdata);
        
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &gdata, &esdata);
        }
        if(forcefrqck)
        {
            double prop= freq_check(&bdata, &gdata, &esdata,afthresh,percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
        update_gwas(&gdata);
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName)+".epi");
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        
        vector<string> set_name;
        vector< vector<string> > snpset;
        vector<int> gene_chr,gene_bp1,gene_bp2;
        if(setlstName!=NULL) sbat_read_snpset(&bdata,setlstName,set_name,gene_chr, gene_bp1,gene_bp2, snpset );
        else if(geneAnnoFileName!=NULL) read_geneAnno(geneAnnoFileName, set_name, gene_chr, gene_bp1, gene_bp2);
             
        unsigned int probNum = esdata._probNum;
               
        cout<<endl<<"Performing MR-IVW analyses..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);

        cis_itvl=cis_itvl*1000;
        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }
        
        string smrfile = string(outFileName)+".msmr";
        FILE* smr=NULL;
        smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("Open error %s\n", smrfile.c_str());
            exit(1);
        }
         
        string outstr="ProbeID\tProbeChr\tProbeName\tProbe_bp\ttopSNP\tA1\tA2\ttopSNP_chr\ttopSNP_bp\tb_top\tse_top\tp_top\tb_ivw\tse_ivw\tp_ivw\tn_iv\th2cis\n";
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }

        FILE* snpfile;
        string snpfilename = string(outFileName) + ".snps";
        string snpstr;
        if (get_snp_effects_flg){
            snpfile = fopen(snpfilename.c_str(), "w");
            if (!(snpfile)) {
                printf("ERROR: open error %s\n", snpfilename.c_str());
                exit(1);
            }
            snpstr="ProbeID\tSNP\tA1\tA2\tbeta.exp\tse.exp\tbeta.out\tse.out\n";
            if(fputs_checked(snpstr.c_str(),snpfile))
            {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }
        long write_count=0;
        map<string, int>::iterator iter;
        SMRWK smrwk;
        
        for(int i=0;i<probNum;i++)
        {
            
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
            int probebp=esdata._epi_bp[i];
            int probechr=esdata._epi_chr[i];
            string probename=esdata._epi_prbID[i];
            string probegene=esdata._epi_gene[i];
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=i;
            // step1: get cis-eQTLs
            printf("\nInitiating the workspace of probe %s for MR-IVW analysis....\n",probename.c_str());
            char* refSNP=NULL;
            bool Flag = false;
            long maxid =fill_smr_wk(&bdata, &gdata, &esdata, &smrwk, refSNP, cis_itvl, Flag);
            if (smrwk.bxz.size() == 0) {
                printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                continue;
            }
            printf("%ld SNPs are included from the cis-region of the probe %s.\n",smrwk.bxz.size(),probename.c_str());
            // step2: get top-SNP
            printf("Checking the top-SNP in the region....\n");
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            VectorXd zsxz;
            zsxz=ei_bxz.array()/ei_sexz.array();
            maxid=max_abs_id(zsxz); // now maxid point to the sig eQTL SNP or ref SNP in the new datastruct(not the raw).
            double pxz_val = pnorm(fabs(zsxz[maxid]));
            //double computing, consistency should be checked
            for(int tid=0;tid<zsxz.size();tid++) {
                if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                {
                    printf("ERROR: zxz not consistent %f vs %f. please report this.\n",zsxz(tid),smrwk.zxz[tid]);
                    exit(EXIT_FAILURE);
                }
            }
            string topsnpname=smrwk.rs[maxid];
            printf("The top SNP of probe %s is %s with p-value %e.\n", probename.c_str(), topsnpname.c_str(),pxz_val);
            if(pxz_val>p_smr){
                printf("WARNING: no SNP passed the p-value threshold %e for MR-IVW analysis for probe %s.\n", p_smr, probename.c_str());
                continue;
            }
            // step3: extract SNPs around the --set-wind around sig SNP
            cout << "Extracting SNPs in the cis-region..." << endl;
            vector<uint32_t> slctId;
            vector<int> slct_bpsnp,slct_snpchr;
            vector<double> slct_bxz, slct_sexz, slct_pxz, slct_byz, slct_seyz, slct_zsxz,slct_zxz, slct_pyz,slct_freq; //slct_zsxz,slct_zxz would be removed one of them
            vector<string> slct_snpName, slct_a1, slct_a2;
            long slct_maxid=-9;
            slctId.swap(smrwk.curId);
            slct_bxz.swap(smrwk.bxz);
            slct_sexz.swap(smrwk.sexz);
            slct_byz.swap(smrwk.byz);
            slct_seyz.swap(smrwk.seyz);
            slct_snpName.swap(smrwk.rs);
            slct_maxid=maxid;
            for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
            slct_zxz.swap(smrwk.zxz);
            slct_pyz.swap(smrwk.pyz);
            slct_bpsnp.swap(smrwk.bpsnp);
            slct_snpchr.swap(smrwk.snpchrom);
            slct_a1.swap(smrwk.allele1);
            slct_a2.swap(smrwk.allele2);
            slct_freq.swap(smrwk.freq);         
            
            int out_raw_id = slctId[slct_maxid];
            double bxz_max = slct_bxz[slct_maxid];
            double sexz_max = slct_sexz[slct_maxid];
            double byz_max = slct_byz[slct_maxid];
            double seyz_max = slct_seyz[slct_maxid];
            double bxy_max = byz_max / bxz_max;
            //double sexy_max = sqrt(pow(seyz_max, 2)/pow(bxz_max,2) + pow(byz_max*sexz_max,2)/pow(bxz_max,4));
            double sexy_max = sqrt(pow(seyz_max, 2)/pow(bxz_max,2));
            double zxy_max = bxy_max / sexy_max;
            double pxy_max = 2*pnorm(fabs(zxy_max));
            
            vector<string> snp4msmr;
            double beta_ivw, se_ivw, p_ivw, h2cis;
            beta_ivw = -9;
            se_ivw = -9;
            h2cis = -9;
            bool p_expo_provided = false;
            int snp_count=smr_ivw_test(&bdata, slctId, slct_snpName, slct_a1, slct_a2, slct_bxz, slct_sexz,slct_pxz, slct_byz, slct_seyz, beta_ivw, se_ivw, p_ivw, h2cis, p_smr, ld_top_multi, ldmatrix, N_ref, trev, snp4msmr,p_expo_provided, get_snp_effects_flg, min_snp, snpfile, snpfilename, probename);
            if(snp_count==-9) continue;
            // output snp set list
            string setstr=probename+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<snp4msmr.size();j++)
            {
                setstr=snp4msmr[j]+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            setstr="end\n";
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            // end of output
                      
            outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + esdata._esi_allele1[out_raw_id] + '\t' + esdata._esi_allele2[out_raw_id] + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t';
            outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t';
            outstr += atos(beta_ivw) + '\t' + atos(se_ivw) + '\t' + dtos(p_ivw) + '\t' + atos(snp_count) + '\t' + atos(h2cis) + '\n';
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
            write_count++;

        }

        cout<<"\nMR-IVW analyses completed.\nAnalysis results of "<<write_count<<" sets have been saved in the file [" + smrfile + "]."<<endl;
        cout<<"SNP sets included in the MR-IVW analyses have been saved in the file [" + setlstfile + "]."<<endl;
        fclose(smr);
        fclose(setlst);
        if (get_snp_effects_flg){
            fclose(snpfile);
        }
    }

void smr_rev_ivw_analysis(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_gwas, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp)
    {
        double theta=0;
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        read_gwas_data(&gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, prbname);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &gdata, &esdata);
        
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &gdata, &esdata);
        }
        if(forcefrqck)
        {
            double prop= freq_check(&bdata, &gdata, &esdata,afthresh,percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
        update_gwas(&gdata);
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName)+".epi");
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        
        vector<string> set_name;
        vector< vector<string> > snpset;
        vector<int> gene_chr,gene_bp1,gene_bp2;
        if(setlstName!=NULL) sbat_read_snpset(&bdata,setlstName,set_name,gene_chr, gene_bp1,gene_bp2, snpset );
        else if(geneAnnoFileName!=NULL) read_geneAnno(geneAnnoFileName, set_name, gene_chr, gene_bp1, gene_bp2);
        
      
        unsigned int probNum = esdata._probNum;
        
        
        cout<<endl<<"Performing multi-SNP based SMR analysis..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);

        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }
        
        string smrfile = string(outFileName)+".msmr";
        FILE* smr=NULL;
        smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("Open error %s\n", smrfile.c_str());
            exit(1);
        }
       
        string outstr="ProbeID\tProbeChr\tProbeName\tProbe_bp\ttopSNP\tA1\tA2\ttopSNP_chr\ttopSNP_bp\tb_top\tse_top\tp_top\tb_ivw\tse_ivw\tp_ivw\tn_iv\th2cis\n";
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }

        FILE* snpfile;
        string snpfilename = string(outFileName) + ".snps";
        string snpstr;
        if (get_snp_effects_flg){
            snpfile = fopen(snpfilename.c_str(), "w");
            if (!(snpfile)) {
                printf("ERROR: open error %s\n", snpfilename.c_str());
                exit(1);
            }
            snpstr="ProbeID\tSNP\tA1\tA2\tbeta.exp\tse.exp\tbeta.out\tse.out\n";
            if(fputs_checked(snpstr.c_str(),snpfile))
            {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }

        long write_count=0;
        map<string, int>::iterator iter;
        SMRWK smrwk;

        for(int i=0;i<probNum;i++)
        {
            
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
            int probebp=esdata._epi_bp[i];
            int probechr=esdata._epi_chr[i];
            string probename=esdata._epi_prbID[i];
            string probegene=esdata._epi_gene[i];
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=i;
            // step1: get cis-eQTLs
            printf("\nInitiating the workspace of probe %s for IVW reverse multi-SNP SMR analysis....\n",probename.c_str());
            long maxid =fill_smr_wk_trans(&bdata, &gdata, &esdata, &smrwk);
            if (smrwk.bxz.size() == 0) {
                printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                continue;
            }
            printf("%ld SNPs are included from the cis-region of the probe %s.\n",smrwk.bxz.size(),probename.c_str());
            // step2: get top-SNP
            printf("Checking the top-SNP in the region...\n");
            Map<VectorXd> ei_pyz(&smrwk.pyz[0],smrwk.pyz.size());
            VectorXd pyz=ei_pyz.array();
            maxid=min_id(pyz); 
            double pyz_val = pyz[maxid];
            
            string topsnpname=smrwk.rs[maxid];
            printf("The top SNP of probe %s is %s with p-value %e.\n", probename.c_str(), topsnpname.c_str(),pyz_val);
            if(pyz_val> p_gwas){
                printf("WARNING: no SNP passed the p-value threshold %e for IVW reverse multi-SNP SMR for probe %s.\n", p_gwas, probename.c_str());
                continue;
            }
            // step3: extract SNPs 
            //else 
            printf("Extracting SNPs ....\n");
            vector<uint32_t> slctId;
            vector<int> slct_bpsnp,slct_snpchr;
            vector<double> slct_bxz, slct_sexz, slct_byz, slct_seyz, slct_pyz,slct_freq; //slct_zsxz,slct_zxz would be removed one of them
            vector<string> slct_snpName,slct_a1, slct_a2;
            long slct_maxid=-9;

            slctId.swap(smrwk.curId);
            slct_bxz.swap(smrwk.bxz);
            slct_sexz.swap(smrwk.sexz);
            slct_byz.swap(smrwk.byz);
            slct_seyz.swap(smrwk.seyz);
            slct_snpName.swap(smrwk.rs);
            slct_maxid=maxid;
            slct_pyz.swap(smrwk.pyz);
            slct_bpsnp.swap(smrwk.bpsnp);
            slct_snpchr.swap(smrwk.snpchrom);
            slct_a1.swap(smrwk.allele1);
            slct_a2.swap(smrwk.allele2);
            slct_freq.swap(smrwk.freq);
            
            int out_raw_id = slctId[slct_maxid];
            double bxz_max = slct_bxz[slct_maxid];
            double sexz_max = slct_sexz[slct_maxid];
            double byz_max = slct_byz[slct_maxid];
            double seyz_max = slct_seyz[slct_maxid];
            // reverse formulas
            double bxy_max =  bxz_max/ byz_max;
            //double sexy_max = sqrt(pow(sexz_max, 2)/pow(byz_max,2) + pow(bxz_max*seyz_max,2)/pow(byz_max,4));
            double sexy_max = sqrt(pow(sexz_max, 2)/pow(byz_max,2));
            double zxy_max = bxy_max / sexy_max;
            double pxy_max = 2*pnorm(fabs(zxy_max));
            
            vector<string> snp4msmr;
            double beta_ivw, se_ivw, p_ivw, h2cis;
            beta_ivw = -9;
            se_ivw = -9;
            h2cis = -9;
            bool p_expo_provided = true;
            // reverse input
            int snp_count=smr_ivw_test(&bdata, slctId, slct_snpName, slct_a1, slct_a2, slct_byz,slct_seyz,slct_pyz, slct_bxz,slct_sexz, beta_ivw, se_ivw, p_ivw, h2cis, p_gwas, ld_top_multi, ldmatrix, N_ref, trev, snp4msmr, p_expo_provided, get_snp_effects_flg, min_snp, snpfile, snpfilename, probename);
            if(snp_count==-9) continue;
            // output snp set list
            string setstr=probename+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<snp4msmr.size();j++)
            {
                setstr=snp4msmr[j]+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            setstr="end\n";
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            // end of output
            outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + esdata._esi_allele1[out_raw_id] + '\t' + esdata._esi_allele2[out_raw_id] + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t';
            outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t';
            outstr += atos(beta_ivw) + '\t' + atos(se_ivw) + '\t' + dtos(p_ivw) + '\t' + atos(snp_count) + '\t' + atos(h2cis) + '\n';
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
            write_count++;
            
        }   

        cout<<"\nReverse MR-IVW completed.\nAnalysis results of "<<write_count<<" probes have been saved in the file [" + smrfile + "]."<<endl;
        cout<<"SNP sets included in reverse MR-IVW analyses have been saved in the file [" + setlstfile + "]."<<endl;
        fclose(smr);
        fclose(setlst);
        if (get_snp_effects_flg){
            fclose(snpfile);
        }      
    }

    void smr_e2e_prbmatch(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp)
    {
        //here eqtlFileName is the outcome and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not supposed to be used together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not supposed to be used together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not supposed to be used together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not supposed to be used together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        eqtlInfo etrait;
        eqtlInfo esdata;
        bInfo bdata;
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(eqtlFileName==NULL) throw("Error: please input QTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&etrait, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait, snplst2exclde);
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&etrait, problstName);
        if(problst2exclde != NULL) exclude_prob(&etrait, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&etrait, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&etrait, oprobe);
        if(oproblst2exclde != NULL) exclude_prob(&etrait, oproblst2exclde);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&etrait, oprobe2rm);
        
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        //if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName); //no need here, already extracted in etrait
        //if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
   
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        //if(snplstName != NULL) extract_snp(&bdata, snplstName);
        //if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &etrait, &esdata);
        // if no snp left after check
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &etrait, &esdata);
        }
        if(forcefrqck)
        {
            double prop=freq_check(&bdata, &etrait, &esdata,afthresh,percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
    
        //the etrait is not updated, so from now on _esi_include should be used always.
        cout<<"Reading exposure QTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(snpproblstName) extract_targets(&esdata, snpproblstName,  prb_snp);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("ERROR: no data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }

        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }   

        string outstr="";
        string smrfile = string(outFileName)+".msmr";
        FILE* smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }
        outstr="Expo_ID\tExpo_Chr\tExpo_Name\tExpo_bp\tOutcome_ID\tOutcome_Chr\tOutcome_Name\tOutcome_bp\ttopSNP\tA1\tA2\ttopSNP_chr\ttopSNP_bp\tb_top\tse_top\tp_top\tb_ivw\tse_ivw\tp_ivw\tn_iv\th2cis\n";
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }

        FILE* snpfile;
        string snpfilename = string(outFileName) + ".snps";
        string snpstr;
        if (get_snp_effects_flg){
            snpfile = fopen(snpfilename.c_str(), "w");
            if (!(snpfile)) {
                printf("ERROR: open error %s\n", snpfilename.c_str());
                exit(1);
            }
            snpstr="Exposure\tSNP\tA1\tA2\tbeta.exp\tse.exp\tbeta.out\tse.out\n";
            if(fputs_checked(snpstr.c_str(),snpfile))
            {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }    
       
        vector<int> includebk=esdata._include;

        logstr= "Number of exposure QTL SNPs:" + atos(esdata._esi_include.size())+ "\n";
        cout<<logstr<<endl;
        //double disp=0.0;

        cis_itvl=cis_itvl*1000;
        long write_count=0;
        SMRWK smrwk;

        unsigned int probNum = etrait._probNum;

        float progr0=0.0 , progr1;
        progress_print(progr0);

        for( int ii=0;ii<probNum;ii++){

            progr1=1.0*ii/probNum;
                if(progr1-progr0-0.01>1e-6 || ii+1==probNum)
                {
                    if(ii+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
               }

            //progress(ii,disp,(int)etrait._probNum);
            //
            string traitname=etrait._epi_prbID[ii];
            string traitgene=etrait._epi_gene[ii];
            cout<<"\nStarting with probe [ "+traitname+" ]..."<<endl;

            // find the matching exposure probe to outcome probe
            int i=-9;
            string exponame;
            string expogene;

            for (int j = 0; j<includebk.size(); j++){
                exponame=esdata._epi_prbID[includebk[j]];
                expogene=esdata._epi_gene[includebk[j]];

                if ((traitgene==expogene) || (traitgene==exponame) || (traitname==expogene) || (traitname==exponame)){
                    i=includebk[j];
                    break;
                }

            }      

            if (i==-9){
                printf("No matching exposure probe found for output probe [ %s ].\n", traitname.c_str());
                continue;
            }

            int expobp=esdata._epi_bp[i];
            int expochr=esdata._epi_chr[i];
            exponame=esdata._epi_prbID[i];
            expogene=esdata._epi_gene[i];

            printf("Matching exposure probe: [ %s ].\n", exponame.c_str());

            gwasData gdata;
            gdata.allele_1.resize(etrait._esi_include.size());
            gdata.allele_2.resize(etrait._esi_include.size());
            gdata.byz.resize(etrait._esi_include.size());
            gdata.seyz.resize(etrait._esi_include.size());
            gdata.freq.resize(etrait._esi_include.size());
            gdata.pvalue.resize(etrait._esi_include.size());
            gdata.splSize.resize(etrait._esi_include.size());
            gdata.snpName.resize(etrait._esi_include.size());
            
            for(int j=0;j<etrait._esi_include.size();j++){
                gdata.seyz[j]=-9;
                gdata.pvalue[j]=-9;
            }
            gdata._include.clear();
            
            int traitchr=etrait._epi_chr[ii];
            int traitbp=etrait._epi_bp[ii];
       
            int count=0;
            if(etrait._rowid.empty())
            {
                for (int j = 0; j<etrait._esi_include.size(); j++)
                {
                    if (fabs(etrait._sexz[ii][etrait._esi_include[j]] + 9) > 1e-6)
                    {
                        gdata.byz[j]=etrait._bxz[ii][etrait._esi_include[j]];
                        gdata.seyz[j]=etrait._sexz[ii][etrait._esi_include[j]];
                        double z=etrait._bxz[ii][etrait._esi_include[j]]/etrait._sexz[ii][etrait._esi_include[j]];
                        double p=pchisq(z*z,1);
                        gdata.pvalue[j]=p;
                        gdata.snpName[j]=etrait._esi_rs[etrait._esi_include[j]];
                        gdata.allele_1[j]=etrait._esi_allele1[etrait._esi_include[j]];
                        gdata.allele_2[j]=etrait._esi_allele2[etrait._esi_include[j]];
                        gdata._include.push_back(etrait._esi_include[j]); // row id selected
                        count++;
                    }
                }
            }
            else
            {
                uint64_t beta_start=etrait._cols[ii<<1];
                uint64_t se_start=etrait._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=etrait._rowid[beta_start+j];
                    int idx=(int)(find(etrait._esi_include.begin(), etrait._esi_include.end(), ge_rowid)-etrait._esi_include.begin());
                    if(idx<etrait._esi_include.size())
                    {
                        gdata.byz[idx]=etrait._val[beta_start+j];
                        gdata.seyz[idx]=etrait._val[se_start+j];
                        double z=etrait._val[beta_start+j]/etrait._val[se_start+j];
                        double p=pchisq(z*z,1);
                        gdata.pvalue[idx]=p;
                        gdata.snpName[idx]=etrait._esi_rs[ge_rowid];
                        gdata.allele_1[idx]=etrait._esi_allele1[ge_rowid];
                        gdata.allele_2[idx]=etrait._esi_allele2[ge_rowid];
                        gdata._include.push_back(idx);
                        count++;
                    }
                    
                }
            }

            gdata.snpNum=gdata.snpName.size();

            // fill eqtl data
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=i;
            //cout<<"\nInitiating the workspace of probe [ "+traitname+" ]..."<<endl;
            //printf("\nInitiating the workspace of probe %s for IVW multi-SNP SMR analysis....\n",exponame.c_str());
            char* refSNP=NULL;
            bool Flag = false;
            long maxid =fill_smr_wk(&bdata, &gdata, &esdata, &smrwk, refSNP, cis_itvl, Flag);      

            if (smrwk.bxz.size() == 0) {
                printf("WARNING: no SNP fetched for probe [ %s ].\n", exponame.c_str());
                continue;
            }
            printf("%ld SNPs are included from the cis-region of the probe [ %s ].\n",smrwk.bxz.size(),exponame.c_str());
            //now if you sepcify reference SNP, maxid point to this SNP, otherwise maxid is -9
            // step2: get top-SNP
            printf("Checking the top-SNP in the region....\n");
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            VectorXd zsxz;
            zsxz=ei_bxz.array()/ei_sexz.array();
            maxid=max_abs_id(zsxz); // now maxid point to the sig eQTL SNP in the new datastruct(not the raw).
            double pxz_val = pnorm(fabs(zsxz[maxid]));
            //double computing, consistency should be checked
            for(int tid=0;tid<zsxz.size();tid++) {
                if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                {
                    //printf("ERROR: zxz not consistent %f vs %f. please report this.\n",zsxz(tid),smrwk.zxz[tid]);
                    exit(EXIT_FAILURE);
                }
            }

            string topsnpname=smrwk.rs[maxid];
            printf("The top SNP of probe [ %s ] is %s with p-value %e.\n", exponame.c_str(), topsnpname.c_str(),pxz_val);
            if(pxz_val>p_smr){
                printf("WARNING: no SNP passed the p-value threshold %e for MR-IVW analysis for probe [ %s ].\n", p_smr, exponame.c_str());
                continue;
            }
            //cout<<maxid<<":"<<topsnpname<<":"<<esdata._esi_rs[smrwk.curId[maxid]]<<":"<<bdata._snp_name[bdata._include[smrwk.curId[maxid]]]<<":"<<gdata.snpName[smrwk.curId[maxid]]<<endl;
            // step3: extract SNPs around the --set-wind around sig (or ref) SNP
            vector<uint32_t> slctId;
            vector<int> slct_bpsnp,slct_snpchr;
            vector<double> slct_bxz, slct_sexz, slct_pxz, slct_byz, slct_seyz, slct_zsxz,slct_zxz, slct_pyz,slct_freq; //slct_zsxz,slct_zxz would be removed one of them
            vector<string> slct_snpName,slct_a1, slct_a2;
            long slct_maxid=-9;
            slctId.swap(smrwk.curId);
            slct_bxz.swap(smrwk.bxz);
            slct_sexz.swap(smrwk.sexz);
            slct_byz.swap(smrwk.byz);
            slct_seyz.swap(smrwk.seyz);
            slct_snpName.swap(smrwk.rs);
            slct_maxid=maxid;
            for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
            slct_zxz.swap(smrwk.zxz);
            slct_pyz.swap(smrwk.pyz);
            slct_bpsnp.swap(smrwk.bpsnp);
            slct_snpchr.swap(smrwk.snpchrom);
            slct_a1.swap(smrwk.allele1);
            slct_a2.swap(smrwk.allele2);
            slct_freq.swap(smrwk.freq);
            
            int out_raw_id = slctId[slct_maxid];
            double bxz_max = slct_bxz[slct_maxid];
            double sexz_max = slct_sexz[slct_maxid];
            double byz_max = slct_byz[slct_maxid];
            double seyz_max = slct_seyz[slct_maxid];
            double bxy_max = byz_max / bxz_max;
            //double sexy_max = sqrt(pow(seyz_max, 2)/pow(bxz_max,2) + pow(byz_max*sexz_max,2)/pow(bxz_max,4));
            double sexy_max = sqrt(pow(seyz_max, 2)/pow(bxz_max,2));
            double zxy_max = bxy_max / sexy_max;
            double pxy_max = 2*pnorm(fabs(zxy_max));
            string a1_max = slct_a1[slct_maxid];
            string a2_max = slct_a2[slct_maxid];
            
            vector<string> snp4msmr;
            double beta_ivw, se_ivw, p_ivw, h2cis;
            beta_ivw = -9;
            se_ivw = -9;
            h2cis = -9;
            bool p_expo_provided = false;
            int snp_count=smr_ivw_test(&bdata, slctId, slct_snpName, slct_a1, slct_a2, slct_bxz,slct_sexz, slct_pxz, slct_byz,slct_seyz, beta_ivw, se_ivw, p_ivw, h2cis, p_smr, ld_top_multi, ldmatrix, N_ref, trev, snp4msmr, p_expo_provided, get_snp_effects_flg, min_snp, snpfile, snpfilename, exponame);
            if(snp_count==-9) continue;
            // output snp set list
            
            string setstr=exponame+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<snp4msmr.size();j++)
            {
                setstr=snp4msmr[j]+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            setstr="end\n";
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            // end of output
            
            outstr = exponame + '\t' + atos(expochr) + '\t' + expogene + '\t' + atos(expobp) + '\t' + traitname + '\t' + atos(traitchr) + '\t' + traitgene + '\t' + atos(traitbp) + '\t';
            outstr += topsnpname + '\t' + a1_max + '\t' + a2_max + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t';
            outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t';
            outstr += atos(beta_ivw) + '\t' + atos(se_ivw) + '\t' + dtos(p_ivw) + '\t' + atos(snp_count) + '\t' + atos(h2cis) + '\n';
            
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
            write_count++;
            
        }

        cout<<"\nMR-IVW analysis on QTL exposure and outcome QTL data (probe match) completed.\nMR analysis results of "<<write_count<<" probes have been saved in the file [" + smrfile + "]."<<endl;
       

        fclose(smr);
        fclose(setlst);
        if (get_snp_effects_flg){
            fclose(snpfile);
        }
        
    }

    void smr_e2e_cis(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool multi_corr_snp_flg, bool get_snp_effects_flg, int min_snp)
    {

        //here eqtlFileName is the outcome and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);

        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not supposed to be used together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not supposed to be used together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not supposed to be used together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not supposed to be used together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        eqtlInfo mdata;
        eqtlInfo esdata;
        bInfo bdata;
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        read_esifile(&mdata, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&mdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&mdata, snplst2exclde);
        read_epifile(&mdata, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&mdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&mdata, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&mdata, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&mdata, oprobe);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&mdata, oprobe2rm);
        
        read_besdfile(&mdata, string(eqtlFileName)+".besd");
        if(mdata._rowid.empty() && mdata._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        //if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName); //no need here, already extracted in etrait
        //if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
   
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        //if(snplstName != NULL) extract_snp(&bdata, snplstName);
        //if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &mdata, &esdata);
        // if no snp left after check
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &mdata, &esdata);
        }
        if(forcefrqck)
        {
            double prop=freq_check(&bdata, &mdata, &esdata,afthresh,percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
    
        //the mdata is not updated, so from now on _esi_include should be used always.
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(snpproblstName) extract_targets(&esdata, snpproblstName,  prb_snp);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("ERROR: no data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }

        string outstr="";
        string smrfile = string(outFileName)+".msmr";
        FILE* smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }

        outstr="Expo_ID\tExpo_Name\tChr\tExpo_bp\tOutcome_ID\tOutcome_Name\tOutcome_bp\tb_ivw\tse_ivw\tp_ivw\tn_iv\n";
        
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }
       
        vector<int> includebk=esdata._include;
        vector<int> includemed=mdata._include;

        logstr= atos(esdata._esi_include.size())+ " SNPs are included in the analysis.\n";
        cout<<logstr<<endl;
        //double disp=0.0;

        cis_itvl=cis_itvl*1000;
        long write_count=0;
        SMRWKMULT smrwk;

        unsigned int probNum = esdata._include.size();
        int medNum = mdata._include.size();

        float progr0=0.0 , progr1;
        progress_print(progr0);

        for( int ii_r=0;ii_r<probNum;ii_r++){

            progr1=1.0*ii_r/probNum;
                if(progr1-progr0-0.01>1e-6 || ii_r+1==probNum)
                {
                    if(ii_r+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }
            
            int ii = includebk[ii_r];

            string exponame=esdata._epi_prbID[ii];
            string expogene=esdata._epi_gene[ii];
        
            int expochr=esdata._epi_chr[ii];
            int expobp=esdata._epi_bp[ii];

            cout<<"\nStarting with exposure [ " + exponame + " ]..."<<endl;

            // find mediator probes within a window of cis_itvl

            vector<string> med_probe_name;
            vector<string> med_probe_ID;
            vector<int> med_probe_bp;
            vector<uint32_t> med_probe_idx;
            int medchr, medbp;

            for (int j = 0; j<medNum; j++){
                medchr=mdata._epi_chr[includemed[j]];
                medbp=mdata._epi_bp[includemed[j]];

                if ((medchr == expochr) && (abs(medbp - expobp) <= cis_itvl)){
                    med_probe_ID.push_back(mdata._epi_prbID[includemed[j]]);
                    med_probe_name.push_back(mdata._epi_gene[includemed[j]]);
                    med_probe_bp.push_back(medbp);
                    med_probe_idx.push_back(includemed[j]);
                }
            }
            
            int medNum_cis = med_probe_idx.size();

            if (medNum_cis < 1){
                cout<<"No mediator probes found for [ " + exponame + " ]..."<<endl;
                continue;
            } else {
                cout<< atos(medNum_cis) + " mediator probes found for [ " + exponame + " ]..."<<endl;
            }
 
            MatrixXd Xm_b_raw(mdata._esi_include.size(), medNum_cis), Xm_se_raw(mdata._esi_include.size(), medNum_cis);

            for(int j=0;j<mdata._esi_include.size();j++){
                for (int k=0;k<medNum_cis;k++){
                    Xm_b_raw(j,k) = -9;
                    Xm_se_raw(j,k) = -9;
                }
            }

            int jj;

            for (int k=0;k<medNum_cis;k++){

                jj = med_probe_idx[k];

                if(mdata._rowid.empty()){

                    for (int j = 0; j<mdata._esi_include.size(); j++){
                        if (fabs(mdata._sexz[jj][j] + 9) > 1e-6){
                            Xm_b_raw(j,k) = mdata._bxz[jj][j];
                            Xm_se_raw(j,k) = mdata._sexz[jj][j];
                        }             
                    }
                } else {

                    uint64_t beta_start=mdata._cols[jj<<1];
                    uint64_t se_start=mdata._cols[1+(jj<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=mdata._rowid[beta_start+j];
                        int idx=(int)(find(mdata._esi_include.begin(), mdata._esi_include.end(), ge_rowid)-mdata._esi_include.begin());
                        if(idx<mdata._esi_include.size())
                        {
                            Xm_b_raw(idx,k) = mdata._val[beta_start+j];
                            Xm_se_raw(idx,k) = mdata._val[se_start+j];
                        }                    
                    }
                }
            }

            // fill data, SNPs within 2*cis_itvl are included
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=ii;
            smrwk.medNum_cis = medNum_cis;
            bool in_med_cis;

            if(esdata._rowid.empty()){
                for (int j = 0; j<esdata._esi_include.size(); j++){

                    int snpbp=esdata._esi_bp[j];
                    int snpchr=esdata._esi_chr[j];

                    in_med_cis = true;
                    for (int k=0;k<medNum_cis;k++){
                        if (Xm_se_raw(j,k) < -8){
                            in_med_cis = false;
                        }
                    }

                    if(snpchr==esdata._epi_chr[ii] && abs(esdata._epi_bp[ii]-snpbp)<=2*cis_itvl && in_med_cis)
                    {
                        smrwk.bxz.push_back(esdata._bxz[ii][j]);
                        smrwk.sexz.push_back(esdata._sexz[ii][j]);
                        smrwk.zxz.push_back(esdata._bxz[ii][j]/esdata._sexz[ii][j]);
                        smrwk.curId.push_back(j); //save snp id of the raw datastruct
                        smrwk.rs.push_back(esdata._esi_rs[j]);
                        smrwk.snpchrom.push_back(esdata._esi_chr[j]);
                        smrwk.allele1.push_back(esdata._esi_allele1[j]);
                        smrwk.allele2.push_back(esdata._esi_allele2[j]);
                        smrwk.bpsnp.push_back(esdata._esi_bp[j]);
                        smrwk.freq.push_back(bdata._mu[bdata._include[j]] / 2);
                    }
                }
            } 
            else {

                uint64_t beta_start=esdata._cols[ii<<1];
                uint64_t se_start=esdata._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                
                
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=esdata._rowid[beta_start+j];
                    int snpbp=esdata._esi_bp[ge_rowid];
                    int snpchr=esdata._esi_chr[ge_rowid];

                    in_med_cis = true;
                    for (int k=0;k<medNum_cis;k++){
                        if (Xm_se_raw(ge_rowid,k) < -8){
                            in_med_cis = false;
                        }
                    }

                    if(snpchr==esdata._epi_chr[ii] && abs(esdata._epi_bp[ii]-snpbp)<=2*cis_itvl && in_med_cis)
                    {
                        smrwk.bxz.push_back(esdata._val[beta_start+j]);
                        smrwk.sexz.push_back(esdata._val[se_start+j]);
                        smrwk.zxz.push_back(esdata._val[beta_start+j]/esdata._val[se_start+j]);
                        smrwk.curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk.rs.push_back(esdata._esi_rs[ge_rowid]);
                        smrwk.snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                        smrwk.allele1.push_back(esdata._esi_allele1[ge_rowid]);
                        smrwk.allele2.push_back(esdata._esi_allele2[ge_rowid]);
                        smrwk.bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                        smrwk.freq.push_back(bdata._mu[bdata._include[ge_rowid]] / 2);
                    }
                }
            }

            int num_cis_snps = smrwk.bxz.size();

            if (num_cis_snps == 0) {
                cout<<"WARNING: no SNP fetched for probe [ " + exponame + " ]."<<endl;
                //printf("WARNING: no SNP fetched for probe [ %s ].\n", exponame.c_str());
                continue;
            }

            smrwk.Xm_b.resize(num_cis_snps, medNum_cis);
            smrwk.Xm_se.resize(num_cis_snps, medNum_cis);
        
            for (int k=0;k<medNum_cis;k++){
                for (int j=0;j<num_cis_snps;j++){
                    smrwk.Xm_b(j, k) =  Xm_b_raw(smrwk.curId[j], k);
                    smrwk.Xm_se(j, k) =  Xm_se_raw(smrwk.curId[j], k);
                }
            }

            cout<< atos(num_cis_snps) + " SNPs are included from the cis-region of the probe [ " + exponame + " ]."<<endl;
            //printf("%ld SNPs are included from the cis-region of the probe [ %s ].\n",smrwk.bxz.size(),exponame.c_str());
            //now if you sepcify reference SNP, maxid point to this SNP, otherwise maxid is -9
            // step2: get top-SNP
            printf("Checking the top-SNP in the region....\n");
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            VectorXd zsxz;
            zsxz=ei_bxz.array()/ei_sexz.array();
            
            int maxidx;
            maxidx=max_abs_id(zsxz); 
            double pxz_val = 2*pnorm(fabs(zsxz[maxidx]));

            //double computing, consistency should be checked
            for(int tid=0;tid<zsxz.size();tid++) {
                if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                {
                    //printf("ERROR: zxz not consistent %f vs %f. please report this.\n",zsxz(tid),smrwk.zxz[tid]);
                    exit(EXIT_FAILURE);
                }
            }

            string topsnpnamex=smrwk.rs[maxidx];
            cout<< "The top SNP of exposure probe [ " + exponame + " ] is " + topsnpnamex + " with p-value " + atos(pxz_val) <<endl;
           
           if((pxz_val>p_smr)){
                logstr="WARNING: no SNP passed the exposure p-value threshold " + atos(p_smr) + ".\n";
                cout<<logstr<<endl;
                continue;
            }     

            vector<string> snp4msmr;

            logstr="Conducting MR analysis for " + atos(medNum_cis) + " outcomes ...";
            cout<<logstr<<endl;


            vector<uint32_t> Id4multi;
            vector<double> zxz4multi;

            double z_x= fabs(qnorm5(p_smr/2));
            double ztmpm;

            for(int j=0;j<smrwk.curId.size();j++)
            {
                double ztmpx=fabs(smrwk.bxz[j]/smrwk.sexz[j]);
            
                if (ztmpx>=z_x){
                     
                    Id4multi.push_back(smrwk.curId[j]);
                    zxz4multi.push_back(ztmpx);
                }
            }

            int snp_count4multi= (int)Id4multi.size();

            if(snp_count4multi < 3){
               
                logstr="No MR analyses performed because less than 3 SNPs passed the p-value threshold.";
                cout<<logstr<<endl;
                continue;
            }

            MatrixXd _X4multi;
            make_XMat(&bdata, Id4multi, _X4multi);
            vector<int> sub_indx4multi;

            double sbat_ld_cutoff=sqrt(ld_top_multi);
            
            sbat_calcu_lambda(_X4multi, snp_count4multi, sbat_ld_cutoff, sub_indx4multi, zxz4multi);
            
            if(snp_count4multi< 3){
                logstr="No MR analyses performed because less than 3 SNPs passed the LD pruning threshold.";
                cout<<logstr<<endl;
                continue;
            }

            logstr=atos(snp_count4multi) + " SNPs passed LD-square threshold of " + atos(ld_top_multi) + " and " + atos(Id4multi.size()-sub_indx4multi.size()) + " SNPs are excluded in the analysis.";
            cout<<logstr<<endl;

            MatrixXd C;
            vector<uint32_t> Id4multi_pruned;
            MatrixXd _X;
            if (ldmatrix)
            {
                for (int j=0; j<snp_count4multi; j++) Id4multi_pruned.push_back(Id4multi[sub_indx4multi[j]]);
                make_XMat(&bdata, Id4multi_pruned, _X);
                cor_calc(C, _X);
            }

            // Calculate the effects to cis-mediators

            logstr="Calculating MR effects...";
            cout<<logstr<<endl;

            double b_ivw = -9, se_ivw = -9, p_ivw = -9, zrev; 
            int n_iv;
            VectorXd b_zx, b_zy, se_zx, se_zy;
            vector<int> snp_e_slct;
            MatrixXd C_sub;          

            for (int k=0; k<medNum_cis; k++){

                snp_e_slct.resize(0);
                for (int j=0; j<snp_count4multi;j++){
                    zrev = (fabs(smrwk.bxz[sub_indx4multi[j]]) - fabs(smrwk.Xm_b(sub_indx4multi[j], k)))/sqrt(pow(smrwk.sexz[sub_indx4multi[j]],2) + pow(smrwk.Xm_se(sub_indx4multi[j], k),2));
                    if (zrev > trev){
                        snp_e_slct.push_back(j);
                    }
                }
                n_iv = snp_e_slct.size();
                if (n_iv < 1) continue;

                b_zx.resize(n_iv), b_zy.resize(n_iv), se_zx.resize(n_iv), se_zy.resize(n_iv);

                for (int j=0; j<n_iv; j++){
                    b_zx[j] = smrwk.bxz[sub_indx4multi[snp_e_slct[j]]];
                    se_zx[j] = smrwk.sexz[sub_indx4multi[snp_e_slct[j]]];
                    b_zy[j] = smrwk.Xm_b(sub_indx4multi[snp_e_slct[j]], k);
                    se_zy[j] = smrwk.Xm_se(sub_indx4multi[snp_e_slct[j]], k);
                }

                if (ldmatrix) {   
                    C_sub.resize(n_iv, n_iv);
                    for (int i=0; i< n_iv; i++){
                        for (int j=i; j< n_iv; j++){
                            C_sub(i,j) = C(snp_e_slct[i], snp_e_slct[j]);
                            C_sub(j,i) = C_sub(i,j);
                        }
                    }
                    mr_ivw_LD(b_zx, b_zy, se_zx, se_zy, C_sub, b_ivw, se_ivw, p_ivw, N_ref);
                } else mr_ivw(b_zx, b_zy, se_zx, se_zy, b_ivw, se_ivw, p_ivw);

                outstr =  exponame + '\t' + expogene + '\t' + atos(expochr) + '\t' + atos(expobp) + '\t'; 
                outstr +=  med_probe_ID[k] + '\t' + med_probe_name[k] + '\t' + atos(med_probe_bp[k]) + '\t'; 
                outstr += atos(b_ivw)  + '\t' + atos(se_ivw) + '\t' + atos(p_ivw) + '\t' + atos(n_iv) + '\n';
                
                if(fputs_checked(outstr.c_str(),smr))
                {
                    printf("ERROR: in writing file %s .\n", smrfile.c_str());
                    exit(EXIT_FAILURE);
                }
                write_count++;
            }          
        }

        cout<<"\nMR-IVW analysis on QTL exposure and outcome QTL data in cis-region completed.\n Results of "<<write_count<<" probes have been saved in the file [" + smrfile + "]."<<endl;     

        fclose(smr);
        
    }
                                             
    int mr_mediation_analysis(bInfo* bdata, SMRWKM* smrwk, double &b_tot_ivw, double &se_tot_ivw, double &p_tot_ivw, double &b_direct, double &se_direct, double &p_direct, double &b_xm, double &se_xm, double &p_xm, double &b_my, double &se_my, double &p_my, int &n_snps_tot, int &n_snps_e, int &n_xm, double p_smr, double p_medsmr, double ld_top_multi, bool ldmatrix, int N_ref, double trev, vector<string> &snp4msmr, bool multi_corr_snp_flg, bool get_snp_effects_flg, int min_snp, FILE* snpfile, string snpfilename, string exponame, string medname)
    {
        string logstr;
        string snpstr;
        logstr="Conducting MR mediation analysis...";
        cout<<logstr<<endl;

        /* step4: Filter out the SNPs with p-smr/p_medsmr threshold and ld-pruning */

        // Id4multi includes SNPs for multivariable regression (exposure & mediator to outcome)
        // Id4mr includes SNPs for MR of exposure to outcome 
        // Id4xmmr includes SNPs for MR of exposure to mediator
        
        vector<uint32_t> Id4multiall, Id4multix, Id4multim, Id4mr, Id4xmmr;

        vector<double> bxz4multiall, bxz4mr, bxz4xmmr;
        vector<double> sexz4multiall, sexz4mr, sexz4xmmr;
        vector<double> bmz4multiall, bmz4mr, bmz4xmmr;
        vector<double> semz4multiall, semz4mr, semz4xmmr;
        vector<double> byz4multiall, byz4mr, seyz4multiall, seyz4mr;
        vector<double> zxz4multix, zmz4multim, zxz4mr, zxz4xmmr, zmz4mr, zyz4mr;
        vector<string> snpnameall, allele1all, allele2all, snpnamemr, allele1mr, allele2mr;
        bool expo_snp, med_snp;

        double z_x= fabs(qnorm5(p_smr/2));
        double z_m= fabs(qnorm5(p_medsmr/2));
        for(int j=0;j<smrwk->curId.size();j++)
        {
            double ztmpx=fabs(smrwk->bxz[j]/smrwk->sexz[j]);
            double ztmpm=fabs(smrwk->bmz[j]/smrwk->semz[j]);

            double zrevx = (fabs(smrwk->bxz[j]) - fabs(smrwk->byz[j]))/sqrt(pow(smrwk->sexz[j],2) + pow(smrwk->seyz[j],2)); 
            double zrevm = (fabs(smrwk->bmz[j]) - fabs(smrwk->byz[j]))/sqrt(pow(smrwk->semz[j],2) + pow(smrwk->seyz[j],2)); 
            double zrevxm = (fabs(smrwk->bxz[j]) - fabs(smrwk->bmz[j]))/sqrt(pow(smrwk->sexz[j],2) + pow(smrwk->semz[j],2));

            expo_snp = false; 
            med_snp = false;
            if ((zrevx > trev) && (ztmpx>=z_x)) expo_snp = true;
            if ((zrevm > trev) && (ztmpm>=z_m)) med_snp = true;

            if ((expo_snp && !med_snp) || (med_snp && !expo_snp) || (expo_snp && med_snp && (zrevxm > trev))){
                
                Id4multiall.push_back(smrwk->curId[j]);
                bxz4multiall.push_back(smrwk->bxz[j]);
                bmz4multiall.push_back(smrwk->bmz[j]);
                byz4multiall.push_back(smrwk->byz[j]);
                sexz4multiall.push_back(smrwk->sexz[j]);
                semz4multiall.push_back(smrwk->semz[j]);
                seyz4multiall.push_back(smrwk->seyz[j]);
                snpnameall.push_back(smrwk->rs[j]);
                allele1all.push_back(smrwk->allele1[j]);
                allele2all.push_back(smrwk->allele2[j]);
                
                if(ztmpx>=z_x){
                    Id4multix.push_back(smrwk->curId[j]);
                    zxz4multix.push_back(smrwk->bxz[j]/smrwk->sexz[j]);
                }

                if(ztmpm>=z_m){
                    Id4multim.push_back(smrwk->curId[j]);
                    zmz4multim.push_back(smrwk->bmz[j]/smrwk->semz[j]);
                }
            } 

            if (expo_snp){
                Id4mr.push_back(smrwk->curId[j]);
                bxz4mr.push_back(smrwk->bxz[j]);
                sexz4mr.push_back(smrwk->sexz[j]);
                zxz4mr.push_back(smrwk->bxz[j]/smrwk->sexz[j]);
                bmz4mr.push_back(smrwk->bmz[j]);
                semz4mr.push_back(smrwk->semz[j]);
                zmz4mr.push_back(smrwk->bmz[j]/smrwk->semz[j]);
                byz4mr.push_back(smrwk->byz[j]);
                seyz4mr.push_back(smrwk->seyz[j]);
                zyz4mr.push_back(smrwk->byz[j]/smrwk->seyz[j]);
                snpnamemr.push_back(smrwk->rs[j]);
                allele1mr.push_back(smrwk->allele1[j]);
                allele2mr.push_back(smrwk->allele2[j]);
            }

            if ((ztmpx>=z_x) && (zrevxm > trev)){
                Id4xmmr.push_back(smrwk->curId[j]);
                bxz4xmmr.push_back(smrwk->bxz[j]);
                bmz4xmmr.push_back(smrwk->bmz[j]);
                zxz4xmmr.push_back(smrwk->bxz[j]/smrwk->sexz[j]);
                sexz4xmmr.push_back(smrwk->sexz[j]);
                semz4xmmr.push_back(smrwk->semz[j]);
            }
        }
        
        int snp_count4multix = (int)Id4multix.size();
        int snp_count4multim = (int)Id4multim.size();
        int snp_count4multi= snp_count4multix + snp_count4multim;
        n_snps_tot=(int)Id4mr.size();
        n_xm=(int)Id4xmmr.size();

        logstr=atos(Id4multiall.size()) + " SNPs with stronger exposure and mediator than outcome effects passed the exposure p-value " + atos(p_smr) + " or mediator p-value " + atos(p_medsmr) + " (MVMR).";
        cout<<logstr<<endl;
        logstr=atos(Id4mr.size()) + " SNPs with stronger exposure than outcome effects passed the exposure p-value threshold " + atos(p_smr) + " (total MR).";
        cout<<logstr<<endl;
        if((snp_count4multi<3) || (n_snps_tot<1)){
            cout<< "Not enough instrumental variables to calculate MR effects." <<endl;
            return -9; // we need at least 3 SNPs to calculate a direct effect
        }
        /* step5: Pruning */
    
        // SNPs for multivariable regression, perform pruning once on exposure-associated SNPs and then with mediator-associated SNPs

        vector<int> sub_indx4multi, sub_indx4multix;
        MatrixXd _X4multi, _X4multix;
        make_XMat(bdata, Id4multix, _X4multix); //_X: one row one individual, one column one SNP 

        double sbat_ld_cutoff=sqrt(ld_top_multi);
        sbat_calcu_lambda(_X4multix, snp_count4multix,  sbat_ld_cutoff, sub_indx4multix, zxz4multix);

        n_snps_e = snp_count4multix;
        if(n_snps_e<3){
            cout << "Less than 3 exposure-associated SNPs. No mediation analysis conducted." << endl;
            return -9;           
        } 

        vector<uint32_t> Id4multi, Id4multix_pruned;
        for (int j = 0; j < sub_indx4multix.size(); j++) Id4multix_pruned.push_back(Id4multix[sub_indx4multix[j]]);

        // add M-SNPs and give highest pruning priority to E-SNPs

        vector<double> bxz4multi, bmz4multi, byz4multi, sexz4multi, semz4multi, seyz4multi, pruning_weights;
        vector<string> snpnamemulti, allele1multi, allele2multi;

        for(int j=0; j<Id4multiall.size();j++)
        {
            med_snp = false;
            expo_snp = find(Id4multix_pruned.begin(), Id4multix_pruned.end(), Id4multiall[j]) != Id4multix_pruned.end();
            
            int midx=(int)(find(Id4multim.begin(), Id4multim.end(), Id4multiall[j])-Id4multim.begin());
            if(midx<Id4multim.size()) med_snp = true;

            if (expo_snp || med_snp){               
                Id4multi.push_back(Id4multiall[j]);
                bxz4multi.push_back(bxz4multiall[j]);
                bmz4multi.push_back(bmz4multiall[j]);
                byz4multi.push_back(byz4multiall[j]);
                sexz4multi.push_back(sexz4multiall[j]);
                semz4multi.push_back(semz4multiall[j]);
                seyz4multi.push_back(seyz4multiall[j]);
                snpnamemulti.push_back(snpnameall[j]);
                allele1multi.push_back(allele1all[j]);
                allele2multi.push_back(allele2all[j]);
                if (expo_snp){
                    pruning_weights.push_back(100);
                }
                else {
                    pruning_weights.push_back(zmz4multim[midx]);
                }
            }
        }

        snp_count4multi = Id4multi.size();
              
        make_XMat(bdata, Id4multi, _X4multi);
        sbat_calcu_lambda(_X4multi, snp_count4multi, sbat_ld_cutoff, sub_indx4multi, pruning_weights);

        vector<uint32_t> Id4multi_pruned;
        for (int j = 0; j < sub_indx4multi.size(); j++) Id4multi_pruned.push_back(Id4multi[sub_indx4multi[j]]);

        logstr=atos(sub_indx4multi.size()) + " SNPs passed LD-square threshold of " + atos(ld_top_multi) + " and " + atos(Id4multiall.size()-sub_indx4multi.size()) + " SNPs are excluded in the multivariable effect calculation.";
        cout<<logstr<<endl;
        logstr=atos(sub_indx4multi.size()) + " SNPs are included in the multivariable effect calculation.";
        cout<<logstr<<endl;

        // SNPs for MR total effect

        vector<int> sub_indx4mr;
        MatrixXd _X4mr;
        make_XMat(bdata, Id4mr, _X4mr); //_X: one row one individual, one column one SNP

        sbat_calcu_lambda(_X4mr, n_snps_tot,  sbat_ld_cutoff, sub_indx4mr, zxz4mr); //the index of slectId, snp_count can change here, pruning is based on eQTL p-values (zyz4smr is not used)
        logstr=atos(sub_indx4mr.size()) + " SNPs passed LD-square threshold of " + atos(ld_top_multi) + " and " + atos(Id4mr.size()-sub_indx4mr.size()) + " SNPs are excluded in the MR total effect calculation.";
        cout<<logstr<<endl;
        vector<uint32_t> Id4mr_pruned;
        for (int j = 0; j < sub_indx4mr.size(); j++) Id4mr_pruned.push_back(Id4mr[sub_indx4mr[j]]);

        // SNPs for MR exposure to mediator
        vector<int> sub_indx4xmmr;
        vector<uint32_t> Id4xmmr_pruned;
        if (n_xm > 0)
        {
            MatrixXd _X4xmmr;
            make_XMat(bdata,Id4xmmr, _X4xmmr); //_X: one row one individual, one column one SNP
            sbat_calcu_lambda(_X4xmmr, n_xm, sbat_ld_cutoff, sub_indx4xmmr, zxz4xmmr); //the index of slectId, snp_count can change here, pruning is based on eQTL p-values (zyz4smr is not used)
            for (int j = 0; j < sub_indx4xmmr.size(); j++) Id4xmmr_pruned.push_back(Id4xmmr[sub_indx4xmmr[j]]);
        }

        /* step5: Calculate MR effects
        //Calculate the direct and mediator effect*/
        
        // Get C matrix
        MatrixXd C;
        if (ldmatrix)
        {
            _X4multi.resize(0,0);
            make_XMat(bdata, Id4multi_pruned, _X4multi);
            cor_calc(C, _X4multi);
        }
        else C.setIdentity(Id4multi_pruned.size(), Id4multi_pruned.size());
        
        // Define regression variables
        MatrixXd X(snp_count4multi, 2);
        VectorXd yz(snp_count4multi), SEs(3*snp_count4multi), vars, betas;

        for(int j=0; j<sub_indx4multi.size();j++)
        {
            yz(j) = byz4multi[sub_indx4multi[j]];
            X(j,0) = bxz4multi[sub_indx4multi[j]];
            X(j,1) = bmz4multi[sub_indx4multi[j]];
            SEs(j) = sexz4multi[sub_indx4multi[j]];
            SEs(snp_count4multi + j) = semz4multi[sub_indx4multi[j]];
            SEs(snp_count4multi*2 + j) = seyz4multi[sub_indx4multi[j]];

            if (get_snp_effects_flg && (snp_count4multi >= min_snp)){
                snpstr = exponame + '\t' + medname + '\t' + snpnamemulti[sub_indx4multi[j]] + '\t' + allele1multi[sub_indx4multi[j]] + '\t' + allele2multi[sub_indx4multi[j]] + '\t';
                snpstr += atos(bxz4multi[sub_indx4multi[j]]) + '\t' + atos(sexz4multi[sub_indx4multi[j]]) + '\t';
                snpstr +=  atos(bmz4multi[sub_indx4multi[j]]) + '\t' + atos(semz4multi[sub_indx4multi[j]]) + '\t';
                snpstr += atos(byz4multi[sub_indx4multi[j]]) + '\t' + atos(seyz4multi[sub_indx4multi[j]]) + '\n'; 

                if(fputs_checked(snpstr.c_str(),snpfile))
                {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
                }
            }           
        }

        mvmr_ivw_LD(X, yz, SEs, C, betas, vars, N_ref);

        b_direct = betas(0);
        b_my = betas(1);

        se_direct = sqrt(vars(0));
        se_my = sqrt(vars(1));

        p_direct = 2*pnorm(fabs(b_direct/se_direct));
        p_my = 2*pnorm(fabs(b_my/se_my));

        cout << "Direct MR effects: b = " << b_direct << "; p = " << p_direct << endl;

        //Calculate the total effect 
        VectorXd bzx(sub_indx4mr.size()), bzy(sub_indx4mr.size()), sezx(sub_indx4mr.size()), sezy(sub_indx4mr.size());

        for(int j=0; j<sub_indx4mr.size();j++){
            bzx[j] = bxz4mr[sub_indx4mr[j]];
            bzy[j] = byz4mr[sub_indx4mr[j]];
            sezx[j] = sexz4mr[sub_indx4mr[j]];
            sezy[j] = seyz4mr[sub_indx4mr[j]];
        }
        
        if (ldmatrix) {
            MatrixXd _X;
            make_XMat(bdata, Id4mr_pruned, _X);
            C.resize(0,0);
            cor_calc(C, _X);      
            mr_ivw_LD(bzx, bzy, sezx, sezy, C, b_tot_ivw, se_tot_ivw, p_tot_ivw, N_ref);
        } else mr_ivw(bzx, bzy, sezx, sezy, b_tot_ivw, se_tot_ivw, p_tot_ivw);

        cout << "Total MR effects: b = " << b_tot_ivw << "; p = " << p_tot_ivw << endl;
        
        // Calculate the exposure to mediator effect

        if (n_xm > 0) {

            bzx.resize(sub_indx4xmmr.size()), bzy.resize(sub_indx4xmmr.size()), sezx.resize(sub_indx4xmmr.size()), sezy.resize(sub_indx4xmmr.size());
           
            for(int j=0; j<sub_indx4xmmr.size();j++)
            {
                bzx[j] = bxz4xmmr[sub_indx4xmmr[j]];
                bzy[j] = bmz4xmmr[sub_indx4xmmr[j]];
                sezx[j] = sexz4xmmr[sub_indx4xmmr[j]];
                sezy[j] = semz4xmmr[sub_indx4xmmr[j]];
            }

            if (ldmatrix) {
                MatrixXd _X;
                make_XMat(bdata, Id4xmmr_pruned, _X);
                C.resize(0,0);
                cor_calc(C, _X);      
                mr_ivw_LD(bzx, bzy, sezx, sezy, C, b_xm, se_xm, p_xm, N_ref);
            } else mr_ivw(bzx, bzy, sezx, sezy, b_xm, se_xm, p_xm);
        }

        for (int j = 0; j < Id4multi_pruned.size(); j++) snp4msmr.push_back(bdata->_snp_name[bdata->_include[Id4multi_pruned[j]]]);
        return snp_count4multi;
    }
           
    int mr_mediation_analysis_cis(bInfo* bdata, SMRWKMULT* smrwk, double &b_tot_ivw, double &se_tot_ivw, double &p_tot_ivw, double &b_top_direct, double &se_top_direct, double &p_top_direct, double &b_direct, double &se_direct, double &p_direct, int &n_snps_tot, int &n_snps_top, int &n_snps_e, double p_smr, double p_medsmr, int &medNum_included, vector<int> &med_sub_indx, int &top_med_indx, double ld_top_multi, bool ldmatrix, int N_ref, double trev, vector<string> &snp4msmr, bool multi_corr_snp_flg, bool get_snp_effects_flg, bool get_snp_effects_top_flg, char* outFileName, int min_snp, string exponame, vector<string> &med_probe_names, vector<string> &med_probe_genes, vector<int> &med_probe_bp, MED_info* med_info, double med_R_thresh, bool uncorr_med, bool exclude_top_SNP, double p_shrinkage, double p_expo_med, string &mediation_status)
    {
        string logstr;
        string snpstr;
        logstr="Conducting MR-mediation analysis...";
        cout<<logstr<<endl;

        int medNum_cis = smrwk->medNum_cis;

        // Get the independent exposure SNPs

        vector<uint32_t> Id4x;
        vector<double> zxz4x;
        vector<int> slct_indx, slct_indx_tmp;
        
        double z_x= fabs(qnorm5(p_smr/2));

        for(int j=0;j<smrwk->curId.size();j++)
        {
            double ztmpx=fabs(smrwk->bxz[j]/smrwk->sexz[j]);
            double zrev = (fabs(smrwk->bxz[j]) - fabs(smrwk->byz[j]))/sqrt(pow(smrwk->sexz[j],2) + pow(smrwk->seyz[j],2)); 
         
            if ((ztmpx>=z_x) && (zrev>trev)){
 
                Id4x.push_back(smrwk->curId[j]);
                slct_indx_tmp.push_back(j);
                zxz4x.push_back(ztmpx);

            }
        }

        // Prune exposure SNPs

        int num_e_SNPs= (int)Id4x.size();

        logstr=atos(num_e_SNPs) + " SNPs with stronger exposure than outcome effects passed the exposure p-value " + atos(p_smr) + ".";
        cout<<logstr<<endl;
   
        double sbat_ld_cutoff=sqrt(ld_top_multi);
        MatrixXd _X;
        make_XMat(bdata, Id4x, _X);
        vector<int> sub_indx4x;     
        sbat_calcu_lambda(_X, num_e_SNPs,  sbat_ld_cutoff, sub_indx4x, zxz4x);
               
        if(num_e_SNPs < 3){
            cout << "Less than 3 independent exposure SNPs found. No MR analysis conducted." << endl;
            mediation_status = "Less than 3 independent E-SNPs.";
            return -9;
        } else {
            logstr=atos(num_e_SNPs) + " SNPs passed LD-square threshold of " + atos(ld_top_multi) + ".";
            cout<<logstr<<endl;
        }

        n_snps_tot = num_e_SNPs;
        VectorXd bzx_totmr(num_e_SNPs), sezx_totmr(num_e_SNPs), bzy_totmr(num_e_SNPs), sezy_totmr(num_e_SNPs);
        vector<uint32_t> Id4x_pruned;

        for (int j=0; j<num_e_SNPs; j++){
            slct_indx.push_back(slct_indx_tmp[sub_indx4x[j]]);
            Id4x_pruned.push_back(Id4x[sub_indx4x[j]]);
            bzx_totmr(j) = smrwk->bxz[slct_indx_tmp[sub_indx4x[j]]];
            sezx_totmr(j) = smrwk->sexz[slct_indx_tmp[sub_indx4x[j]]];
            bzy_totmr(j) = smrwk->byz[slct_indx_tmp[sub_indx4x[j]]];
            sezy_totmr(j) = smrwk->seyz[slct_indx_tmp[sub_indx4x[j]]];
        } 

        MatrixXd Ctot;
        if (ldmatrix)
        {
            make_XMat(bdata, Id4x_pruned, _X);
            cor_calc(Ctot, _X);
        }

        VectorXd bzx_totmr_extop(num_e_SNPs-1), sezx_totmr_extop(num_e_SNPs-1); 
        VectorXd bzy_totmr_extop(num_e_SNPs-1), sezy_totmr_extop(num_e_SNPs-1);
        MatrixXd Ctot_extop(num_e_SNPs-1, num_e_SNPs-1);
        if (exclude_top_SNP){
            int top_idx=0;
            double top_val=0;
            for (int j=0; j<num_e_SNPs; j++){
                if (fabs(bzx_totmr(j)/sezx_totmr(j)) > top_val){
                    top_idx = j;
                    top_val = fabs(bzx_totmr(j)/sezx_totmr(j));
                }
            }

            // Fill in effect sizes
            int j_new=0;
            for (int j=0; j<num_e_SNPs; j++){
                if (j==top_idx) continue;
                bzx_totmr_extop(j_new) = bzx_totmr(j);
                sezx_totmr_extop(j_new) = sezx_totmr(j);
                bzy_totmr_extop(j_new) = bzy_totmr(j);
                sezy_totmr_extop(j_new) = sezy_totmr(j);
                int ji_new =0;
                for (int ji=0; ji<num_e_SNPs; ji++){
                    if (ji==top_idx) continue;
                    Ctot_extop(j_new,ji_new) = Ctot(j,ji);
                    ji_new++;
                }
                j_new++;
            }
        }

        // Based on these SNPs calculate MR effects on mediator.

        double b_xm, se_xm, p_xm, zrevm;
        int n_xm;
        VectorXd expo_xz, med_mz, expo_sexz, med_semz;
        vector<double> cor_med_val, cor_med_val_slct; 
        vector<int> snp_e_slct, med_indx;
        MatrixXd Cxm;

        medNum_included = 0;

        for (int k=0; k<medNum_cis; k++){

            snp_e_slct.resize(0);
            for (int j=0; j<num_e_SNPs;j++){
                zrevm = (fabs(bzx_totmr[j]) - fabs(smrwk->Xm_b(slct_indx[j], k)))/sqrt(pow(sezx_totmr[j],2) + pow(smrwk->Xm_se(slct_indx[j], k),2));
                if (zrevm > trev){
                    snp_e_slct.push_back(j);
                }
            }
            n_xm = snp_e_slct.size();
            if (n_xm < 1){
                cor_med_val.push_back(0);
                continue;
            }

            expo_xz.resize(n_xm), med_mz.resize(n_xm), expo_sexz.resize(n_xm), med_semz.resize(n_xm);

            for (int j=0; j<n_xm; j++){
                expo_xz[j] = bzx_totmr[snp_e_slct[j]];
                expo_sexz[j] = sezx_totmr[snp_e_slct[j]];
                med_mz[j] = smrwk->Xm_b(slct_indx[snp_e_slct[j]], k);
                med_semz[j] = smrwk->Xm_se(slct_indx[snp_e_slct[j]], k);
            }

            if (ldmatrix) {   
                Cxm.resize(n_xm, n_xm);
                for (int i=0; i< n_xm; i++){
                    for (int j=i; j< n_xm; j++){
                        Cxm(i,j) = Ctot(snp_e_slct[i], snp_e_slct[j]);
                        Cxm(j,i) = Cxm(i,j);
                    }
                }
                mr_ivw_LD(expo_xz, med_mz, expo_sexz, med_semz, Cxm, b_xm, se_xm, p_xm, N_ref);
            } else mr_ivw(expo_xz, med_mz, expo_sexz, med_semz, b_xm, se_xm, p_xm);

            if (p_xm <= p_expo_med){
                medNum_included++;
                med_indx.push_back(k);
                cor_med_val.push_back(fabs(b_xm/se_xm));
                cor_med_val_slct.push_back(fabs(b_xm/se_xm));

                if (exclude_top_SNP){
                    if (n_xm > 1){
                        int top_idx=0;
                        double top_val=0;
                        for (int j=0; j<n_xm; j++){
                            if (fabs(expo_xz(j)/expo_sexz(j)) > top_val){
                                top_idx = j;
                                top_val = fabs(expo_xz(j)/expo_sexz(j));
                            }
                        }

                        VectorXd expo_xz_extop(n_xm-1), expo_sexz_extop(n_xm-1); 
                        VectorXd med_mz_extop(n_xm-1), med_semz_extop(n_xm-1);
                        MatrixXd Cxm_extop(n_xm-1, n_xm-1);

                        // Fill in effect sizes
                        int j_new=0;
                        for (int j=0; j<n_xm; j++){
                            if (j==top_idx) continue;
                            expo_xz_extop(j_new) = expo_xz(j);
                            expo_sexz_extop(j_new) = expo_sexz(j);
                            med_mz_extop(j_new) = med_mz(j);
                            med_semz_extop(j_new) = med_semz(j);
                            int ji_new =0;
                            for (int ji=0; ji<n_xm; ji++){
                                if (ji==top_idx) continue;
                                Cxm_extop(j_new,ji_new) = Cxm(j,ji);
                                ji_new++;
                            }
                            j_new++;
                        }
                        mr_ivw_LD(expo_xz_extop, med_mz_extop, expo_sexz_extop, med_semz_extop, Cxm_extop, b_xm, se_xm, p_xm, N_ref);
                        n_xm = n_xm-1;
                    }
                    else {
                        n_xm = 0;
                        b_xm = 0; 
                        se_xm = 0;
                        p_xm = 1;
                    }
                }

                med_info->Med_ID.push_back(med_probe_names[k]); 
                med_info->Gene_ID.push_back(med_probe_genes[k]);
                med_info->Med_bp.push_back(med_probe_bp[k]);  
                med_info->b_xm.push_back(b_xm);
                med_info->se_xm.push_back(se_xm);
                med_info->p_xm.push_back(p_xm);
                med_info->n_xm.push_back(n_xm);
            } else {
                cor_med_val.push_back(0);
            }
        }

        logstr="From " + atos(medNum_cis) + " mediators in cis " + atos(medNum_included) + " are causally-associated to the exposure probe.";
        cout<<logstr<<endl;

        if (medNum_included < 1){
            mediation_status = "No causally-associated mediator found.";
            return -9;
        }

        // Find the top-correlated mediator

        top_med_indx = max_element(cor_med_val.begin(), cor_med_val.end()) - cor_med_val.begin(); 

        // get the SNPs significantly associated to the included mediators

        slct_indx_tmp.resize(0);
        vector<int> slct_indx_tmp_top;
        bool plei_SNP, med_SNP;
        double z_m = fabs(qnorm5(p_medsmr/2));
        double ztmpx, ztmpm, zrevx;
        int snp_count4multi = 0;
        n_snps_top = 0;
        n_snps_e = 0;
        vector<double> zmz4multi, zmz4top;
        vector<int> slct_indx_e, slct_indx_etop;
        vector<uint32_t> Id4multi, Id4top;

        for(int j=0;j<smrwk->curId.size();j++){

            plei_SNP = false;
            med_SNP = false;

            ztmpx=fabs(smrwk->bxz[j]/smrwk->sexz[j]);
            zrevx = (fabs(smrwk->bxz[j]) - fabs(smrwk->byz[j]))/sqrt(pow(smrwk->sexz[j],2) + pow(smrwk->seyz[j],2)); 

            if ((ztmpx>=z_x) && (zrevx>trev)){
                if (find(slct_indx.begin(), slct_indx.end(), j) != slct_indx.end()){
                    zrevm = (fabs(smrwk->bxz[j]) - fabs(smrwk->Xm_b(j, top_med_indx)))/sqrt(pow(smrwk->sexz[j],2) + pow(smrwk->Xm_se(j,top_med_indx),2));
                    if (zrevm > trev){
                        n_snps_top++;
                        slct_indx_tmp_top.push_back(j);
                        zmz4top.push_back(100);
                        slct_indx_etop.push_back(n_snps_top-1);
                        Id4top.push_back(smrwk->curId[j]);
                    }
                    for (int k=0; k<medNum_included; k++){  
                        zrevm = (fabs(smrwk->bxz[j]) - fabs(smrwk->Xm_b(j, med_indx[k])))/sqrt(pow(smrwk->sexz[j],2) + pow(smrwk->Xm_se(j, med_indx[k]),2));        
                        if (zrevm < trev){
                            plei_SNP = true;
                            break;
                        }
                    }
                    if (!plei_SNP){
                        snp_count4multi++;
                        n_snps_e++;
                        slct_indx_tmp.push_back(j);
                        //zmz4multi.push_back(100);
                        slct_indx_e.push_back(snp_count4multi-1);
                        Id4multi.push_back(smrwk->curId[j]);
                    }
                }
            } else {
                ztmpm=fabs(smrwk->Xm_b(j, top_med_indx)/smrwk->Xm_se(j, top_med_indx));
                zrevm = (fabs(smrwk->Xm_b(j, top_med_indx)) - fabs(smrwk->byz[j]))/sqrt(pow(smrwk->seyz[j],2) + pow(smrwk->Xm_se(j, top_med_indx),2));
                if ((ztmpm >= z_m) && (zrevm > trev)){
                    n_snps_top++;
                    slct_indx_tmp_top.push_back(j);
                    zmz4top.push_back(ztmpm);
                    Id4top.push_back(smrwk->curId[j]);
                }
                for (int k=0; k<medNum_included; k++){
                    ztmpm=fabs(smrwk->Xm_b(j, med_indx[k])/smrwk->Xm_se(j, med_indx[k]));
                    if (ztmpm >= z_m) med_SNP = true;
                    zrevm = (fabs(smrwk->Xm_b(j, med_indx[k])) - fabs(smrwk->byz[j]))/sqrt(pow(smrwk->seyz[j],2) + pow(smrwk->Xm_se(j, med_indx[k]),2));
                    if (zrevm < trev) plei_SNP = true;
                }
                if (med_SNP && !plei_SNP){
                    snp_count4multi++;
                    slct_indx_tmp.push_back(j);
                    //zmz4multi.push_back(ztmpm);
                    Id4multi.push_back(smrwk->curId[j]);
                }
            }
        } 

        if(n_snps_e < 3){
            cout << "Less than 3 instrumental SNPs had larger exposure than mediator effects.\nNo MR analysis conducted." << endl;
            mediation_status = "Less than 3 E-SNPs.";
            return -9;
        }

        // Rank the SNPs

        MatrixXd Xm_rank(slct_indx_tmp.size(), med_indx.size());
        vector<double> tmp_z(slct_indx_tmp.size());
        for (int k=0; k<med_indx.size(); k++){
            for (int j=0; j<slct_indx_tmp.size(); j++){
                tmp_z[j] = fabs(smrwk->Xm_b(slct_indx_tmp[j], med_indx[k])/smrwk->Xm_se(slct_indx_tmp[j], med_indx[k]));
            }
            int r = 0;
            vector<int> sorted_z = sort_indexes(tmp_z);
            for (int i=0; i<sorted_z.size();i++) {
                Xm_rank(sorted_z[i], k) = r;
                r++;
            }
        }

        // SNP rank is sum of ranks along the mediator axis
        VectorXd snp_rank(slct_indx_tmp.size());
        snp_rank = Xm_rank.rowwise().sum();
        vector<double> snp_rank_double;
        for (int j=0; j<snp_rank.size();j++) snp_rank_double.push_back(snp_rank[j]);

        int max_rank = slct_indx_tmp.size()* med_indx.size(); // assigned to E-SNPs 
        for (int j= 0; j<slct_indx_e.size(); j++) snp_rank_double[slct_indx_e[j]] = max_rank;

        // Rank-based pruning for all cis-mediator SNPs

        make_XMat(bdata, Id4multi, _X);
        vector<int> sub_indx4multi_raw, sub_indx4multi_tmp, sub_indx4multi;     
        sbat_calcu_lambda(_X, snp_count4multi, sbat_ld_cutoff, sub_indx4multi_raw, snp_rank_double);
        for (int j=0; j<sub_indx4multi_raw.size(); j++) sub_indx4multi_tmp.push_back(slct_indx_tmp[sub_indx4multi_raw[j]]);

        // Find uncorrelated mediators if required

        if (uncorr_med){

            // Calculate correlation matrix among mediators
            MatrixXd X_raw(snp_count4multi, med_indx.size());

            for (int j=0; j<snp_count4multi; j++){
                for (int k=0; k<med_indx.size(); k++){
                    X_raw(j,k) = smrwk->Xm_b(sub_indx4multi_tmp[j], med_indx[k]);
                }
            }
                     
            MatrixXd C_med;
            cor_calc(C_med, X_raw);
            
            // prune the mediators; pruning priority is the correlation to the exposure
            int m = med_indx.size();  
            int qi = 0;    
            vector<int> rm_ID1, med_indx_raw;
            rm_cor_sbat(C_med, med_R_thresh, m, rm_ID1, cor_med_val_slct);
            //Create new index
            for (int i=0 ; i<m ; i++) {
                if (rm_ID1.size() == 0){
                    med_indx_raw.push_back(i);
                } else {
                    if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                    else med_indx_raw.push_back(i);
                }
            }
            medNum_included = med_indx_raw.size();
            cout << medNum_included << " mediators, uncorrelated among each other below R = " << med_R_thresh << " were found." << endl;

            for (int k=0; k < medNum_included; k++){
                med_sub_indx.push_back(med_indx[med_indx_raw[k]]);
            }

            // select SNPs associated to the uncorrelated mediators

            for(int j=0;j<sub_indx4multi_tmp.size();j++){

                ztmpx=fabs(smrwk->bxz[sub_indx4multi_tmp[j]]/smrwk->sexz[sub_indx4multi_tmp[j]]);
                if (ztmpx>=z_x){
                    sub_indx4multi.push_back(sub_indx4multi_tmp[j]);
                }
                else {
                    for (int k=0; k<medNum_included; k++){
                        ztmpm=fabs(smrwk->Xm_b(sub_indx4multi_tmp[j], med_sub_indx[k])/smrwk->Xm_se(sub_indx4multi_tmp[j], med_sub_indx[k]));
                        if (ztmpm >= z_m){
                            sub_indx4multi.push_back(sub_indx4multi_tmp[j]);
                            break;
                        }                        
                    }
                }
            }           
        } else {
            for (int k=0; k < med_indx.size(); k++){
                med_sub_indx.push_back(med_indx[k]);
            }
            for(int j=0;j<sub_indx4multi_tmp.size();j++){
                sub_indx4multi.push_back(sub_indx4multi_tmp[j]);
            }
        }

        snp_count4multi = sub_indx4multi.size();

        if(snp_count4multi < (medNum_included + 3)){
            cout << "Less than " << atos(medNum_included + 3) << " instrumental SNPs were found for both the exposure and mediators.\nNo MR analysis conducted." << endl;
            mediation_status = "Not enough instrumental variables.";
            return -9;
        } else {
            logstr=atos(snp_count4multi) + " SNPs passed LD-square threshold of " + atos(ld_top_multi) + " and are included in the cis-mediation analysis.";
            cout<<logstr<<endl;
        }

        MatrixXd X(snp_count4multi, medNum_included + 1), X_z(snp_count4multi, medNum_included + 1);
        VectorXd yz(snp_count4multi), SEs((medNum_included + 2)*snp_count4multi);
        vector<uint32_t> Id4cis_pruned;
        vector<int> snp_indx_e;

        for (int j=0; j<snp_count4multi; j++){
            Id4cis_pruned.push_back(smrwk->curId[sub_indx4multi[j]]);
            if (fabs(smrwk->bxz[sub_indx4multi[j]]/smrwk->sexz[sub_indx4multi[j]])> z_x) snp_indx_e.push_back(j);
            X(j,0) = smrwk->bxz[sub_indx4multi[j]];
            yz(j) = smrwk->byz[sub_indx4multi[j]];
            SEs(j) = smrwk->sexz[sub_indx4multi[j]];
            SEs((medNum_included+1)*snp_count4multi+j) = smrwk->seyz[sub_indx4multi[j]];
            X_z(j,0) = fabs(smrwk->bxz[sub_indx4multi[j]]/smrwk->sexz[sub_indx4multi[j]]);
            snp4msmr.push_back(smrwk->rs[sub_indx4multi[j]]);
            for (int k=0; k<medNum_included; k++){
                X(j,k+1) = smrwk->Xm_b(sub_indx4multi[j], med_sub_indx[k]);
                SEs((k+1)*snp_count4multi+j) = smrwk->Xm_se(sub_indx4multi[j], med_sub_indx[k]);
                X_z(j,k+1) = fabs(smrwk->Xm_b(sub_indx4multi[j], med_sub_indx[k])/smrwk->Xm_se(sub_indx4multi[j], med_sub_indx[k]));
            }
        }

        MatrixXd Ccis;
        if (ldmatrix)
        {
            make_XMat(bdata, Id4cis_pruned, _X);
            cor_calc(Ccis, _X);
        }
        else Ccis.setIdentity(Id4cis_pruned.size(), Id4cis_pruned.size());

        // find idx of top SNP if exclude_top_SNP activated
        MatrixXd X_extop(snp_count4multi-1, medNum_included + 1);
        VectorXd yz_extop(snp_count4multi-1), SEs_extop((medNum_included + 2)*(snp_count4multi-1));
        MatrixXd Ccis_extop(snp_count4multi-1, snp_count4multi-1);
        if (exclude_top_SNP){
            int top_idx=0;
            double top_val=0;
            for (int j=0; j<snp_count4multi; j++){
                if (X_z(j,0) > top_val){
                    top_idx = j;
                    top_val = X_z(j,0);
                }
            }

            // Fill in effect sizes
            int j_new=0;
            for (int j=0; j<snp_count4multi; j++){
                if (j==top_idx) continue;
                X_extop(j_new,0) = X(j,0);
                yz_extop(j_new) = yz(j);
                SEs_extop(j_new) = SEs(j);
                SEs_extop((medNum_included+1)*(snp_count4multi-1)+j_new) = SEs((medNum_included+1)*(snp_count4multi-1)+j);
                for (int k=0; k<medNum_included; k++){
                    X_extop(j_new,k+1) = X(j,k+1);
                    SEs_extop((k+1)*(snp_count4multi-1)+j_new) = SEs((k+1)*(snp_count4multi-1)+j);
                }
                int ji_new =0;
                for (int ji=0; ji<snp_count4multi; ji++){
                    if (ji==top_idx) continue;
                    Ccis_extop(j_new,ji_new) = Ccis(j,ji);
                    ji_new++;
                }
                j_new++;
            }
        }

        // Pruning for top-mediator SNPs
        make_XMat(bdata, Id4top, _X);
        vector<int> sub_indx4top;     
        sbat_calcu_lambda(_X, n_snps_top, sbat_ld_cutoff, sub_indx4top, zmz4top);

        MatrixXd Xtop(n_snps_top, 2);
        VectorXd yztop(n_snps_top), SEstop(3*n_snps_top);
        vector<uint32_t> Id4top_pruned;

        for (int j=0; j<n_snps_top; j++){
            Id4top_pruned.push_back(Id4top[sub_indx4top[j]]);
            Xtop(j,0) = smrwk->bxz[slct_indx_tmp_top[sub_indx4top[j]]];        
            Xtop(j,1) = smrwk->Xm_b(slct_indx_tmp_top[sub_indx4top[j]], top_med_indx);
            yztop(j) = smrwk->byz[slct_indx_tmp_top[sub_indx4top[j]]];
            SEstop(j) = smrwk->sexz[slct_indx_tmp_top[sub_indx4top[j]]];
            SEstop(n_snps_top+j) = smrwk->Xm_se(slct_indx_tmp_top[sub_indx4top[j]], top_med_indx);
            SEstop(2*n_snps_top+j) = smrwk->seyz[slct_indx_tmp_top[sub_indx4top[j]]];
        }

        MatrixXd Ctop;
        if (ldmatrix)
        {
            make_XMat(bdata, Id4top_pruned, _X);
            cor_calc(Ctop, _X);
        }
        else Ctop.setIdentity(Id4top_pruned.size(), Id4top_pruned.size());

        // Calculate MR effects

        cout << "Calculating MR effects ..." << endl;

        // total exposure effect

        if (!exclude_top_SNP){
            if (ldmatrix) {    
                mr_ivw_LD(bzx_totmr, bzy_totmr, sezx_totmr, sezy_totmr, Ctot, b_tot_ivw, se_tot_ivw, p_tot_ivw, N_ref);
            } else mr_ivw(bzx_totmr, bzy_totmr, sezx_totmr, sezy_totmr, b_tot_ivw, se_tot_ivw, p_tot_ivw);       

            cout << "Total MR effects: b = " << b_tot_ivw << "; p = " << p_tot_ivw << endl;
        } else {
            if (ldmatrix) {    
                mr_ivw_LD(bzx_totmr_extop, bzy_totmr_extop, sezx_totmr_extop, sezy_totmr_extop, Ctot_extop, b_tot_ivw, se_tot_ivw, p_tot_ivw, N_ref);
            } else mr_ivw(bzx_totmr_extop, bzy_totmr_extop, sezx_totmr_extop, sezy_totmr_extop, b_tot_ivw, se_tot_ivw, p_tot_ivw);       

            cout << "Total MR effects (without top SNP): b = " << b_tot_ivw << "; p = " << p_tot_ivw << endl;
        }
        
        // Calculate the direct effect for top mediator

        VectorXd betas, vars;

        if (n_snps_top > 3){

            mvmr_ivw_LD(Xtop, yztop, SEstop, Ctop, betas, vars, N_ref);
            b_top_direct = betas(0);
            se_top_direct = sqrt(vars[0]);
            p_top_direct = 2*pnorm(fabs(b_top_direct/se_top_direct));

            cout << "Top-mediator mediation MR effects: b = " << b_top_direct << "; p = " << p_top_direct << endl;
        } else {
            cout << "Less than 3 instrumental SNPs were found for both the exposure and and top mediator.\nNo top-MR analysis conducted." << endl;
        }

        // Calculate the direct effect for cis-mediators

        n_snps_e = snp_indx_e.size();

        // Shrink non-significant effect sizes if required
        if (p_shrinkage < 1){
            MatrixXd X_E(n_snps_e, medNum_included + 1);
            for (int j=0; j< n_snps_e; j++){
                for (int k=0;k<(medNum_included+1);k++){
                    X_E(j,k) = X(snp_indx_e[j],k);
                }
            }
            double z_shrinkage = fabs(qnorm5(p_shrinkage/2)); 
            X = (X_z.array() < z_shrinkage).select(0,X);
            for (int j=0; j< n_snps_e; j++){
                for (int k=0;k<(medNum_included+1);k++){
                    X(snp_indx_e[j],k) = X_E(j,k);
                }
            }
        }
            
        betas.resize(0), vars.resize(0);

        if (!exclude_top_SNP){
            mvmr_ivw_LD(X, yz, SEs, Ccis, betas, vars, N_ref);
            b_direct = betas(0);
            se_direct = sqrt(vars[0]);
            p_direct = 2*pnorm(fabs(b_direct/se_direct));

            for (int k=0; k<medNum_included; k++){
                med_info->b_my.push_back(betas(k+1));
                med_info->se_my.push_back(sqrt(vars[k+1]));
                med_info->p_my.push_back(2*pnorm(fabs(betas(k+1)/sqrt(vars[k+1]))));
                med_info->n_my.push_back(snp_count4multi);
            }

            cout << "Cis-mediation MR effects: b = " << b_direct << "; p = " << p_direct << endl;
        } else {
            mvmr_ivw_LD(X_extop, yz_extop, SEs_extop, Ccis, betas, vars, N_ref);
            b_direct = betas(0);
            se_direct = sqrt(vars[0]);
            p_direct = 2*pnorm(fabs(b_direct/se_direct));

            for (int k=0; k<medNum_included; k++){
                med_info->b_my.push_back(betas(k+1));
                med_info->se_my.push_back(sqrt(vars[k+1]));
                med_info->p_my.push_back(2*pnorm(fabs(betas(k+1)/sqrt(vars[k+1]))));
                med_info->n_my.push_back(snp_count4multi);
            }

            cout << "Cis-mediation MR effects (without top SNP): b = " << b_direct << "; p = " << p_direct << endl;
        }

        FILE* snpfile;
        string snpfilename = string(outFileName) + "." + exponame + ".snps";

        if (get_snp_effects_flg){
            snpfile = fopen(snpfilename.c_str(), "w");
            if (!(snpfile)) {
                printf("ERROR: open error %s\n", snpfilename.c_str());
                exit(1);
            }
            snpstr="SNP\tA1\tA2\tbeta.out\tse.out\tbeta." + exponame + "\tse." + exponame;

            for (int k=0; k<medNum_included; k++){
                snpstr += "\tbeta." + med_probe_names[med_sub_indx[k]] + "\tse." + med_probe_names[med_sub_indx[k]]; 
            }

            snpstr += "\n";
            
            if(fputs_checked(snpstr.c_str(),snpfile))
            {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
            }

            for (int j=0; j<snp_count4multi; j++){
                snpstr = smrwk->rs[sub_indx4multi[j]] + '\t' + smrwk->allele1[sub_indx4multi[j]] + '\t' + smrwk->allele2[sub_indx4multi[j]] + '\t';
                snpstr += atos(smrwk->byz[sub_indx4multi[j]]) + '\t' + atos(smrwk->seyz[sub_indx4multi[j]]) + '\t' + atos(smrwk->bxz[sub_indx4multi[j]]) + '\t' + atos(smrwk->sexz[sub_indx4multi[j]]);

                for (int k=0; k<medNum_included; k++){
                    snpstr += "\t" + atos(smrwk->Xm_b(sub_indx4multi[j], med_sub_indx[k])) + '\t' + atos(smrwk->Xm_se(sub_indx4multi[j], med_sub_indx[k]));
                }
           
                snpstr += "\n"; 

                if(fputs_checked(snpstr.c_str(),snpfile))
                {
                    printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                    exit(EXIT_FAILURE);
                }
            }

            fclose(snpfile);
        }

        FILE* topsnpfile;
        string topsnpfilename = string(outFileName) + "." + exponame + ".top.snps";

        if (get_snp_effects_top_flg){
            topsnpfile = fopen(topsnpfilename.c_str(), "w");
            if (!(topsnpfile)) {
                printf("ERROR: open error %s\n", topsnpfilename.c_str());
                exit(1);
            }
            snpstr="SNP\tA1\tA2\tbeta.out\tse.out\tbeta." + exponame + "\tse." + exponame + "\tbeta." + med_probe_names[top_med_indx] + "\tse." + med_probe_names[top_med_indx] + "\n";
       
            if(fputs_checked(snpstr.c_str(),topsnpfile))
            {
                printf("ERROR: in writing file %s .\n", topsnpfilename.c_str());
                exit(EXIT_FAILURE);
            }

            for (int j=0; j<n_snps_top; j++){
                snpstr = smrwk->rs[slct_indx_tmp_top[sub_indx4top[j]]] + '\t' + smrwk->allele1[slct_indx_tmp_top[sub_indx4top[j]]] + '\t' + smrwk->allele2[slct_indx_tmp_top[sub_indx4top[j]]] + '\t';
                snpstr += atos(smrwk->byz[slct_indx_tmp_top[sub_indx4top[j]]]) + '\t' + atos(smrwk->seyz[slct_indx_tmp_top[sub_indx4top[j]]]) + '\t' + atos(smrwk->bxz[slct_indx_tmp_top[sub_indx4top[j]]]) + '\t';
                snpstr += atos(smrwk->sexz[slct_indx_tmp_top[sub_indx4top[j]]]) + '\t' + atos(smrwk->Xm_b(slct_indx_tmp_top[sub_indx4top[j]], top_med_indx)) + '\t';
                snpstr += atos(smrwk->Xm_se(slct_indx_tmp_top[sub_indx4top[j]], top_med_indx)) + '\n';

                if(fputs_checked(snpstr.c_str(),topsnpfile))
                {
                    printf("ERROR: in writing file %s .\n", topsnpfilename.c_str());
                    exit(EXIT_FAILURE);
                }
            }

            fclose(topsnpfile);
        }

        return snp_count4multi;
    }


    void smr_mediation_prbmatch(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName, char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool multi_corr_snp_flg, bool get_snp_effects_flg, int min_snp)
    {
        //here eqtlFileName is the outcome and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not supposed to be used together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not supposed to be used together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not supposed to be used together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not supposed to be used together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        eqtlInfo mdata;
        eqtlInfo esdata;
        bInfo bdata;
        gwasData gdata;
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if((eqtlFileName==NULL) || (eqtlFileName2==NULL)) throw("Error: please input exposure and mediator QTL summary data by the flag --eqtl-summary.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        read_esifile(&mdata, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&mdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&mdata, snplst2exclde);
        read_epifile(&mdata, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&mdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&mdata, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&mdata, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&mdata, oprobe);
        if(oproblst2exclde != NULL) exclude_prob(&mdata, oproblst2exclde);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&mdata, oprobe2rm);
        
        read_besdfile(&mdata, string(eqtlFileName)+".besd");
        if(mdata._rowid.empty() && mdata._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        //if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName); //no need here, already extracted in etrait
        //if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_gwas_data(&gdata, gwasFileName);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        //if(snplstName != NULL) extract_snp(&bdata, snplstName);
        //if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &mdata, &esdata, &gdata);
        // if no snp left after check
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &mdata, &esdata, &gdata);
        }
        if(forcefrqck)
        {
            double prop=freq_check(&bdata, &mdata, &esdata, &gdata, afthresh, percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
    
        //mdata is not updated, so from now on _esi_include should be used always.
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(snpproblstName) extract_targets(&esdata, snpproblstName,  prb_snp);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("ERROR: no data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }
        
        update_gwas(&gdata);

        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }

        string outstr="";
        string smrfile = string(outFileName)+".msmr";
        FILE* smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("ERROR: open error %s\n", smrfile.c_str());
            exit(1);
        }

        outstr="Expo_ID\tExpo_Name\tChr\tExpo_bp\tMed_ID\tMed_Name\tMed_bp\tb_ivw_tot\tse_ivw_tot\tp_ivw_tot\t";
        outstr+="b_direct\tse_direct\tp_direct\tb_xm\tse_xm\tp_xm\tb_my\tse_my\tp_my\t";
        outstr+="n_snps_tot\tn_snps_med\tn_snps_e\tn_xm\n";

        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }

        FILE* snpfile;
        string snpfilename = string(outFileName) + ".snps";
        string snpstr;
        if (get_snp_effects_flg){
            snpfile = fopen(snpfilename.c_str(), "w");
            if (!(snpfile)) {
                printf("ERROR: open error %s\n", snpfilename.c_str());
                exit(1);
            }
            snpstr="Exposure\tMediator\tSNP\tA1\tA2\tbeta.exp\tse.exp\tbeta.med\tse.med\tbeta.out\tse.out\n";
            if(fputs_checked(snpstr.c_str(),snpfile))
            {
                printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }
       
        vector<int> includebk=esdata._include;

        logstr= atos(esdata._esi_include.size())+ " SNPs are included in the analysis.\n";
        cout<<logstr<<endl;
        //double disp=0.0;

        cis_itvl=cis_itvl*1000;
        long write_count=0;
        SMRWKM smrwk;

        unsigned int probNum = mdata._probNum;

        float progr0=0.0 , progr1;
        progress_print(progr0);

        for( int ii=0;ii<probNum;ii++){

            progr1=1.0*ii/probNum;
                if(progr1-progr0-0.01>1e-6 || ii+1==probNum)
                {
                    if(ii+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }

            //progress(ii,disp,(int)etrait._probNum);
            //
            string medname=mdata._epi_prbID[ii];
            string medgene=mdata._epi_gene[ii];
            //printf("\nStarting with mediator [ ]...\n");
            cout<<"\nStarting with mediator [ " + medname + " ]..."<<endl;

            // find the matching probe to mediator qtl data
            int jj=-9;
            string exponame;
            string expogene;

            for (int j = 0; j<includebk.size(); j++){
                exponame=esdata._epi_prbID[includebk[j]];
                expogene=esdata._epi_gene[includebk[j]];

                if ((medgene==expogene) || (medgene==exponame) || (medname==expogene) || (medname==exponame)){
                    jj=includebk[j];
                    break;
                }

            }      

            if (jj==-9){
                cout<<"No matching exposure probe found for mediator [ " + medname + " ]..."<<endl;
                continue;
            }

            int expobp=esdata._epi_bp[jj];
            int expochr=esdata._epi_chr[jj];
            exponame=esdata._epi_prbID[jj];
            expogene=esdata._epi_gene[jj];
            cout<<"Matching exposure probe: [ " + exponame + " ]..."<<endl;

            gwasData meddata;
            meddata.allele_1.resize(mdata._esi_include.size());
            meddata.allele_2.resize(mdata._esi_include.size());
            meddata.byz.resize(mdata._esi_include.size());
            meddata.seyz.resize(mdata._esi_include.size());
            meddata.freq.resize(mdata._esi_include.size());
            meddata.pvalue.resize(mdata._esi_include.size());
            meddata.splSize.resize(mdata._esi_include.size());
            meddata.snpName.resize(mdata._esi_include.size());
            
            for(int j=0;j<mdata._esi_include.size();j++){
                meddata.seyz[j]=-9;
                meddata.pvalue[j]=-9;
            }
            meddata._include.clear();
            
            int medchr=mdata._epi_chr[ii];
            int medbp=mdata._epi_bp[ii];
       
            int count=0;
            if(mdata._rowid.empty())
            {
                for (int j = 0; j<mdata._esi_include.size(); j++)
                {
                    if (fabs(mdata._sexz[ii][mdata._esi_include[j]] + 9) > 1e-6)
                    {
                        meddata.byz[j]=mdata._bxz[ii][mdata._esi_include[j]];
                        meddata.seyz[j]=mdata._sexz[ii][mdata._esi_include[j]];
                        double z=mdata._bxz[ii][mdata._esi_include[j]]/mdata._sexz[ii][mdata._esi_include[j]];
                        double p=pchisq(z*z,1);
                        meddata.pvalue[j]=p;
                        meddata.snpName[j]=mdata._esi_rs[mdata._esi_include[j]];
                        meddata.allele_1[j]=mdata._esi_allele1[mdata._esi_include[j]];
                        meddata.allele_2[j]=mdata._esi_allele2[mdata._esi_include[j]];
                        meddata._include.push_back(mdata._esi_include[j]); // row id selected
                        count++;
                    }
                }
            }
            else
            {
                uint64_t beta_start=mdata._cols[ii<<1];
                uint64_t se_start=mdata._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=mdata._rowid[beta_start+j];
                    int idx=(int)(find(mdata._esi_include.begin(), mdata._esi_include.end(), ge_rowid)-mdata._esi_include.begin());
                    if(idx<mdata._esi_include.size())
                    {
                        meddata.byz[idx]=mdata._val[beta_start+j];
                        meddata.seyz[idx]=mdata._val[se_start+j];
                        double z=mdata._val[beta_start+j]/mdata._val[se_start+j];
                        double p=pchisq(z*z,1);
                        meddata.pvalue[idx]=p;
                        meddata.snpName[idx]=mdata._esi_rs[ge_rowid];
                        meddata.allele_1[idx]=mdata._esi_allele1[ge_rowid];
                        meddata.allele_2[idx]=mdata._esi_allele2[ge_rowid];
                        meddata._include.push_back(idx);
                        count++;
                    }
                    
                }
            }

            meddata.snpNum=meddata.snpName.size();

            // fill eqtl data
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=jj;
            fill_smr_wk(&bdata, &gdata, &esdata, &meddata, &smrwk, cis_itvl);      

            if (smrwk.bxz.size() == 0) {
                cout<<"WARNING: no SNP fetched for probe [ " + exponame + " ]."<<endl;
                continue;
            }
            cout<< atos(smrwk.bxz.size()) + " SNPs are included from the cis-region of the probe [ " + exponame + " ]."<<endl;

            // step2: get top-SNP
            printf("Checking the top-SNP in the region....\n");
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            VectorXd zsxz;
            zsxz=ei_bxz.array()/ei_sexz.array();
            Map<VectorXd> ei_bmz(&smrwk.bmz[0],smrwk.bmz.size());
            Map<VectorXd> ei_semz(&smrwk.semz[0],smrwk.semz.size());
            VectorXd zsmz;
            zsmz=ei_bmz.array()/ei_semz.array();
            int maxidx, maxidm;
            maxidx=max_abs_id(zsxz); 
            maxidm=max_abs_id(zsmz); 
            double pxz_val = pnorm(fabs(zsxz[maxidx]));
            double pmz_val = pnorm(fabs(zsmz[maxidm]));
            //double computing, consistency should be checked
            for(int tid=0;tid<zsxz.size();tid++) {
                if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                {
                    //printf("ERROR: zxz not consistent %f vs %f. please report this.\n",zsxz(tid),smrwk.zxz[tid]);
                    exit(EXIT_FAILURE);
                }
            }

            string topsnpnamex=smrwk.rs[maxidx];
            string topsnpnamem=smrwk.rs[maxidm];
            cout<< "The top SNP of exposure probe [ " + exponame + " ] is " + topsnpnamex + " with p-value " + atos(pxz_val) <<endl;
            cout<< "The top SNP of mediator probe [ " + medname + " ] is " + topsnpnamem + " with p-value " + atos(pmz_val) <<endl;
            if((pxz_val>p_smr) && (pmz_val > p_medsmr)){
                logstr="WARNING: no SNP passed the exposure p-value threshold " + atos(p_smr) + " or the mediator p-value threshold " + atos(p_medsmr) + ".\n";
                cout<<logstr<<endl;
                continue;
            }     

            vector<string> snp4msmr;
            double b_tot_ivw = -9, se_tot_ivw = -9, p_tot_ivw = -9, b_direct = -9, se_direct = -9, p_direct = -9, b_xm = -9, se_xm = -9, p_xm = -9, b_my = -9, se_my = -9, p_my = -9; 
            int n_snps_tot, n_snps_e, n_xm;

            int snp_count=mr_mediation_analysis(&bdata, &smrwk, b_tot_ivw, se_tot_ivw, p_tot_ivw, b_direct, se_direct, p_direct, b_xm, se_xm, p_xm, b_my, se_my, p_my, n_snps_tot, n_snps_e, n_xm, p_smr, p_medsmr, ld_top_multi, ldmatrix, N_ref, trev, snp4msmr, multi_corr_snp_flg, get_snp_effects_flg, min_snp, snpfile, snpfilename, exponame, medname);
            if(snp_count==-9){
                logstr="No mediation analysis performed because less than 3 SNPs passed the p-value or LD-square thresholds.";
                cout<<logstr<<endl;
                continue;
            }
            
            // output snp set list

            string setstr=exponame+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<snp4msmr.size();j++)
            {
                setstr=snp4msmr[j]+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            setstr="end\n";
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            // end of output

            outstr =  exponame + '\t' + expogene + '\t' + atos(medchr) + '\t' + atos(expobp) + '\t' + medname + '\t' + medgene + '\t' + atos(medbp) + '\t'; 
            outstr += atos(b_tot_ivw)  + '\t' + atos(se_tot_ivw) + '\t' + dtos(p_tot_ivw) + '\t';
            outstr += atos(b_direct)  + '\t' + atos(se_direct) + '\t' + dtos(p_direct) + '\t';
            outstr += atos(b_xm) + '\t' + atos(se_xm) + '\t' + dtos(p_xm) + '\t';
            outstr += atos(b_my) + '\t' + atos(se_my) + '\t' + dtos(p_my) + '\t';
            outstr += atos(n_snps_tot) + '\t' + atos(snp_count) + '\t' + atos(n_snps_e) + '\t' + atos(n_xm) + '\n';
            
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
            write_count++;
            
        }

        cout<<"\nMR-IVW mediation analysis on QTL exposure through QTL mediator data completed.\n Results of "<<write_count<<" probes have been saved in the file [" + smrfile + "]."<<endl;     

        fclose(smr);
        fclose(setlst);
        if (get_snp_effects_flg){
            fclose(snpfile);
        }
        
    }

    void smr_multi_mediation_cis(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl,char* traitlstName, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, bool cis2all, double ld_top_multi, bool ldmatrix, double trev, double afthresh,double percenthresh, bool multi_corr_snp_flg, bool get_snp_effects_flg, bool get_snp_effects_top_flg, int min_snp, bool incmp_expo, double med_R_thresh, bool uncorr_med, bool exclude_top_SNP, double p_shrinkage, double p_expo_med)
    {

        //here eqtlFileName is the outcome and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not supposed to be used together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not supposed to be used together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not supposed to be used together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not supposed to be used together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        eqtlInfo mdata;
        eqtlInfo esdata;
        bInfo bdata;
        gwasData gdata;
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if((eqtlFileName==NULL) || (eqtlFileName2==NULL)) throw("Error: please input exposure and mediator QTL summary data by the flag --eqtl-summary.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        read_esifile(&mdata, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&mdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&mdata, snplst2exclde);
        read_epifile(&mdata, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&mdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&mdata, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&mdata, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&mdata, oprobe);
        if(oproblst2exclde != NULL) exclude_prob(&mdata, oproblst2exclde);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&mdata, oprobe2rm);
        
        read_besdfile(&mdata, string(eqtlFileName)+".besd");
        if(mdata._rowid.empty() && mdata._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        read_gwas_data(&gdata, gwasFileName);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        //if(snplstName != NULL) extract_snp(&bdata, snplstName);
        //if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &mdata, &esdata, &gdata); 
        // if no snp left after check
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &mdata, &esdata, &gdata);
        }
        if(forcefrqck)
        {
            double prop=freq_check(&bdata, &mdata, &esdata, &gdata, afthresh, percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
    
        //mdata is not updated, so from now on _esi_include should be used always.
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("ERROR: no data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }
        
        update_gwas(&gdata);

        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }
        
        string expolstfile = string(outFileName)+".expoProbes.info";
        FILE* expolst=NULL;
        expolst = fopen(expolstfile.c_str(), "w");
        if (!(expolst)) {
            printf("Open error %s\n", expolstfile.c_str());
            exit(1);
        }

        string expostr="Expo_ID\tExpo_Name\tMediation_status\n";
        if(fputs_checked(expostr.c_str(),expolst))
        {
            printf("ERROR: in writing file %s .\n", expolstfile.c_str());
            exit(EXIT_FAILURE);
        }
        string mediation_status;

        string medlstfile = string(outFileName)+".medProbes.info";
        FILE* medlst=NULL;
        medlst = fopen(medlstfile.c_str(), "w");
        if (!(medlst)) {
            printf("Open error %s\n", medlstfile.c_str());
            exit(1);
        }

        string medstr="Expo_ID\tExpo_Name\tChr\tExpo_bp\tMed_ID\tMed_Name\tMed_bp\tb_xm\tse_xm\tp_xm\tn_xm\tb_my\tse_my\tp_my\tn_my\n";
        if(fputs_checked(medstr.c_str(),medlst))
        {
            printf("ERROR: in writing file %s .\n", medlstfile.c_str());
            exit(EXIT_FAILURE);
        }

        string outstr="";
        string smrfile = string(outFileName)+".msmr";
        FILE* smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("ERROR: open error %s\n", smrfile.c_str());
            exit(1);
        }
        
        outstr="Expo_ID\tExpo_Name\tChr\tExpo_bp\tTop_Med_ID\tTop_Med_Name\tb_ivw_tot\tse_ivw_tot\tp_ivw_tot\t";
        outstr+="b_direct_top\tse_direct_top\tp_direct_top\tb_direct_cis\tse_direct_cis\tp_direct_cis\tn_snps_tot\t";
        outstr+="n_snps_top\tn_snps_cis\tn_snps_e\tnum_med_cis\tnum_med_cis_present\n";
        
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }
       
        vector<int> includebk=esdata._include;
        vector<int> includemed=mdata._include;

        logstr= atos(esdata._esi_include.size())+ " SNPs are included in the analysis.\n";
        cout<<logstr<<endl;
        //double disp=0.0;

        cis_itvl=cis_itvl*1000;
        long write_count=0;
        SMRWKMULT smrwk;
        MED_info med_info;

        unsigned int probNum = esdata._include.size();
        int medNum = mdata._include.size();

        float progr0=0.0 , progr1;
        progress_print(progr0);

        for( int ii_r=0;ii_r<probNum;ii_r++){

            progr1=1.0*ii_r/probNum;
                if(progr1-progr0-0.01>1e-6 || ii_r+1==probNum)
                {
                    if(ii_r+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }
            
            int ii = includebk[ii_r];

            string exponame=esdata._epi_prbID[ii];
            string expogene=esdata._epi_gene[ii];
        
            int expochr=esdata._epi_chr[ii];
            int expobp=esdata._epi_bp[ii];

            cout<<"\nStarting with exposure [ " + exponame + " ]..."<<endl;

            // find mediator probes within a window of cis_itvl

            vector<string> med_probe_names;
            vector<string> med_probe_genes;
            vector<int> med_probe_bp;
            vector<uint32_t> med_probe_idx;
            int medchr, medbp;

            for (int j = 0; j<medNum; j++){
                medchr=mdata._epi_chr[includemed[j]];
                medbp=mdata._epi_bp[includemed[j]];

                if ((medchr == expochr) && (abs(medbp - expobp) <= cis_itvl)){
                    med_probe_names.push_back(mdata._epi_prbID[includemed[j]]);
                    med_probe_genes.push_back(mdata._epi_gene[includemed[j]]);
                    med_probe_bp.push_back(medbp);
                    med_probe_idx.push_back(includemed[j]);
                }
            }
            
            int medNum_cis = med_probe_idx.size();

            if (medNum_cis < 1){
                cout<<"No mediator probes found for [ " + exponame + " ]..."<<endl;
                mediation_status = "No mediator in the cis-region.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            } else {
                cout<< atos(medNum_cis) + " mediator probes found for [ " + exponame + " ]..."<<endl;
            }
 
            MatrixXd Xm_b_raw(mdata._esi_include.size(), medNum_cis), Xm_se_raw(mdata._esi_include.size(), medNum_cis);

            for(int j=0;j<mdata._esi_include.size();j++){
                for (int k=0;k<medNum_cis;k++){
                    Xm_b_raw(j,k) = -9;
                    Xm_se_raw(j,k) = -9;
                }
            }

            int jj;

            for (int k=0;k<medNum_cis;k++){

                jj = med_probe_idx[k];

                if(mdata._rowid.empty()){

                    for (int j = 0; j<mdata._esi_include.size(); j++){
                        if (fabs(mdata._sexz[jj][mdata._esi_include[j]] + 9) > 1e-6){
                            Xm_b_raw(j,k) = mdata._bxz[jj][mdata._esi_include[j]];
                            Xm_se_raw(j,k) = mdata._sexz[jj][mdata._esi_include[j]];
                        }             
                    }
                }
                else {
                    uint64_t beta_start=mdata._cols[jj<<1];
                    uint64_t se_start=mdata._cols[1+(jj<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=mdata._rowid[beta_start+j];
                        int idx=(int)(find(mdata._esi_include.begin(), mdata._esi_include.end(), ge_rowid)-mdata._esi_include.begin());
                        if(idx<mdata._esi_include.size())
                        {
                            Xm_b_raw(idx,k) = mdata._val[beta_start+j];
                            Xm_se_raw(idx,k) = mdata._val[se_start+j];
                        }                    
                    }
                }
            }

            // fill data, SNPs within 2*cis_itvl are included
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=ii;
            smrwk.medNum_cis = medNum_cis;
            init_med_info(&med_info);
            bool in_med_cis;

            if(esdata._rowid.empty()){
                for (int j = 0; j<esdata._esi_include.size(); j++){

                    int snpbp=esdata._esi_bp[j];
                    int snpchr=esdata._esi_chr[j];

                    in_med_cis = true;
                    for (int k=0;k<medNum_cis;k++){
                        if (Xm_se_raw(j,k) < -8){
                            in_med_cis = false;
                        }
                    }

                    if(snpchr==esdata._epi_chr[ii] && abs(esdata._epi_bp[ii]-snpbp)<=2*cis_itvl && gdata.seyz[j]+9>1e-6 && in_med_cis)
                    {
                        smrwk.bxz.push_back(esdata._bxz[ii][j]);
                        smrwk.sexz.push_back(esdata._sexz[ii][j]);
                        smrwk.zxz.push_back(esdata._bxz[ii][j]/esdata._sexz[ii][j]);
                        smrwk.byz.push_back(gdata.byz[j]);
                        smrwk.seyz.push_back(gdata.seyz[j]);
                        smrwk.pyz.push_back(gdata.pvalue[j]);
                        smrwk.curId.push_back(j); //save snp id of the raw datastruct
                        smrwk.rs.push_back(esdata._esi_rs[j]);
                        smrwk.snpchrom.push_back(esdata._esi_chr[j]);
                        smrwk.allele1.push_back(esdata._esi_allele1[j]);
                        smrwk.allele2.push_back(esdata._esi_allele2[j]);
                        smrwk.bpsnp.push_back(esdata._esi_bp[j]);
                        smrwk.freq.push_back(bdata._mu[bdata._include[j]] / 2);
                    }
                }
            } 
            else {

                uint64_t beta_start=esdata._cols[ii<<1];
                uint64_t se_start=esdata._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=esdata._rowid[beta_start+j];
                    int snpbp=esdata._esi_bp[ge_rowid];
                    int snpchr=esdata._esi_chr[ge_rowid];

                    in_med_cis = true;
                    for (int k=0;k<medNum_cis;k++){
                        if (Xm_se_raw(ge_rowid,k) < -8){
                            in_med_cis = false;
                        }
                    }

                    if(snpchr==esdata._epi_chr[ii] && abs(esdata._epi_bp[ii]-snpbp)<=2*cis_itvl && gdata.seyz[ge_rowid]+9>1e-6 && in_med_cis)
                    {
                        smrwk.bxz.push_back(esdata._val[beta_start+j]);
                        smrwk.sexz.push_back(esdata._val[se_start+j]);
                        smrwk.zxz.push_back(esdata._val[beta_start+j]/esdata._val[se_start+j]);
                        smrwk.byz.push_back(gdata.byz[ge_rowid]);
                        smrwk.seyz.push_back(gdata.seyz[ge_rowid]);
                        smrwk.pyz.push_back(gdata.pvalue[ge_rowid]);
                        smrwk.curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk.rs.push_back(esdata._esi_rs[ge_rowid]);
                        smrwk.snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                        smrwk.allele1.push_back(esdata._esi_allele1[ge_rowid]);
                        smrwk.allele2.push_back(esdata._esi_allele2[ge_rowid]);
                        smrwk.bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                        smrwk.freq.push_back(bdata._mu[bdata._include[ge_rowid]] / 2);
                    }
                }
            }

            int num_cis_snps = smrwk.bxz.size();

            if (num_cis_snps == 0) {
                cout<<"WARNING: no SNP fetched for probe [ " + exponame + " ]." <<endl;
                mediation_status = "No SNPs fetched for exposure probe.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }

            smrwk.Xm_b.resize(num_cis_snps, medNum_cis);
            smrwk.Xm_se.resize(num_cis_snps, medNum_cis);
        
            for (int k=0;k<medNum_cis;k++){
                for (int j=0;j<num_cis_snps;j++){
                    smrwk.Xm_b(j, k) =  Xm_b_raw(smrwk.curId[j], k);
                    smrwk.Xm_se(j, k) =  Xm_se_raw(smrwk.curId[j], k);
                }
            }

            cout<< atos(num_cis_snps) + " SNPs are included from the cis-region of the probe [ " + exponame + " ]."<<endl;
            
            // step2: get top-SNP
            printf("Checking the top-SNP in the region....\n");
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            VectorXd zsxz;
            zsxz=ei_bxz.array()/ei_sexz.array();
            
            int maxidx;
            maxidx=max_abs_id(zsxz); 
            double pxz_val = pnorm(fabs(zsxz[maxidx]));

            //double computing, consistency should be checked
            for(int tid=0;tid<zsxz.size();tid++) {
                if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                {
                    exit(EXIT_FAILURE);
                }
            }

            string topsnpnamex=smrwk.rs[maxidx];
            cout<< "The top SNP of exposure probe [ " + exponame + " ] is " + topsnpnamex + " with p-value " + atos(pxz_val) <<endl;
           
           if((pxz_val>p_smr)){
                //printf("WARNING: no SNP passed the exposure p-value threshold %e or the mediator p-value threshold %e .\n", p_smr, p_medsmr);
                logstr="WARNING: no SNP passed the exposure p-value threshold " + atos(p_smr) + ".\n";
                cout<<logstr<<endl;
                mediation_status = "No SNP passed the exposure p-value.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }     

            vector<string> snp4msmr;
            double b_tot_ivw = -9, se_tot_ivw = -9, p_tot_ivw = -9, b_top_direct = -9, se_top_direct = -9, p_top_direct = -9, b_direct = -9, se_direct = -9, p_direct = -9;
            int medNum_included; 
            vector<int> med_sub_indx;
            int n_snps_tot = -9, n_snps_top = -9, top_med_indx = -9, n_snps_e = -9;
            mediation_status = "Mediation analysis conducted.";

            int snp_count=mr_mediation_analysis_cis(&bdata, &smrwk, b_tot_ivw, se_tot_ivw, p_tot_ivw, b_top_direct, se_top_direct, p_top_direct, b_direct, se_direct, p_direct, n_snps_tot, n_snps_top, n_snps_e, p_smr, p_medsmr, medNum_included, med_sub_indx, top_med_indx, ld_top_multi, ldmatrix, N_ref, trev, snp4msmr, multi_corr_snp_flg, get_snp_effects_flg, get_snp_effects_top_flg, outFileName, min_snp, exponame, med_probe_names, med_probe_genes, med_probe_bp, &med_info, med_R_thresh, uncorr_med, exclude_top_SNP, p_shrinkage, p_expo_med, mediation_status);

            expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
            if(fputs_checked(expostr.c_str(),expolst))
            {
                printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                exit(EXIT_FAILURE);
            }
            if(snp_count==-9) continue;
            
            // output snp set list

            string setstr=exponame+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<snp4msmr.size();j++)
            {
                setstr=snp4msmr[j]+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            setstr="end\n";
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            // end of output

            //output mediator information

            for(int k=0;k<medNum_included;k++)
            {
                medstr = exponame + '\t' + expogene + '\t' + atos(expochr) + '\t' + atos(expobp) + '\t';
                medstr += med_info.Med_ID[k] + '\t' + med_info.Gene_ID[k] + '\t' + atos(med_info.Med_bp[k]) + '\t'; 
                medstr += atos(med_info.b_xm[k]) + '\t' + atos(med_info.se_xm[k]) + '\t' + atos(med_info.p_xm[k]) + '\t' + atos(med_info.n_xm[k]) + '\t';
                medstr += atos(med_info.b_my[k]) + '\t' + atos(med_info.se_my[k]) + '\t' + atos(med_info.p_my[k]) + '\t' + atos(med_info.n_my[k]) + '\n';
                if(fputs_checked(medstr.c_str(),medlst))
                {
                    printf("ERROR: in writing file %s .\n", medlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }

            // end of output
            if (exclude_top_SNP){
                n_snps_tot = n_snps_tot-1;
                snp_count = snp_count-1;
                n_snps_e = n_snps_e-1;
            }        

            outstr =  exponame + '\t' + expogene + '\t' + atos(expochr) + '\t' + atos(expobp) + '\t'; 
            outstr += med_probe_names[top_med_indx] + '\t' + med_probe_genes[top_med_indx] + '\t';
            outstr += atos(b_tot_ivw)  + '\t' + atos(se_tot_ivw) + '\t' + atos(p_tot_ivw) + '\t';
            outstr += atos(b_top_direct)  + '\t' + atos(se_top_direct) + '\t' + atos(p_top_direct) + '\t';
            outstr += atos(b_direct)  + '\t' + atos(se_direct) + '\t' + atos(p_direct) + '\t';
            outstr += atos(n_snps_tot) + '\t' + atos(n_snps_top) + '\t' + atos(snp_count) + '\t' + atos(n_snps_e) + '\t';
            outstr += atos(medNum_included) + '\t' + atos(medNum_cis) + '\n';
            
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
            write_count++;
            
        }

        cout<<"\nMR analysis on QTL risk and mediator data completed.\n Results of "<<write_count<<" probes have been saved in the file [" + smrfile + "]."<<endl;     

        fclose(smr);
        fclose(setlst);
        fclose(medlst);
        fclose(expolst);
    }



    void smr_multi_mediation_trans(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, char* gwasFileName, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde, char* eproblst2exclde, double p_smr, double p_medsmr, int cis_itvl, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double ld_top_multi, bool ldmatrix, double trev, char* snpproblstName,double afthresh,double percenthresh, bool get_snp_effects_flg, int min_snp, bool incmp_expo, char* core_med_lst, char* qtl_med_lst, char* med_corr_snps_lst, int num_med, double med_R_thresh, bool uncorr_med, double p_shrinkage, double p_expo_med)
    {
        cout << "Starting MR trans-mediation. Be patient and allocate enough resources." << endl;
        //here eqtlFileName is the mediator and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not supposed to be used together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not supposed to be used together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not supposed to be used together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not supposed to be used together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        eqtlInfo mdata;
        eqtlInfo esdata;
        bInfo bdata;
        gwasData gdata;
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if((eqtlFileName==NULL) | (eqtlFileName2==NULL)) throw("Error: please input exposure and mediator QTL summary data by the flag --eqtl-summary.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(core_med_lst==NULL) throw("Error: please input core mediator file by the flag --core-mediators.");
        if(qtl_med_lst==NULL) throw("Error: please input mediator-QTL file by the flag --mediator-qtls.");
        if(med_corr_snps_lst==NULL && uncorr_med) throw("Error: please input file with SNPs to calculate the correlation among core mediators by the flag --mediator-corr-snps.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        
        // first read the exposure QTL data
        cout<<"Reading exposure QTL summary data..."<<endl;
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        //if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName); 
        //if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }        
        read_esifile(&mdata, string(eqtlFileName)+".esi");
        read_gwas_data(&gdata, gwasFileName);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        int N_ref = bdata._indi_num;
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        if (!incmp_expo){
            cout << "WARNING: Currently, the analysis does not support --incomplete-expo to be inactive. --incomplete-expo will be enabled." << endl;
            incmp_expo = true;
        }
        allele_check(&bdata, &mdata, &esdata, &gdata, incmp_expo); 
        // if no snp left after check
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &mdata, &esdata, &gdata);
        }
        if(forcefrqck)
        {
            double prop=freq_check(&bdata, &mdata, &esdata, &gdata, afthresh, percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
    
        //mdata is not updated, so from now on _esi_include should be used always.
        cout<<"Reading mediator QTL summary data..."<<endl;
        read_epifile(&mdata, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&mdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&mdata, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&mdata, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&mdata, oprobe);
        //if(snpproblstName) extract_targets(&esdata, snpproblstName,  prb_snp);
        if(oproblst2exclde != NULL) exclude_prob(&mdata, oproblst2exclde);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&mdata, oprobe2rm);
        
        /* mdata besd is not fully read into memory 
        read_besdfile(&mdata, string(eqtlFileName)+".besd");
        if(mdata._rowid.empty() && mdata._bxz.empty())
        {
            printf("ERROR: no data included in the analysis.\n");
            exit(EXIT_FAILURE);
        } */

        uint64_t m_epinum = mdata._probNum;
        uint64_t m_esinum = mdata._snpNum;
        
        update_gwas(&gdata);

        vector<string> core_mediators;
        read_msglist(core_med_lst, core_mediators, "core mediators");
        cout << core_mediators.size() << " core mediators were found." << endl;

        vector<string> med_qtls;
        read_msglist(qtl_med_lst, med_qtls, "mediator qtls");
        cout << med_qtls.size() << " mediator qtls were found." << endl;

        // intersection with included mediator SNPs
        vector<string> med_qtls_name;
        vector<uint64_t> med_qtls_included;
        vector<uint32_t> med_qtls_curId;
        vector<string> trans_rsvec, trans_allele1vec, trans_allele2vec;

        uint64_t tmp_idx; 
        vector<int>::iterator cur_iter;
        for (int i = 0; i < med_qtls.size(); i++){
            if (mdata._snp_name_map.find(med_qtls[i]) != mdata._snp_name_map.end()){
                tmp_idx = (uint64_t)mdata._snp_name_map[med_qtls[i]];
                cur_iter = find(mdata._esi_include.begin(), mdata._esi_include.end(), tmp_idx);
                if (cur_iter != mdata._esi_include.end()){
                    med_qtls_name.push_back(med_qtls[i]);
                    med_qtls_included.push_back((uint64_t)tmp_idx);
                    med_qtls_curId.push_back((uint32_t)(cur_iter - mdata._esi_include.begin()));
                    trans_rsvec.push_back(med_qtls[i]);
                    trans_allele1vec.push_back(mdata._esi_allele1[tmp_idx]);
                    trans_allele2vec.push_back(mdata._esi_allele2[tmp_idx]);
                }                
            } 
        } 

        cout << "Of these, " << med_qtls_included.size() << " were found in the mediator QTL dataset." << endl;

        // read SNP list to calculate correlation between core-mediators

        FILE * fptr;
        string besd_file = string(eqtlFileName)+".besd";
        fptr = fopen(besd_file.c_str(),"rb");
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        int infoLen=sizeof(uint32_t);
        if(indicator==DENSE_FILE_TYPE_3) infoLen=RESERVEDUNITS*sizeof(int);

        if(indicator!=DENSE_FILE_TYPE_1 && indicator!=DENSE_FILE_TYPE_3) {
            cout << "ERROR: the mediator BESD file is not in dense BESD format." << endl;
            exit(EXIT_FAILURE);
        }
        float* tmp=(float*)malloc(sizeof(float)*m_esinum<<1);

        MatrixXd C_core_med;

        if (uncorr_med){

            vector<string> med_corr_snps;
            read_msglist(med_corr_snps_lst, med_corr_snps, "mediator correlation-qtls");
            cout << med_corr_snps.size() << " mediator correlation-qtls were found." << endl;

            // calculate the correlation among core-mediators

            cout << "Calculating the correlation among core-mediators..." << endl;

            vector<uint64_t> med_corr_snps_idx, core_mediators_idx;

            for (int i = 0; i < med_corr_snps.size(); i++) med_corr_snps_idx.push_back((uint64_t)mdata._snp_name_map[med_corr_snps[i]]);
            for (int i = 0; i < core_mediators.size(); i++) core_mediators_idx.push_back((uint64_t)mdata._probe_name_map[core_mediators[i]]);

            // get the effect sizes
            MatrixXd X_mat(med_corr_snps_idx.size(), core_mediators_idx.size());

            for (int k = 0; k < core_mediators_idx.size(); k++){

                fseek(fptr,((core_mediators_idx[k]*m_esinum)<<3)+infoLen, SEEK_SET);
                fread(tmp, sizeof(float), m_esinum*2,fptr);

                for (int j = 0; j< med_corr_snps_idx.size(); j++){
                    X_mat(j,k) = tmp[med_corr_snps_idx[j]];
                }
            }
            cor_calc(C_core_med, X_mat); 
        }       

        // open output files

        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }

        string medlstfile = string(outFileName)+".medProbes.info";
        FILE* medlst=NULL;
        medlst = fopen(medlstfile.c_str(), "w");
        if (!(medlst)) {
            printf("Open error %s\n", medlstfile.c_str());
            exit(1);
        }

        string medstr="Expo_ID\tExpo_Name\tExpo_Chr\tExpo_bp\tMed_ID\tMed_Name\tMed_Chr\tMed_bp\tb_xm\tse_xm\tp_xm\tn_xm\tb_my\tse_my\tp_my\tn_my\n";
        if(fputs_checked(medstr.c_str(),medlst))
        {
            printf("ERROR: in writing file %s .\n", medlstfile.c_str());
            exit(EXIT_FAILURE);
        }

        string expolstfile = string(outFileName)+".expoProbes.info";
        FILE* expolst=NULL;
        expolst = fopen(expolstfile.c_str(), "w");
        if (!(expolst)) {
            printf("Open error %s\n", expolstfile.c_str());
            exit(1);
        }

        string expostr="Expo_ID\tExpo_Name\tMediation_status\n";
        if(fputs_checked(expostr.c_str(),expolst))
        {
            printf("ERROR: in writing file %s .\n", expolstfile.c_str());
            exit(EXIT_FAILURE);
        }
        string mediation_status;

        string outstr="";
        string smrfile = string(outFileName)+".msmr";
        FILE* smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("ERROR: open error %s\n", smrfile.c_str());
            exit(1);
        }

        outstr="Expo_ID\tExpo_Name\tFoc_Med_ID\tFoc_Med_Name\tChr\tExpo_bp\tb_ivw_tot\tse_ivw_tot\tp_ivw_tot\t";
        outstr+="b_direct_foc\tse_direct_foc\tp_direct_foc\tb_direct_trans\tse_direct_trans\tp_direct_trans\t";
        outstr+="n_snps_tot\tn_snps_foc\tn_snps_trans\tn_snps_e\tnum_med_trans\n";
        
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }
       
        vector<int> includeexp=esdata._include;
        vector<int> includemed=mdata._include;

        logstr= atos(esdata._esi_include.size())+ " SNPs are included in the analysis.\n";
        cout<<logstr<<endl;
        //double disp=0.0;

        cis_itvl=cis_itvl*1000;
        long write_count=0;
        SMRWKMULTTRANS smrwk;
        MED_info med_info;

        unsigned int probNum = esdata._include.size();
        int medNum = mdata._include.size();

        float progr0=0.0 , progr1;
        progress_print(progr0);

        for (int ii_r=0;ii_r<probNum;ii_r++){

            progr1=1.0*ii_r/probNum;
                if(progr1-progr0-0.01>1e-6 || ii_r+1==probNum)
                {
                    if(ii_r+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }
            
            int ii = includeexp[ii_r];

            string exponame=esdata._epi_prbID[ii];
            string expogene=esdata._epi_gene[ii];
        
            int expochr=esdata._epi_chr[ii];
            int expobp=esdata._epi_bp[ii];

            cout<<"\nStarting with exposure [ " + exponame + " ]..."<<endl;

            init_smr_wk(&smrwk);
            init_med_info(&med_info);
            smrwk.cur_prbidx=ii;

            // step1: get cis-exposure SNPs (e-SNPs); they should be in mediator SNP set as well
            printf("\nInitiating the cis-region of the exposure probe....\n", exponame.c_str());

            if(esdata._rowid.empty())
            {
                for (int j_r = 0; j_r<esdata._esi_include.size() ; j_r++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    if ((esdata._esi_include[j_r] > -9) && (mdata._esi_include[j_r] > -9)){
                        int j = esdata._esi_include[j_r];
                        if (fabs(esdata._sexz[ii][j] + 9) > 1e-6)
                        {
                            int snpbp=esdata._esi_bp[j];
                            int snpchr=esdata._esi_chr[j];
                            if(snpchr==expochr && fabs(expobp-snpbp)<=cis_itvl && gdata.seyz[j_r]+9>1e-6)
                            {
                                smrwk.bxz.push_back(esdata._bxz[ii][j]);
                                smrwk.sexz.push_back(esdata._sexz[ii][j]);
                                smrwk.zxz.push_back(fabs(esdata._bxz[ii][j]/esdata._sexz[ii][j]));
                                smrwk.byz.push_back(gdata.byz[j_r]);
                                smrwk.seyz.push_back(gdata.seyz[j_r]);
                                smrwk.pyz.push_back(gdata.pvalue[j_r]);
                                smrwk.curId.push_back(j_r);
                                smrwk.rs.push_back(esdata._esi_rs[j]);
                                smrwk.snpchrom.push_back(esdata._esi_chr[j]);
                                smrwk.allele1.push_back(esdata._esi_allele1[j]);
                                smrwk.allele2.push_back(esdata._esi_allele2[j]);
                                smrwk.bpsnp.push_back(esdata._esi_bp[j]);
                                smrwk.freq.push_back(bdata._mu[bdata._include[j_r]] / 2); 
                            }
                        }
                    }
                }   
            }
            else{
                uint64_t beta_start=esdata._cols[ii<<1];
                uint64_t se_start=esdata._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=esdata._rowid[beta_start+j];
                    uint32_t idx=(uint32_t)(find(esdata._esi_include.begin(), esdata._esi_include.end(), ge_rowid)-esdata._esi_include.begin());
                    if (idx < esdata._esi_include.size()){
                        int snpbp=esdata._esi_bp[ge_rowid];
                        int snpchr=esdata._esi_chr[ge_rowid];
                        if((snpchr==expochr) && (abs(expobp-snpbp)<=cis_itvl) && (mdata._esi_include[idx] > -9))
                        {
                            smrwk.bxz.push_back(esdata._val[beta_start+j]);
                            smrwk.sexz.push_back(esdata._val[se_start+j]);
                            smrwk.zxz.push_back(fabs(esdata._val[beta_start+j]/esdata._val[se_start+j]));
                            smrwk.byz.push_back(gdata.byz[idx]);
                            smrwk.seyz.push_back(gdata.seyz[idx]);
                            smrwk.pyz.push_back(gdata.pvalue[idx]);
                            smrwk.curId.push_back(idx);
                            smrwk.rs.push_back(esdata._esi_rs[ge_rowid]);
                            smrwk.snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                            smrwk.allele1.push_back(esdata._esi_allele1[ge_rowid]);
                            smrwk.allele2.push_back(esdata._esi_allele2[ge_rowid]);
                            smrwk.bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                            smrwk.freq.push_back(bdata._mu[bdata._include[idx]] / 2);
                        }
                    }
                }
            }

            int num_cisSNPs = smrwk.curId.size();
            if (num_cisSNPs < 3){
                cout << "Less than 3 cis-exposure SNPs. No MR analysis conducted." << endl;
                mediation_status = "Less than 3 significant E-SNPs for total MR calculation.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }

            // determine significant cis-exposure QTLs

            vector<uint32_t> sigexpoId;
            vector<double> sigexpobxz, sigexposexz, sigexpozxz, sigexpobyz, sigexposeyz;
            vector<string> sigsnpnames;
            vector<int> sigexpobpsnp;
   
            double z_smr= fabs(qnorm5(p_smr/2));

            for(int j=0;j<smrwk.curId.size();j++)
            {
                double ztmp = smrwk.zxz[j]; // zxz is in absolute value
                double zrev = (fabs(smrwk.bxz[j]) - fabs(smrwk.byz[j]))/sqrt(pow(smrwk.sexz[j],2) + pow(smrwk.seyz[j],2));
                if((ztmp>=z_smr) && (zrev>trev))
                {
                    sigexpoId.push_back(smrwk.curId[j]);
                    sigexpobxz.push_back(smrwk.bxz[j]);
                    sigexposexz.push_back(smrwk.sexz[j]);
                    sigexpozxz.push_back(ztmp); 
                    sigsnpnames.push_back(smrwk.rs[j]);
                    sigexpobyz.push_back(smrwk.byz[j]);
                    sigexposeyz.push_back(smrwk.seyz[j]);
                    sigexpobpsnp.push_back(smrwk.bpsnp[j]);
                }
            }
            int num_sigexpoSNPs=(int)sigexpoId.size();
            if(num_sigexpoSNPs<3) {
                cout << "Less than 3 significant cis-exposure SNPs. No MR analysis conducted." << endl;
                mediation_status = "Less than 3 significant E-SNPs for total MR calculation.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }

            // prune exposure SNPs
        
            vector<int> sub_indx;
            MatrixXd _X;
            make_XMat(&bdata, sigexpoId, _X); //_X: one row one individual, one column one SNP

            double sbat_ld_cutoff=sqrt(ld_top_multi);
            sbat_calcu_lambda(_X, num_sigexpoSNPs,  sbat_ld_cutoff, sub_indx, sigexpozxz); // num_sigexpoSNPs can change here
            cout << sub_indx.size() << " significant independent cis-exposure QTLs were found." << endl;

            int num_e_SNPs = sub_indx.size();
            if(num_e_SNPs<3) {
                cout << "Less than 3 independent significant cis-exposure SNPs. No MR analysis conducted." << endl;
                mediation_status = "Less than 3 E-SNPs.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }
            
            // get SNP names of independent SNPs & initialize and fill vectors for total MR calculation

            vector<uint32_t> e_Id;
            vector<string> e_snpnames;
            vector<int> e_bpsnp;

            VectorXd bzx_totmr(num_e_SNPs), sezx_totmr(num_e_SNPs), bzy_totmr(num_e_SNPs), sezy_totmr(num_e_SNPs);
            
            for (int j=0; j<sub_indx.size(); j++){
                e_Id.push_back(sigexpoId[sub_indx[j]]);
                e_snpnames.push_back(sigsnpnames[sub_indx[j]]);
                e_bpsnp.push_back(sigexpobpsnp[sub_indx[j]]);

                bzx_totmr(j) = sigexpobxz[sub_indx[j]];
                sezx_totmr(j) = sigexposexz[sub_indx[j]];
                bzy_totmr(j) = sigexpobyz[sub_indx[j]];
                sezy_totmr(j) = sigexposeyz[sub_indx[j]];
            }

            MatrixXd Ctot;
            if (ldmatrix)
            {
                make_XMat(&bdata, e_Id, _X);
                cor_calc(Ctot, _X);
            } else Ctot.setIdentity(num_e_SNPs, num_e_SNPs);

            // calculate default standard error for exposure SNPs (for missing values in mediation regression)
            double default_expo_se = sezx_totmr.mean();

            // step2: calculate correlations between exposure and core mediators using the e-SNPs

            // check if there is a focal mediator
            int jj=-9;
            smrwk.medFocal = false;
            smrwk.medFocal_core = false; // focal mediator that is also in the core mediator set

            string medname;
            string medgene;

            for (int j = 0; j<includemed.size(); j++){
                medname=mdata._epi_prbID[includemed[j]];
                medgene=mdata._epi_gene[includemed[j]];

                if ((medgene==expogene) || (medgene==exponame) || (medname==expogene) || (medname==exponame)){
                    jj=includemed[j];
                    break;
                }
            }      

            if (jj>-9){
                smrwk.medFocal = true;
                cout<<"Focal mediator found : [ " + medname + " ]"<<endl;
                if (find(core_mediators.begin(), core_mediators.end(), medname) != core_mediators.end()) smrwk.medFocal_core = true;
            } else {
                cout<<"No focal mediator found. Continuing with trans-mediators..." << endl;
            }           

            // get BESD indexes
            vector<uint64_t> e_snp_idx;
            vector<uint64_t> core_probe_idx;

            vector<string> cur_core_mediators; // in addition to core mediators, cur_core_mediatiors might contain a focal mediator
            vector<int> cur_core_idx;

            for (int k=0; k<core_mediators.size(); k++){
                cur_core_mediators.push_back(core_mediators[k]);
                cur_core_idx.push_back(k);
            }

            // put the focal mediator in first position
            if (smrwk.medFocal){
                if (smrwk.medFocal_core) {
                    int idx=(int)(find(cur_core_mediators.begin(), cur_core_mediators.end(), medname)-cur_core_mediators.begin());
                    swap(cur_core_mediators[0], cur_core_mediators[idx]);
                    swap(cur_core_idx[0], cur_core_idx[idx]);
                } else {
                    cur_core_mediators.insert(cur_core_mediators.begin(), medname);   
                    cur_core_idx.insert(cur_core_idx.begin(), cur_core_idx.size());      
                }
            } 

            vector<string> e_rsvec, e_allele1vec, e_allele2vec;

            for (int i = 0; i < e_snpnames.size(); i++){
                e_snp_idx.push_back((uint64_t)mdata._snp_name_map[e_snpnames[i]]);
                e_rsvec.push_back(e_snpnames[i]);
                e_allele1vec.push_back(mdata._esi_allele1[e_snp_idx[i]]);
                e_allele2vec.push_back(mdata._esi_allele2[e_snp_idx[i]]);
            }
            for (int i = 0; i < cur_core_mediators.size(); i++) core_probe_idx.push_back((uint64_t)mdata._probe_name_map[cur_core_mediators[i]]);

            // get the effect sizes
            MatrixXd e_Xm_b(e_snp_idx.size(), core_probe_idx.size());
            MatrixXd e_Xm_se(e_snp_idx.size(), core_probe_idx.size());
           
            fseek(fptr,0L,SEEK_SET);
            
            tmp=(float*)malloc(sizeof(float)*m_esinum<<1);
            cout << "Reading mediator QTL data of the independent cis-exposure QTLs ..." << endl;

            for (int k = 0; k < core_probe_idx.size(); k++){

                fseek(fptr,((core_probe_idx[k]*m_esinum)<<3)+infoLen, SEEK_SET);
                fread(tmp, sizeof(float), m_esinum*2,fptr);

                for (int j = 0; j< e_snp_idx.size(); j++){
                    e_Xm_b(j,k) = tmp[e_snp_idx[j]];
                    e_Xm_se(j,k) = tmp[e_snp_idx[j] + m_esinum];
                }
            }

            // calculate the correlations
            double cor_thresh = fabs(qnorm5(p_expo_med/2));
            double cor_tval, zrevm, b_xm, se_xm, p_xm, cor;
            vector<double> cor_med_val, cor_med_tval, secor_med_val; 
            vector<int> ncor_med_val;
            vector<uint64_t> core_probe_idx_cor, core_probe_idx_slct;
            vector<uint32_t> Id4em_all, Id4em_x; // only relevant if focal mediator with m-SNPs
            VectorXd expo_xz, expo_sexz, med_mz, med_semz;
            MatrixXd Cxm;
            int num_medPrbs = 0;
            int n_xm;

            vector<string> med_prbs_raw, med_prbs;
            vector<int> med_indx_raw, med_indx, cur_core_idx_slct, snp_e_slct;

            cout << "Calculate the causal effect of the exposure on the mediators ..." << endl;

            for (int k=0; k<cur_core_mediators.size(); k++){

                snp_e_slct.resize(0);

                if (smrwk.medFocal && (k == 0)){

                    for (int j=0; j<num_e_SNPs;j++){
                        zrevm = (fabs(bzx_totmr[j]) - fabs(e_Xm_b(j,k)))/sqrt(pow(sezx_totmr[j],2) + pow(e_Xm_se(j,k),2));
                        if (zrevm > trev){
                            snp_e_slct.push_back(j);
                            Id4em_all.push_back(e_Id[j]);
                            Id4em_x.push_back(e_Id[j]);
                        }
                    }                  
                    n_xm = snp_e_slct.size();
                    if (n_xm < 1){
                        b_xm = -9;
                        se_xm = -9;
                    } else {
                        expo_xz.resize(n_xm), med_mz.resize(n_xm), expo_sexz.resize(n_xm), med_semz.resize(n_xm);
                        for (int j=0; j<n_xm; j++){
                            expo_xz[j] = bzx_totmr[snp_e_slct[j]];
                            expo_sexz[j] = sezx_totmr[snp_e_slct[j]];
                            med_mz[j] = e_Xm_b(snp_e_slct[j],k);
                            med_semz[j] = e_Xm_se(snp_e_slct[j],k);
                        }
                        if (ldmatrix) {   
                            Cxm.resize(n_xm, n_xm);
                            for (int i=0; i< n_xm; i++){
                                for (int j=i; j< n_xm; j++){
                                    Cxm(i,j) = Ctot(snp_e_slct[i], snp_e_slct[j]);
                                    Cxm(j,i) = Cxm(i,j);
                                }
                            }
                            mr_ivw_LD(expo_xz, med_mz, expo_sexz, med_semz, Cxm, b_xm, se_xm, p_xm, N_ref);
                        } else mr_ivw(expo_xz, med_mz, expo_sexz, med_semz, b_xm, se_xm, p_xm);
                    }

                    num_medPrbs+=1;
                    med_prbs_raw.push_back(cur_core_mediators[k]);
                    med_indx_raw.push_back(k);
                    cur_core_idx_slct.push_back(cur_core_idx[k]);
                    cor_med_tval.push_back(1000); 
                    cor_med_val.push_back(b_xm); 
                    secor_med_val.push_back(se_xm);
                    ncor_med_val.push_back(n_xm);
                    core_probe_idx_cor.push_back((uint64_t)core_probe_idx[k]);
                    
                } else {

                    for (int j=0; j<num_e_SNPs;j++){
                        zrevm = (fabs(bzx_totmr[j]) - fabs(e_Xm_b(j,k)))/sqrt(pow(sezx_totmr[j],2) + pow(e_Xm_se(j,k),2));
                        if (zrevm > trev){
                            snp_e_slct.push_back(j);
                        }
                    }                  
                    n_xm = snp_e_slct.size();
                    if (n_xm < 1) continue;
                    
                    expo_xz.resize(n_xm), med_mz.resize(n_xm), expo_sexz.resize(n_xm), med_semz.resize(n_xm);
                    for (int j=0; j<n_xm; j++){
                        expo_xz[j] = bzx_totmr[snp_e_slct[j]];
                        expo_sexz[j] = sezx_totmr[snp_e_slct[j]];
                        med_mz[j] = e_Xm_b(snp_e_slct[j],k);
                        med_semz[j] = e_Xm_se(snp_e_slct[j],k);
                    }
                    if (ldmatrix) {   
                        Cxm.resize(n_xm, n_xm);
                        for (int i=0; i< n_xm; i++){
                            for (int j=i; j< n_xm; j++){
                                Cxm(i,j) = Ctot(snp_e_slct[i], snp_e_slct[j]);
                                Cxm(j,i) = Cxm(i,j);
                            }
                        }
                        mr_ivw_LD(expo_xz, med_mz, expo_sexz, med_semz, Cxm, b_xm, se_xm, p_xm, N_ref);
                    } else mr_ivw(expo_xz, med_mz, expo_sexz, med_semz, b_xm, se_xm, p_xm);

                    cor_tval = fabs(b_xm/se_xm);

                    if (cor_tval >= cor_thresh){

                        num_medPrbs+=1;
                        med_prbs_raw.push_back(cur_core_mediators[k]);
                        med_indx_raw.push_back(k);
                        cur_core_idx_slct.push_back(cur_core_idx[k]);
                        cor_med_tval.push_back(cor_tval); 
                        cor_med_val.push_back(b_xm); 
                        secor_med_val.push_back(se_xm);
                        ncor_med_val.push_back(n_xm);
                        core_probe_idx_cor.push_back((uint64_t)core_probe_idx[k]);
                    }
                }
            }

            cout << num_medPrbs << " mediator probes passed the correlation threshold and are included in the trans-mediation calculation." << endl;

            if(num_medPrbs<1) {
                cout << "Less than 1 mediator probe found. No MR analysis conducted." << endl;
                mediation_status = "No causally-associated mediator probe found.";
                expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
                if(fputs_checked(expostr.c_str(),expolst))
                {
                    printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }

            if ((num_med > -9) && (num_medPrbs > num_med) && !uncorr_med){
                cout << "Finding the top " << num_med << " correlated mediators ..." << endl;
                
                vector<int> sorted_corr = sort_indexes(cor_med_tval);
                reverse(sorted_corr.begin(), sorted_corr.end());

                for (int k=0; k < num_med; k++){
                    med_prbs.push_back(med_prbs_raw[sorted_corr[k]]);
                    med_indx.push_back(med_indx_raw[sorted_corr[k]]);
                    med_info.Med_ID.push_back(med_prbs_raw[sorted_corr[k]]);
                    med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[sorted_corr[k]]]);
                    med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[sorted_corr[k]]]);
                    med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[sorted_corr[k]]]);
                    med_info.b_xm.push_back(cor_med_val[sorted_corr[k]]);
                    med_info.se_xm.push_back(secor_med_val[sorted_corr[k]]);
                    med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[sorted_corr[k]]/secor_med_val[sorted_corr[k]])));
                    med_info.n_xm.push_back(ncor_med_val[sorted_corr[k]]);
                    core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[sorted_corr[k]]);
                }
                num_medPrbs = num_med;
            }
            else if (uncorr_med){
                vector<double> z_cor;
                MatrixXd C_med;
                if (smrwk.medFocal && !smrwk.medFocal_core){
                    C_med.resize(cur_core_idx_slct.size()-1, cur_core_idx_slct.size()-1);
                    for (int i=1; i<cur_core_idx_slct.size(); i++){
                        for (int j=1; j<cur_core_idx_slct.size(); j++){
                            C_med(i-1,j-1) = C_core_med(cur_core_idx_slct[i], cur_core_idx_slct[j]);
                        }
                        z_cor.push_back(cor_med_tval[i]);
                    }
                } else {
                    C_med.resize(cur_core_idx_slct.size(), cur_core_idx_slct.size());
                    for (int i=0; i<cur_core_idx_slct.size(); i++){
                        for (int j=0; j<cur_core_idx_slct.size(); j++){
                            C_med(i,j) = C_core_med(cur_core_idx_slct[i], cur_core_idx_slct[j]);
                        }
                        z_cor.push_back(cor_med_tval[i]);
                    }
                }

                // prune the mediators
                int m = z_cor.size();  
                int qi = 0;    
                vector<int> rm_ID1, med_sub_indx;
                rm_cor_sbat(C_med, med_R_thresh, m, rm_ID1, z_cor);
                //Create new index
                for (int i=0 ; i<m ; i++) {
                    if (rm_ID1.size() == 0){
                        med_sub_indx.push_back(i);
                    } else {
                        if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                        else med_sub_indx.push_back(i);
                    }
                }

                if (smrwk.medFocal && !smrwk.medFocal_core){
                    num_medPrbs = med_sub_indx.size()+1; 
                } else {
                    num_medPrbs = med_sub_indx.size();
                }
                
                cout << num_medPrbs << " uncorrelated mediators were found." << endl;

                // Select top N mediators if required
                if ((num_med > -9) && (num_medPrbs > num_med)){
                    cout << "Finding the top " << num_med << " uncorrelated mediators ..." << endl;
                    
                    vector<double> z_cor_slct;
                    for (int i=1; i<med_sub_indx.size(); i++) z_cor_slct.push_back(z_cor[med_sub_indx[i]]);
                    vector<int> sorted_corr = sort_indexes(z_cor_slct);
                    reverse(sorted_corr.begin(), sorted_corr.end());

                    if (smrwk.medFocal && !smrwk.medFocal_core){
                        med_prbs.push_back(med_prbs_raw[0]);
                        med_indx.push_back(med_indx_raw[0]);
                        med_info.Med_ID.push_back(med_prbs_raw[0]);
                        med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[0]]);
                        med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[0]]);
                        med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[0]]);
                        med_info.b_xm.push_back(cor_med_val[0]);
                        med_info.se_xm.push_back(secor_med_val[0]);
                        med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[0]/secor_med_val[0])));
                        med_info.n_xm.push_back(ncor_med_val[0]);
                        core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[0]);
                       
                        for (int k=0; k<num_med-1; k++){
                            med_prbs.push_back(med_prbs_raw[med_sub_indx[sorted_corr[k]]+1]);
                            med_indx.push_back(med_indx_raw[med_sub_indx[sorted_corr[k]]+1]);
                            med_info.Med_ID.push_back(med_prbs_raw[med_sub_indx[sorted_corr[k]]+1]);
                            med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]+1]]);
                            med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]+1]]);
                            med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]+1]]);
                            med_info.b_xm.push_back(cor_med_val[med_sub_indx[sorted_corr[k]]+1]);
                            med_info.se_xm.push_back(secor_med_val[med_sub_indx[sorted_corr[k]]+1]);
                            med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[med_sub_indx[sorted_corr[k]]+1]/secor_med_val[med_sub_indx[sorted_corr[k]]+1])));
                            med_info.n_xm.push_back(ncor_med_val[med_sub_indx[sorted_corr[k]]+1]);
                            core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]+1]);
                        }
                    } else {
                        for (int k=0; k<num_med; k++){
                            med_prbs.push_back(med_prbs_raw[med_sub_indx[sorted_corr[k]]]);
                            med_indx.push_back(med_indx_raw[med_sub_indx[sorted_corr[k]]]);
                            med_info.Med_ID.push_back(med_prbs_raw[med_sub_indx[sorted_corr[k]]]);
                            med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]]]);
                            med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]]]);
                            med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]]]);
                            med_info.b_xm.push_back(cor_med_val[med_sub_indx[sorted_corr[k]]]);
                            med_info.se_xm.push_back(secor_med_val[med_sub_indx[sorted_corr[k]]]);
                            med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[med_sub_indx[sorted_corr[k]]]/secor_med_val[med_sub_indx[sorted_corr[k]]])));
                            med_info.n_xm.push_back(ncor_med_val[med_sub_indx[sorted_corr[k]]]);
                            core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[med_sub_indx[sorted_corr[k]]]);
                        }
                    }
                    num_medPrbs = num_med;
                } else {
                    if (smrwk.medFocal && !smrwk.medFocal_core){
                        med_prbs.push_back(med_prbs_raw[0]);
                        med_indx.push_back(med_indx_raw[0]);
                        med_info.Med_ID.push_back(med_prbs_raw[0]);
                        med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[0]]);
                        med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[0]]);
                        med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[0]]);
                        med_info.b_xm.push_back(cor_med_val[0]);
                        med_info.se_xm.push_back(secor_med_val[0]);
                        med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[0]/secor_med_val[0])));
                        med_info.n_xm.push_back(ncor_med_val[0]);
                        core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[0]);
                        for (int k=0; k<num_medPrbs; k++){
                            med_prbs.push_back(med_prbs_raw[med_sub_indx[k]+1]);
                            med_indx.push_back(med_indx_raw[med_sub_indx[k]+1]);
                            med_info.Med_ID.push_back(med_prbs_raw[med_sub_indx[k]+1]);
                            med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[med_sub_indx[k]+1]]);
                            med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[med_sub_indx[k]+1]]);
                            med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[med_sub_indx[k]+1]]);
                            med_info.b_xm.push_back(cor_med_val[med_sub_indx[k]+1]);
                            med_info.se_xm.push_back(secor_med_val[med_sub_indx[k]+1]);
                            med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[med_sub_indx[k]+1]/secor_med_val[med_sub_indx[k]+1])));
                            med_info.n_xm.push_back(ncor_med_val[med_sub_indx[k]+1]);
                            core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[med_sub_indx[k]+1]);
                        }
                    } else {
                        for (int k=0; k<num_medPrbs; k++){
                            med_prbs.push_back(med_prbs_raw[med_sub_indx[k]]);
                            med_indx.push_back(med_indx_raw[med_sub_indx[k]]);
                            med_info.Med_ID.push_back(med_prbs_raw[med_sub_indx[k]]);
                            med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[med_sub_indx[k]]]);
                            med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[med_sub_indx[k]]]);
                            med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[med_sub_indx[k]]]);
                            med_info.b_xm.push_back(cor_med_val[med_sub_indx[k]]);
                            med_info.se_xm.push_back(secor_med_val[med_sub_indx[k]]);
                            med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[med_sub_indx[k]]/secor_med_val[med_sub_indx[k]])));
                            med_info.n_xm.push_back(ncor_med_val[med_sub_indx[k]]);
                            core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[med_sub_indx[k]]);
                        }
                    }
                }
            } else {
                for (int k=0; k< num_medPrbs; k++){
                    med_prbs.push_back(med_prbs_raw[k]);
                    med_indx.push_back(med_indx_raw[k]);
                    med_info.Med_ID.push_back(med_prbs_raw[k]);
                    med_info.Gene_ID.push_back(mdata._epi_gene[(int)core_probe_idx_cor[k]]);
                    med_info.Med_Chr.push_back(mdata._epi_chr[(int)core_probe_idx_cor[k]]);
                    med_info.Med_bp.push_back(mdata._epi_bp[(int)core_probe_idx_cor[k]]);
                    med_info.b_xm.push_back(cor_med_val[k]);
                    med_info.se_xm.push_back(secor_med_val[k]);
                    med_info.p_xm.push_back(2*pnorm(fabs(cor_med_val[k]/secor_med_val[k])));
                    med_info.n_xm.push_back(ncor_med_val[k]);
                    core_probe_idx_slct.push_back((uint64_t)core_probe_idx_cor[k]);
                }
            }

            // step3: get mediator QTLs and prune them

            cout << "Reading mediator QTLs ..." << endl;

            // get BESD indexes
            vector<uint64_t> med_probe_idx;

            for (int i=0; i < med_prbs.size(); i++) med_probe_idx.push_back((uint64_t)mdata._probe_name_map[med_prbs[i]]);

            // get the effect sizes
            MatrixXd trans_Xm_b(med_qtls_included.size(), med_probe_idx.size());
            MatrixXd trans_Xm_se(med_qtls_included.size(), med_probe_idx.size());
            MatrixXd trans_Xm_z(med_qtls_included.size(), med_probe_idx.size());

            for (int k=0; k < med_probe_idx.size(); k++){

                fseek(fptr,((med_probe_idx[k]*m_esinum)<<3)+infoLen, SEEK_SET);
                //fread(&tmp[0], sizeof(float), 2*m_esinum,fptr);
                fread(tmp, sizeof(float), 2*m_esinum,fptr);

                for (int j = 0; j< med_qtls_included.size(); j++){
                    trans_Xm_b(j,k) = tmp[med_qtls_included[j]];
                    trans_Xm_se(j,k) = tmp[med_qtls_included[j] + m_esinum];
                    trans_Xm_z(j,k) = fabs(trans_Xm_b(j,k)/trans_Xm_se(j,k));
                }
            }

            // select cis-qtls of selected mediators
            cout << "Select cis-QTLs ..." << endl;

            vector<int> trans_Xm_idx_slct_tmp;
            double z_m = fabs(qnorm5(p_medsmr/2));
            double ztmpm;
            int medchr, medbp, snpchr, snpbp;

            for (int j=0; j<med_qtls_included.size();j++){
                if (find(sigexpoId.begin(), sigexpoId.end(), med_qtls_curId[j]) != sigexpoId.end()) trans_Xm_idx_slct_tmp.push_back(j);
                else {
                    snpchr = mdata._esi_chr[(int)med_qtls_included[j]];
                    snpbp = mdata._esi_bp[(int)med_qtls_included[j]];
                    for (int k=0; k<num_medPrbs; k++){ 
                        ztmpm=trans_Xm_z(j,k);
                        medchr = mdata._epi_chr[(int)core_probe_idx_slct[k]];
                        medbp = mdata._epi_bp[(int)core_probe_idx_slct[k]];
                        if ((ztmpm >= z_m) && (medchr == snpchr) && (abs(medbp - snpbp) <= cis_itvl)){
                            trans_Xm_idx_slct_tmp.push_back(j);
                            break;
                        }
                    }
                }
            }

            // filter out SNPs with pleiotropic effects
            cout << "Remove pleiotropic cis-QTLs ..." << endl;
            
            vector<int> trans_Xm_idx_slct;
            vector<int> e_snp_Xm_idx_slct; // idx of trans-SNPs that are also independent e-SNPs (and not focal m-SNPs)
            vector<int> m_snp_Xm_idx_slct; // idx of trans-SNPs that are focal m-SNPs
            vector<double> zmz_foc_m; // only relevant if focal mediator with m-SNPs
            vector<int> trans_Xm_chr, trans_Xm_bp;
            vector<int> rm_eID;
            vector<uint32_t> trans_curId_slct, Id4em_m; // Id4em_m ... : only relevant if focal mediator with m-SNPs
            uint32_t curId; 
            int midx, eidx;
            double ztmpx, byztmp, seyztmp;

            bool plei_SNP, med_SNP, slct_SNP;

            for (int j=0; j<trans_Xm_idx_slct_tmp.size();j++){
  
                plei_SNP = false;
                med_SNP = false;
                slct_SNP = false;
                midx = trans_Xm_idx_slct_tmp[j];
                curId = (uint32_t)med_qtls_curId[midx];
                byztmp = gdata.byz[curId];
                seyztmp = gdata.seyz[curId];

                if (find(sigexpoId.begin(), sigexpoId.end(), curId) != sigexpoId.end()){ // sigexpoId already checked for reverse causality from outcome to exposure
                    int idx = (int)(find(e_Id.begin(), e_Id.end(), curId) - e_Id.begin());
                    if (idx != e_Id.size()){
                        for (int k=0; k<num_medPrbs; k++){   
                            zrevm = (fabs(bzx_totmr[idx])-fabs(trans_Xm_b(midx,k)))/sqrt(pow(trans_Xm_se(midx,k),2) + pow(sezx_totmr[idx],2));
                            if (zrevm < trev){
                                plei_SNP = true;
                                break;
                            }
                        }
                        if (!plei_SNP){
                            trans_Xm_idx_slct.push_back(midx);
                            trans_curId_slct.push_back(curId);
                            trans_Xm_chr.push_back(mdata._esi_chr[(int)med_qtls_included[midx]]);
                            trans_Xm_bp.push_back(mdata._esi_bp[(int)med_qtls_included[midx]]);
                            e_snp_Xm_idx_slct.push_back(trans_Xm_idx_slct.size()-1);
                            slct_SNP = true;
                            if (smrwk.medFocal_core & (trans_Xm_z(midx,0)>= z_m)){
                                m_snp_Xm_idx_slct.push_back(trans_Xm_idx_slct.size()-1);
                                Id4em_all.push_back(curId);
                                Id4em_m.push_back(curId);
                                zmz_foc_m.push_back(trans_Xm_z(midx,0));                              
                            }                          
                        } else { // remove pleiotropic SNP
                            rm_eID.push_back(idx);
                        }               
                    }
                } else {
                    if (smrwk.medFocal_core){
                        ztmpm=trans_Xm_z(midx,0);
                        zrevm = (fabs(trans_Xm_b(midx,0)) - fabs(byztmp))/sqrt(pow(trans_Xm_se(midx,0),2) + pow(seyztmp,2));
                        if ((ztmpm >= z_m) && (zrevm > trev)){
                            eidx = (int)(find(smrwk.curId.begin(), smrwk.curId.end(), curId) - smrwk.curId.begin());
                            if (eidx != smrwk.curId.size()){
                                trans_Xm_idx_slct.push_back(midx);
                                trans_curId_slct.push_back(curId);
                                trans_Xm_chr.push_back(mdata._esi_chr[(int)med_qtls_included[midx]]);
                                trans_Xm_bp.push_back(mdata._esi_bp[(int)med_qtls_included[midx]]);
                                m_snp_Xm_idx_slct.push_back(trans_Xm_idx_slct.size()-1);
                                slct_SNP = true;
                                Id4em_all.push_back(curId);
                                Id4em_m.push_back(curId);
                                zmz_foc_m.push_back(trans_Xm_z(midx,0));
                            }
                        }
                    }
                    if (!slct_SNP){
                        for (int k=0; k<num_medPrbs; k++){ 
                            ztmpm=trans_Xm_z(midx,k);
                            if (ztmpm >= z_m) med_SNP = true;
                            zrevm = (fabs(trans_Xm_b(midx,k)) - fabs(byztmp))/sqrt(pow(trans_Xm_se(midx,k),2) + pow(seyztmp,2));
                            if (zrevm < trev) plei_SNP = true;
                        }   
                        if (med_SNP && !plei_SNP){
                            trans_Xm_idx_slct.push_back(midx);
                            trans_curId_slct.push_back(curId);
                            trans_Xm_chr.push_back(mdata._esi_chr[(int)med_qtls_included[midx]]);
                            trans_Xm_bp.push_back(mdata._esi_bp[(int)med_qtls_included[midx]]);
                        }
                    }              
                }
            }    

            // define e_Ids that are not pleiotropic
            vector<uint32_t> e_Id_notpleio; 

            for (int i=0; i<e_Id.size(); i++){
                if (find(rm_eID.begin(), rm_eID.end(), i) == rm_eID.end()) e_Id_notpleio.push_back(e_Id[i]);
            }

            int num_med_qtls = trans_Xm_idx_slct.size();

            int num_focalSNPs = m_snp_Xm_idx_slct.size();;

            if (smrwk.medFocal_core && (m_snp_Xm_idx_slct.size() < 1)){
                smrwk.medFocal_core = false;
            }
           
            // Rank-based mediator pruning: for each mediator rank SNPs according to z-score
            MatrixXd trans_Xm_rank(trans_Xm_idx_slct.size(), med_probe_idx.size());
            vector<double> tmp_z(trans_Xm_idx_slct.size());
            for (int k=0; k<med_probe_idx.size(); k++){
                for (int j=0; j<trans_Xm_idx_slct.size(); j++) tmp_z[j] = trans_Xm_z(trans_Xm_idx_slct[j],k);
                int r = 0;
                vector<int> sorted_z = sort_indexes(tmp_z);
                for (int i=0; i<sorted_z.size();i++) {
                    trans_Xm_rank(sorted_z[i],k) = r;
                    r++;
                }
            } 

            // SNP rank is sum of ranks along the mediator axis
            VectorXd trans_snp_rank(trans_Xm_idx_slct.size());
            trans_snp_rank = trans_Xm_rank.rowwise().sum();

            int max_rank = trans_Xm_idx_slct.size()* med_probe_idx.size(); // assigned to E-SNPs 

            // Initialize design matrices, _em is for the focal mediation using E-SNPs & M-SNPs
            MatrixXd X_em, X, C_em, C, X_se, X_z;
            VectorXd yz_em, yz;
            VectorXd SEs_em, SEs;
            int num_SNPs, num_em_SNPs;
            vector<int> snp_indx_e, snp_indx_em;
            vector<string> rsvec, allele1vec, allele2vec;
            bool trans_med = false; // no trans-mediation performed if the number of SNPs is lower than the number of mediators
            mediation_status = "Not enough instrumental variables.";

            if (!smrwk.medFocal_core)
            {
                
                // attribute max_rank to E-SNPs
                for (int i=0; i<e_snp_Xm_idx_slct.size(); i++) trans_snp_rank[e_snp_Xm_idx_slct[i]] = max_rank;

                sub_indx.resize(0);
                vector<uint32_t> trans_snp_curId;

                vector<uint32_t> trans_chr_snp_curId;
                vector<double> trans_chr_snp_rank_all;
                vector<int> sub_chr_indx;
                //vector<int> sub_indx_tmp;
                vector<int> chr_indx;
                //vector<uint32_t> trans_snp_curId_chr;
                MatrixXd _Xchr;
                int num_SNPs_chr;

                // prune by chromosome
                
                for (int c=1; c<23; c++){
                    trans_chr_snp_curId.resize(0);
                    trans_chr_snp_rank_all.resize(0);
                    sub_chr_indx.resize(0);
                    chr_indx.resize(0);

                    for (int i=0; i<trans_snp_rank.size(); i++){
                        if (trans_Xm_chr[i] == c){
                            trans_chr_snp_curId.push_back((uint32_t)trans_curId_slct[i]); 
                            trans_chr_snp_rank_all.push_back(trans_snp_rank[i]);
                            chr_indx.push_back(i);
                        }
                    }
                    if (c == expochr){
                        // add the remaining E-SNPs
                        for (int i=0; i<e_Id.size(); i++){
                            if (find(trans_chr_snp_curId.begin(), trans_chr_snp_curId.end(), e_Id[i]) == trans_chr_snp_curId.end()){
                                trans_chr_snp_curId.push_back((uint32_t)e_Id[i]);
                                trans_chr_snp_rank_all.push_back(max_rank);
                                chr_indx.push_back(trans_snp_rank.size() + i);
                            }
                        }
                    }
                    num_SNPs_chr = trans_chr_snp_curId.size();
                    if (num_SNPs_chr > 0){
                        make_XMat(&bdata, trans_chr_snp_curId, _Xchr); //_X: one row one individual, one column one SNP
                        sbat_calcu_lambda(_Xchr, num_SNPs_chr, sbat_ld_cutoff, sub_chr_indx, trans_chr_snp_rank_all); // num_SNPs can change here
                        for (int j=0; j<sub_chr_indx.size(); j++){
                            sub_indx.push_back(chr_indx[sub_chr_indx[j]]);
                            trans_snp_curId.push_back((uint32_t)trans_chr_snp_curId[sub_chr_indx[j]]);  
                        }
                    }
                } 

                // add non-pleiotropic E-SNPs (if not already among trans_snp_curId)

                for(int j=0; j<e_Id_notpleio.size();j++)
                {
                    int idx=(int)(find(trans_snp_curId.begin(), trans_snp_curId.end(), e_Id_notpleio[j]) - trans_snp_curId.begin());
                    if (idx == trans_snp_curId.size()){
                        trans_snp_curId.push_back((uint32_t)e_Id_notpleio[j]);
                        sub_indx.push_back(trans_snp_rank.size());
                    }
                }

                // fill in effect sizes matrix for mediation regression
                num_SNPs = trans_snp_curId.size();

                if (num_SNPs > (num_medPrbs + 3)){
                    cout << num_SNPs << " independent SNPs are included in the trans-mediation calculation." << endl;
                    trans_med = true;

                    X.resize(num_SNPs, med_probe_idx.size() + 1);
                    X_se.resize(num_SNPs, med_probe_idx.size() + 1);
                    X_z.resize(num_SNPs, med_probe_idx.size() + 1);
                    yz.resize(num_SNPs);
                    SEs.resize((med_probe_idx.size() + 2)*num_SNPs);
                    if (ldmatrix){
                        make_XMat(&bdata, trans_snp_curId, _X);
                        cor_calc(C, _X);
                    } else C.setIdentity(num_SNPs, num_SNPs); 
                    rsvec.resize(num_SNPs);
                    allele1vec.resize(num_SNPs);
                    allele2vec.resize(num_SNPs);

                    for(int j=0; j<sub_indx.size();j++)
                    {                   
                        // trait data
                        yz(j) = gdata.byz[trans_snp_curId[j]];
                        SEs(num_SNPs*(med_probe_idx.size() + 1) + j) = gdata.seyz[trans_snp_curId[j]];

                        // exposure data            
                        int idx = (int)(find(smrwk.curId.begin(), smrwk.curId.end(), trans_snp_curId[j]) - smrwk.curId.begin());
                        if (idx != smrwk.curId.size()){
                            X(j, 0) = smrwk.bxz[idx];
                            X_se(j,0) = smrwk.sexz[idx];
                            X_z(j, 0) = fabs(X(j, 0)/X_se(j,0));
                            SEs(j) = smrwk.sexz[idx];
                            if (find(e_Id.begin(), e_Id.end(), trans_snp_curId[j]) != e_Id.end()){
                                snp_indx_e.push_back(j);
                                snp_indx_em.push_back(j);
                            }
                        } else {
                            X(j,0) = 0;
                            X_se(j,0) = default_expo_se;
                            X_z(j, 0) = fabs(X(j, 0)/X_se(j,0));
                            SEs(j) = default_expo_se;
                        }

                        // mediator data

                        if (sub_indx[j] < trans_snp_rank.size()){
                            // SNP Info
                            rsvec[j] = trans_rsvec[trans_Xm_idx_slct[sub_indx[j]]];
                            allele1vec[j] = trans_allele1vec[trans_Xm_idx_slct[sub_indx[j]]];
                            allele2vec[j] = trans_allele2vec[trans_Xm_idx_slct[sub_indx[j]]];

                            for(int k=0; k<med_probe_idx.size(); k++){                         
                                X(j,k+1) = trans_Xm_b(trans_Xm_idx_slct[sub_indx[j]], k);
                                X_se(j,k+1) = trans_Xm_se(trans_Xm_idx_slct[sub_indx[j]], k);
                                X_z(j,k+1) = trans_Xm_z(trans_Xm_idx_slct[sub_indx[j]], k);
                                SEs(num_SNPs*(k+1) + j) = trans_Xm_se(trans_Xm_idx_slct[sub_indx[j]], k);
                            }
                        } else {
                            int idx = (int)(find(e_Id.begin(), e_Id.end(), trans_snp_curId[j]) - e_Id.begin());
                            // SNP Info
                            rsvec[j] = e_rsvec[idx];
                            allele1vec[j] = e_allele1vec[idx];
                            allele2vec[j] = e_allele2vec[idx];

                            for(int k=0; k<med_probe_idx.size(); k++){
                                X(j,k+1) = e_Xm_b(idx, med_indx[k]);
                                X_se(j,k+1) = e_Xm_se(idx, med_indx[k]);
                                X_z(j,k+1) = fabs(X(j,k+1)/X_se(j,k+1));
                                SEs(num_SNPs*(k+1) + j) = e_Xm_se(idx, med_indx[k]);
                            }
                        }
                    }
                } else {
                    cout << num_SNPs << " independent SNPs were found which is less than the number of mediators." << endl;
                    cout << "No trans-mediation performed." << endl;
                }
            } else {

                // trans-SNPs will be outside a window of 1000kB around m-SNPs and e-SNPs
                
                // get window boundaries
                VectorXd em_snps_bp(num_e_SNPs + num_focalSNPs);
                for (int i=0; i< num_e_SNPs; i++) em_snps_bp(i) = e_bpsnp[i];
                for (int i=0; i< num_focalSNPs; i++) em_snps_bp(num_e_SNPs+i) = trans_Xm_bp[m_snp_Xm_idx_slct[i]];

                int lower_bp = em_snps_bp.minCoeff() - 1000000;
                int upper_bp = em_snps_bp.maxCoeff() + 1000000;

                sub_indx.resize(0);
                vector<uint32_t> trans_snp_curId;
                vector<uint32_t> trans_snp_curId_all;
                vector<uint32_t> trans_chr_snp_curId;
                vector<double> trans_chr_snp_rank_all;
                vector<int> sub_chr_indx;
                vector<int> chr_indx;
                MatrixXd _Xchr;
                int num_SNPs_chr;

                // prune by chromosome

                for (int c=1; c<23; c++){

                    trans_chr_snp_curId.resize(0);
                    trans_chr_snp_rank_all.resize(0);
                    sub_chr_indx.resize(0);
                    chr_indx.resize(0);

                    for (int i=0; i<trans_snp_rank.size(); i++){
                        if (c == expochr){
                            if ((trans_Xm_chr[i] == c) && ((trans_Xm_bp[i] < lower_bp) || (trans_Xm_bp[i] > upper_bp))){
                                trans_chr_snp_curId.push_back((uint32_t)trans_curId_slct[i]); 
                                trans_chr_snp_rank_all.push_back(trans_snp_rank[i]);
                                chr_indx.push_back(i);
                            }
                        } else if (trans_Xm_chr[i] == c){
                            trans_chr_snp_curId.push_back((uint32_t)trans_curId_slct[i]); 
                            trans_chr_snp_rank_all.push_back(trans_snp_rank[i]);
                            chr_indx.push_back(i);
                        }
                    }
                    
                    num_SNPs_chr = trans_chr_snp_curId.size();
                    if (num_SNPs_chr > 0){
                        make_XMat(&bdata, trans_chr_snp_curId, _Xchr); //_X: one row one individual, one column one SNP
                        sbat_calcu_lambda(_Xchr, num_SNPs_chr, sbat_ld_cutoff, sub_chr_indx, trans_chr_snp_rank_all); // num_SNPs can change here
                        for (int j=0; j<sub_chr_indx.size(); j++){
                            sub_indx.push_back(chr_indx[sub_chr_indx[j]]);
                            trans_snp_curId.push_back((uint32_t)trans_chr_snp_curId[sub_chr_indx[j]]); 
                            trans_snp_curId_all.push_back((uint32_t)trans_chr_snp_curId[sub_chr_indx[j]]);
                        }
                    }
                }

                // Add M-SNPs to E-SNPs and prune 

                vector<uint32_t> em_snp_curId;
                vector<double> em_snp_rank;
                bool expo_snp, med_snp;

                for(int j=0; j<Id4em_all.size();j++)
                {
                    med_snp = false;
                    expo_snp = find(Id4em_x.begin(), Id4em_x.end(), Id4em_all[j]) != Id4em_x.end();
                    
                    int midx=(int)(find(Id4em_m.begin(), Id4em_m.end(), Id4em_all[j]) - Id4em_m.begin());
                    if(midx<Id4em_m.size()) med_snp = true;

                    if (expo_snp || med_snp){               
                        em_snp_curId.push_back((uint32_t)Id4em_all[j]);
                        if (expo_snp){
                            em_snp_rank.push_back(100);
                        }
                        else {
                            em_snp_rank.push_back(zmz_foc_m[midx]);
                        }
                    }
                }

                num_em_SNPs = em_snp_curId.size();        
                make_XMat(&bdata, em_snp_curId, _X);
                vector<int> em_sub_indx;
                sbat_calcu_lambda(_X, num_em_SNPs, sbat_ld_cutoff, em_sub_indx, em_snp_rank);
                
                vector<uint32_t> em_snp_curId_slct;
                for (int j = 0; j < em_sub_indx.size(); j++){
                    em_snp_curId_slct.push_back((uint32_t)em_snp_curId[em_sub_indx[j]]);
                    trans_snp_curId_all.push_back((uint32_t)em_snp_curId[em_sub_indx[j]]);
                } 
                num_em_SNPs = em_sub_indx.size();

                if (num_em_SNPs > 2){
                    cout << num_em_SNPs << " independent SNPs are included in the focal mediator MR calculation." << endl;
                } else {
                    cout << num_em_SNPs << " independent SNPs were found which is less than 3." << endl;
                    cout << "No focal mediator MR calculation performed." << endl;
                    smrwk.medFocal = false;
                }

                num_SNPs = em_sub_indx.size() + sub_indx.size();

                // fill in effect sizes matrix for mediation regression

                if (num_SNPs > (med_probe_idx.size() + 3)){
                    cout << num_SNPs << " independent SNPs are included in the trans-mediation calculation." << endl;
                    trans_med = true;
                } else {
                    cout << num_SNPs << " independent SNPs were found which is less than the number of mediators." << endl;
                    cout << "No trans-mediation performed." << endl;
                }

                // Get C matrix
                if (ldmatrix){
                    make_XMat(&bdata, trans_snp_curId_all, _X);
                    cor_calc(C, _X);
                    make_XMat(&bdata, em_snp_curId_slct, _X);
                    cor_calc(C_em, _X);
                } else {
                    C.setIdentity(num_SNPs, num_SNPs); 
                    C_em.setIdentity(num_em_SNPs, num_em_SNPs);
                }
                
                X_em.resize(num_em_SNPs, 2);
                yz_em.resize(num_em_SNPs);
                SEs_em.resize(3*num_em_SNPs);

                X.resize(num_SNPs, med_probe_idx.size() + 1);
                X_se.resize(num_SNPs, med_probe_idx.size() + 1);
                X_z.resize(num_SNPs, med_probe_idx.size() + 1);
                yz.resize(num_SNPs);
                SEs.resize((med_probe_idx.size() + 2)*num_SNPs);

                rsvec.resize(num_SNPs);
                allele1vec.resize(num_SNPs);
                allele2vec.resize(num_SNPs);

                // fill in EM data
                for(int j=0; j<num_em_SNPs;j++)
                {
                    // trait data
                    yz_em(j) = gdata.byz[(int)em_snp_curId_slct[j]];
                    SEs_em(num_em_SNPs*2 + j) = gdata.seyz[(int)em_snp_curId_slct[j]];

                    yz(j) = gdata.byz[(int)em_snp_curId_slct[j]];
                    SEs(num_SNPs*(med_probe_idx.size() + 1) + j) = gdata.seyz[(int)em_snp_curId_slct[j]];

                    // exposure data       
                    int idx = (int)(find(smrwk.curId.begin(), smrwk.curId.end(), em_snp_curId_slct[j]) - smrwk.curId.begin());   
                    X_em(j, 0) = smrwk.bxz[idx];
                    SEs_em(j) = smrwk.sexz[idx];

                    X(j, 0) = smrwk.bxz[idx];
                    X_se(j,0) = smrwk.sexz[idx];
                    X_z(j,0) = fabs(X(j, 0)/X_se(j,0));
                    SEs(j) = smrwk.sexz[idx];

                    if (find(e_Id.begin(), e_Id.end(), em_snp_curId_slct[j]) != e_Id.end()){
                        snp_indx_e.push_back(j);
                    }
                    snp_indx_em.push_back(j);

                    // mediator data
                    idx = (int)(find(e_Id.begin(), e_Id.end(), em_snp_curId_slct[j]) - e_Id.begin());
                    if (idx < e_Id.size()){
                        // SNP Info
                        rsvec[j] = e_rsvec[idx];
                        allele1vec[j] = e_allele1vec[idx];
                        allele2vec[j] = e_allele2vec[idx];

                        X_em(j,1) = e_Xm_b(idx, 0);
                        SEs_em(num_em_SNPs + j) = e_Xm_se(idx, 0);
                        for(int k=0; k<med_probe_idx.size(); k++){
                            X(j,k+1) = e_Xm_b(idx, med_indx[k]);
                            X_se(j,k+1) = e_Xm_se(idx, med_indx[k]);
                            X_z(j,k+1) = fabs(X(j,k+1)/X_se(j,k+1));
                            SEs(num_SNPs*(k+1) + j) = e_Xm_se(idx, med_indx[k]);
                        }
                    } else {
                        idx = (int)(find(trans_curId_slct.begin(), trans_curId_slct.end(), em_snp_curId_slct[j]) - trans_curId_slct.begin());
                        // SNP Info
                        rsvec[j] = trans_rsvec[trans_Xm_idx_slct[idx]];
                        allele1vec[j] = trans_allele1vec[trans_Xm_idx_slct[idx]];
                        allele2vec[j] = trans_allele2vec[trans_Xm_idx_slct[idx]];
                        
                        X_em(j,1) = trans_Xm_b(trans_Xm_idx_slct[idx], 0);
                        SEs_em(num_em_SNPs + j) = trans_Xm_se(trans_Xm_idx_slct[idx], 0);
                        for(int k=0; k<med_probe_idx.size(); k++){
                            X(j,k+1) = trans_Xm_b(trans_Xm_idx_slct[idx], k);
                            X_se(j,k+1) = trans_Xm_se(trans_Xm_idx_slct[idx], k);
                            X_z(j,k+1) = trans_Xm_z(trans_Xm_idx_slct[idx], k);
                            SEs(num_SNPs*(k+1) + j) = trans_Xm_se(trans_Xm_idx_slct[idx], k);
                        }
                    }
                }

                // fill in trans data

                for(int j=0; j<sub_indx.size();j++)
                {
                    // trait data
                    yz(num_em_SNPs+j) = gdata.byz[trans_snp_curId[j]];
                    SEs(num_SNPs*(med_probe_idx.size() + 1) + num_em_SNPs+j) = gdata.seyz[trans_snp_curId[j]];

                    // exposure data         
                    int idx = (int)(find(smrwk.curId.begin(), smrwk.curId.end(), trans_snp_curId[j])- smrwk.curId.begin());
                    if (idx < smrwk.curId.size()){
                        X(num_em_SNPs+j, 0) = smrwk.bxz[idx];
                        X_se(num_em_SNPs+j, 0) = smrwk.sexz[idx];
                        SEs(num_em_SNPs+j) = smrwk.sexz[idx];
                    } else {
                        X(num_em_SNPs+j,0) = 0;
                        X_se(num_em_SNPs+j,0) = default_expo_se;
                        SEs(num_em_SNPs+j) = default_expo_se;
                    }
                    X_z(num_em_SNPs+j,0) = fabs(X(num_em_SNPs+j, 0)/X_se(num_em_SNPs+j,0));

                    // mediator data
                    for(int k=0; k<med_probe_idx.size(); k++){
                        X(num_em_SNPs+j,k+1) = trans_Xm_b(trans_Xm_idx_slct[sub_indx[j]], k);
                        X_se(num_em_SNPs+j,k+1) = trans_Xm_se(trans_Xm_idx_slct[sub_indx[j]], k);
                        X_z(num_em_SNPs+j,k+1) = trans_Xm_z(trans_Xm_idx_slct[sub_indx[j]], k);
                        SEs(num_SNPs*(k+1) + num_em_SNPs+j) = trans_Xm_se(trans_Xm_idx_slct[sub_indx[j]], k);
                    }
                    // SNP Info
                    rsvec[num_em_SNPs+j] = trans_rsvec[trans_Xm_idx_slct[sub_indx[j]]];
                    allele1vec[num_em_SNPs+j] = trans_allele1vec[trans_Xm_idx_slct[sub_indx[j]]];
                    allele2vec[num_em_SNPs+j] = trans_allele2vec[trans_Xm_idx_slct[sub_indx[j]]];
                }
            }

            // SNP effects and phenotypic correlation matrix

            if (get_snp_effects_flg && trans_med)
            {
                cout << "Saving design matrix ... " << endl;
                FILE* snpfile;
                string snpfilename = string(outFileName) + "." + exponame + ".snps";

                snpfile = fopen(snpfilename.c_str(), "w");
                if (!(snpfile)) {
                    printf("ERROR: open error %s\n", snpfilename.c_str());
                    exit(1);
                }
                string snpstr="SNP\tA1\tA2\tbeta.out\tse.out\tbeta." + exponame + "\tse." + exponame;

                for(int k=0;k<num_medPrbs;k++){
                    snpstr += "\tbeta." + med_prbs[k] + "\tse." + med_prbs[k];
                }

                snpstr += "\n";
                
                if(fputs_checked(snpstr.c_str(),snpfile))
                {
                    printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                    exit(EXIT_FAILURE);
                }

                for (int j=0; j<num_SNPs; j++){
                    snpstr = rsvec[j] + '\t' + allele1vec[j] + '\t' + allele2vec[j] + '\t';
                    snpstr += atos(yz[j]) + '\t' + atos(SEs(num_SNPs*(med_probe_idx.size() + 1)+j));

                    for (int k=0; k<num_medPrbs+1; k++){
                        snpstr += '\t' + atos(X(j,k)) + '\t' + atos(X_se(j,k));
                    }
            
                    snpstr += "\n"; 

                    if(fputs_checked(snpstr.c_str(),snpfile))
                    {
                        printf("ERROR: in writing file %s .\n", snpfilename.c_str());
                        exit(EXIT_FAILURE);
                    }
                }

                fclose(snpfile);              

            }

            // Effect shrinkage: set non-significant effects to zero except for E-SNPs; keep those in a separate matrix and re-add them to the matrix
            int num_E = snp_indx_e.size();

            if (num_E < 3){
                cout << "Less than 3 SNPs had a higher exposure than mediator effect. No trans-mediation analysis performed." << endl;
                mediation_status = "Less than 3 E-SNPs.";
                trans_med = false;  
            } 

            MatrixXd X_E(num_E, med_probe_idx.size() + 1);
            VectorXd y_E(num_E);
            for (int j=0; j< num_E; j++){
                for (int k=0;k<(med_probe_idx.size()+1);k++){
                    X_E(j,k) = X(snp_indx_e[j],k);
                    y_E(j) = yz(snp_indx_e[j]);
                }
            }

            double z_shrinkage = fabs(qnorm5(p_shrinkage/2)); 
            X = (X_z.array() < z_shrinkage).select(0,X);

            for (int j=0; j< num_E; j++){
                for (int k=0;k<(med_probe_idx.size()+1);k++){
                    X(snp_indx_e[j],k) = X_E(j,k);
                }
            }
            

            // step4: calculate MR effects

            cout << "Calculating MR effects ..." << endl;

            // 1. total exposure effect

            double b_tot_ivw, se_tot_ivw, p_tot_ivw;
            if (ldmatrix) mr_ivw_LD(bzx_totmr, bzy_totmr, sezx_totmr, sezy_totmr, Ctot, b_tot_ivw, se_tot_ivw, p_tot_ivw, N_ref);
            else mr_ivw(bzx_totmr, bzy_totmr, sezx_totmr, sezy_totmr, b_tot_ivw, se_tot_ivw, p_tot_ivw);
            cout << "Total MR effects: b = " << b_tot_ivw << "; p = " << p_tot_ivw << endl;

            // 2. Direct effect through focal mediator

            double b_direct_med = -9, se_direct_med = -9, p_direct_med = -9;
            int n_snps_med;

            if (smrwk.medFocal){

                if (smrwk.medFocal_core){
                    n_snps_med = num_em_SNPs;
                } else {
                    n_snps_med = num_e_SNPs;
                    X_em.resize(num_e_SNPs, 2);
                    yz_em.resize(num_e_SNPs);
                    SEs_em.resize(3*num_e_SNPs);
                    C_em = Ctot;

                    for(int j=0; j<num_e_SNPs;j++)
                    {
                        yz_em(j) = bzy_totmr[j];
                        X_em(j,0) = bzx_totmr[j];
                        X_em(j,1) = e_Xm_b(j, 0);

                        SEs_em(j) = sezx_totmr[j];
                        SEs_em(num_e_SNPs + j) = e_Xm_se(j, 0);
                        SEs_em(num_e_SNPs*2 + j) = sezy_totmr[j];
                    }
                }

                VectorXd betas, vars;
                 
                mvmr_ivw_LD(X_em, yz_em, SEs_em, C_em, betas, vars, N_ref);

                b_direct_med = betas[0];
                se_direct_med = sqrt(vars[0]);
                p_direct_med = 2*pnorm(fabs(b_direct_med/se_direct_med));

                cout << "Focal mediator MR effects: b = " << b_direct_med << "; p = " << p_direct_med << endl;

            }

            // Direct effect through trans mediators

            double b_direct_trans = -9, se_direct_trans = -9, p_direct_trans = -9;

            if (trans_med){
                mediation_status = "Mediation analysis conducted.";
                VectorXd betas, vars;
                mvmr_ivw_LD(X, yz, SEs, C, betas, vars, N_ref);
                b_direct_trans = betas(0);
                se_direct_trans = sqrt(vars[0]);
                p_direct_trans = 2*pnorm(fabs(b_direct_trans/se_direct_trans));
                cout << "Trans-mediation MR effects: b = " << b_direct_trans << "; p = " << p_direct_trans << endl;

                for (int k=0; k<num_medPrbs; k++){
                    med_info.b_my.push_back(betas(k+1));
                    med_info.se_my.push_back(sqrt(vars[k+1]));
                    med_info.p_my.push_back(2*pnorm(fabs(betas(k+1)/sqrt(vars[k+1]))));
                    med_info.n_my.push_back(num_SNPs);
                }

            }

            expostr = exponame + '\t' + expogene + '\t' + mediation_status + '\n';
            if(fputs_checked(expostr.c_str(),expolst))
            {
                printf("ERROR: in writing file %s .\n", expolstfile.c_str());
                exit(EXIT_FAILURE);
            }

            // output snp set list

            string setstr=exponame+'\n';
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<num_e_SNPs;j++)
            {
                setstr=e_snpnames[j]+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            setstr="end\n";
            if(fputs_checked(setstr.c_str(),setlst))
            {
                printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                exit(EXIT_FAILURE);
            }
            // end of output

            //output mediator information

            if (trans_med){
                for(int k=0;k<num_medPrbs;k++)
                {
                    medstr = exponame + '\t' + expogene + '\t' + atos(expochr) + '\t' + atos(expobp) + '\t';
                    medstr += med_info.Med_ID[k] + '\t' + med_info.Gene_ID[k] + '\t' + atos(med_info.Med_Chr[k]) + '\t' + atos(med_info.Med_bp[k]) + '\t'; 
                    medstr += atos(med_info.b_xm[k]) + '\t' + atos(med_info.se_xm[k]) + '\t' + atos(med_info.p_xm[k]) + '\t' + atos(med_info.n_xm[k]) + '\t';
                    medstr += atos(med_info.b_my[k]) + '\t' + atos(med_info.se_my[k]) + '\t' + atos(med_info.p_my[k]) + '\t' + atos(med_info.n_my[k]) + '\n';
                    if(fputs_checked(medstr.c_str(),medlst))
                    {
                        printf("ERROR: in writing file %s .\n", medlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
            }

            // end of output 

            outstr =  exponame + '\t' + expogene + '\t' + (smrwk.medFocal ? medname : "NA") + '\t' + (smrwk.medFocal ? medgene : "NA") + '\t';
            outstr += atos(expochr) + '\t' + atos(expobp) + '\t'; 
            outstr += atos(b_tot_ivw)  + '\t' + atos(se_tot_ivw) + '\t' + atos(p_tot_ivw) + '\t';
            outstr += (smrwk.medFocal ? atos(b_direct_med) : "NA")  + '\t' + (smrwk.medFocal ? atos(se_direct_med) : "NA") + '\t' + (smrwk.medFocal ? atos(p_direct_med) : "NA") + '\t';
            outstr += (trans_med ? atos(b_direct_trans) : "NA") + '\t' + (trans_med ? atos(se_direct_trans) : "NA") + '\t' + (trans_med ? atos(p_direct_trans) : "NA") + '\t';
            outstr += atos(num_e_SNPs) + '\t' + (jj > 0 ? atos(n_snps_med) : "NA") + '\t' + atos(num_SNPs) + '\t' + atos(num_E) + '\t'+ atos(num_medPrbs) + '\n';
            
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
            write_count++;
            
        }

        cout<<"\nMR mediation analysis on cis-exposure and trans-mediators completed.\n Results of "<<write_count<<" probes have been saved in the file [" + smrfile + "]."<<endl;     

        fclose(smr);
        fclose(medlst); 
        fclose(setlst);   
        fclose(fptr); 
        fclose(expolst);      
    }

    void calcul_cormat(char* outFileName, char* eqtlFileName, char* snplstName, char* problstName, char* genelistName, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl)
    {
        string logstr;     
        
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading QTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
        
        uint64_t epinum = eqtlinfo._probNum;
        uint64_t esinum = eqtlinfo._snpNum;

        vector<uint64_t> e_snp_idx;
        vector<uint64_t> e_prb_idx;

        for (int i = 0; i < eqtlinfo._esi_include.size(); i++) e_snp_idx.push_back(eqtlinfo._snp_name_map[eqtlinfo._esi_rs[eqtlinfo._esi_include[i]]]);
        for (int i = 0; i < eqtlinfo._include.size(); i++) e_prb_idx.push_back(eqtlinfo._probe_name_map[eqtlinfo._epi_prbID[eqtlinfo._include[i]]]);

        // get the effect sizes
        
        FILE * fptr;
        string besd_file = string(eqtlFileName)+".besd";
        fptr = fopen(besd_file.c_str(),"rb");
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        int infoLen=sizeof(uint32_t);
        if(indicator==DENSE_FILE_TYPE_3) infoLen=RESERVEDUNITS*sizeof(int);

        if(indicator!=DENSE_FILE_TYPE_1 && indicator!=DENSE_FILE_TYPE_3) {
            cout << "ERROR: the BESD file is not in dense BESD format." << endl;
            exit(EXIT_FAILURE);
        }
        
        MatrixXd C, X_mat(eqtlinfo._esi_include.size(), eqtlinfo._include.size());
        float* tmp=(float*)malloc(sizeof(float)*esinum<<1);
        cout<<"Getting SNP effect sizes..."<<endl;

        for (int i = 0; i < eqtlinfo._include.size(); i++){

            fseek(fptr,((e_prb_idx[i]*esinum)<<3)+infoLen, SEEK_SET);
            fread(tmp, sizeof(float), esinum*2,fptr);

            for (int j = 0; j< eqtlinfo._esi_include.size(); j++){

                double se=tmp[e_snp_idx[j] + esinum];
                if(fabs(se+9)<1e-6) X_mat(j,i) = 0;
                else X_mat(j,i) = tmp[e_snp_idx[j]];
            }
        }

        fclose(fptr);
        cout<<"Calculating correlation matrix..."<<endl;
        cor_calc(C, X_mat); 

        cout<<"Writing to output file..."<<endl;
        FILE* corfile;
        string corfilename = string(outFileName) + ".cor";
        
        corfile = fopen(corfilename.c_str(), "w");
        if (!(corfile)) {
            printf("ERROR: open error %s\n", corfilename.c_str());
            exit(1);
        }

        string outstr="id";

        for(int i=0;i<eqtlinfo._include.size();i++){
            outstr += "\t" + eqtlinfo._epi_prbID[eqtlinfo._include[i]];
        }
        outstr += "\n";
        
        if(fputs_checked(outstr.c_str(),corfile))
        {
            printf("ERROR: in writing file %s .\n", corfilename.c_str());
            exit(EXIT_FAILURE);
        }

        for (int i=0; i<eqtlinfo._include.size(); i++){
            outstr = eqtlinfo._epi_prbID[eqtlinfo._include[i]];
            for (int j=0; j<eqtlinfo._include.size(); j++){
                outstr += '\t' + atos(C(i,j));
            }   
            outstr += "\n"; 

            if(fputs_checked(outstr.c_str(),corfile))
            {
                printf("ERROR: in writing file %s .\n", corfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }

        fclose(corfile);              

    }

    void calcul_cormat(char* outFileName, char* eqtlFileName, char* eqtlFileName2, char* snplstName, char* problstName, char* genelistName, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl)
    {
        eqtlInfo eqtlinfo;
        cout<<"Reading QTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
            read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._valNum!=0){
                cout << "ERROR: the BESD file is not in dense BESD format." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
     
        eqtlInfo eqtlinfo2;
        cout<<"Reading second QTL summary data..."<<endl;

        read_esifile(&eqtlinfo2, string(eqtlFileName2)+".esi");
        allele_check(&eqtlinfo, &eqtlinfo2);
        read_epifile(&eqtlinfo2, string(eqtlFileName2)+".epi");

        cout << "Getting SNP effect sizes of second QTL summary data ..." << endl;

        uint64_t epinum = eqtlinfo2._probNum;
        uint64_t esinum = eqtlinfo2._snpNum;

        vector<uint64_t> e_snp_idx;
        vector<uint64_t> e_prb_idx;

        for (int i = 0; i < eqtlinfo2._esi_include.size(); i++) e_snp_idx.push_back(eqtlinfo2._snp_name_map[eqtlinfo2._esi_rs[eqtlinfo2._esi_include[i]]]);
        for (int i = 0; i < eqtlinfo2._include.size(); i++) e_prb_idx.push_back(eqtlinfo2._probe_name_map[eqtlinfo2._epi_prbID[eqtlinfo2._include[i]]]);

        FILE * fptr;
        string besd_file = string(eqtlFileName2)+".besd";
        fptr = fopen(besd_file.c_str(),"rb");
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        int infoLen=sizeof(uint32_t);
        if(indicator==DENSE_FILE_TYPE_3) infoLen=RESERVEDUNITS*sizeof(int);

        if(indicator!=DENSE_FILE_TYPE_1 && indicator!=DENSE_FILE_TYPE_3) {
            cout << "ERROR: the BESD file is not in dense BESD format." << endl;
            exit(EXIT_FAILURE);
        }

        MatrixXd X_mat(eqtlinfo2._esi_include.size(), eqtlinfo2._include.size());
        float* tmp=(float*)malloc(sizeof(float)*esinum<<1);

        for (int i = 0; i < eqtlinfo2._include.size(); i++){

            fseek(fptr,((e_prb_idx[i]*esinum)<<3)+infoLen, SEEK_SET);
            fread(tmp, sizeof(float), esinum*2,fptr);

            for (int j = 0; j< eqtlinfo2._esi_include.size(); j++){

                double se=tmp[e_snp_idx[j] + esinum];
                if(fabs(se+9)<1e-6) X_mat(j,i) = 0;
                else X_mat(j,i) = tmp[e_snp_idx[j]];
            }
        }
        fclose(fptr);

        cout << "Calculating pairwise correlations ..." << endl;

        MatrixXd C(eqtlinfo2._include.size(), eqtlinfo._include.size());
        vector <double> x,y;

        for (int i = 0; i < eqtlinfo2._include.size(); i++){
            for (int j = 0; j < eqtlinfo._include.size(); j++){
                x.clear(), y.clear();
                for (int k=0; k<eqtlinfo2._esi_include.size(); k++){
                    x.push_back(X_mat(k,i));
                    double se=eqtlinfo._sexz[j][k];
                    if(fabs(se+9)<1e-6) y.push_back(0);
                    else y.push_back(eqtlinfo._bxz[j][k]);
                }
                C(i,j) = cor(x,y);
            }
        }

        cout<<"Writing to output file..."<<endl;
        FILE* corfile;
        string corfilename = string(outFileName) + ".cor";
        
        corfile = fopen(corfilename.c_str(), "w");
        if (!(corfile)) {
            printf("ERROR: open error %s\n", corfilename.c_str());
            exit(1);
        }

        string outstr="id";

        for(int i=0;i<eqtlinfo._include.size();i++){
            outstr += "\t" + eqtlinfo._epi_prbID[eqtlinfo._include[i]];
        }
        outstr += "\n";
        
        if(fputs_checked(outstr.c_str(),corfile))
        {
            printf("ERROR: in writing file %s .\n", corfilename.c_str());
            exit(EXIT_FAILURE);
        }

        for (int i=0; i<eqtlinfo2._include.size(); i++){
            outstr = eqtlinfo2._epi_prbID[eqtlinfo2._include[i]];
            for (int j=0; j<eqtlinfo._include.size(); j++){
                outstr += '\t' + atos(C(i,j));
            }   
            outstr += "\n"; 

            if(fputs_checked(outstr.c_str(),corfile))
            {
                printf("ERROR: in writing file %s .\n", corfilename.c_str());
                exit(EXIT_FAILURE);
            }
        }

        fclose(corfile);       

    }

}