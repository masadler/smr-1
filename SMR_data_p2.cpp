//
//  SMR_data_p2.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 10/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "SMR_data_p2.h"
FILE* techeQTLfile=NULL;
namespace SMRDATA
{
    double adjSE(double beta,double p)
    {
        if(p==1)
        {
            return 1e10;
        }
        else
        {
            double z2=qchisq(p, 1);
            return fabs(beta/sqrt(z2));
        }
    }
    // below for txt 2 besd
    int read_probeinfolst(vector<probeinfolst> &prbiflst,char* syllabusName)
    {
        ifstream flptr(syllabusName);
        if (!flptr) {
            printf("ERROR: can not open the file %s to read.\n",syllabusName);
            exit(EXIT_FAILURE);
        }
        printf( "Reading eQTL probe information from \"%s\".\n", syllabusName);
        map<string, int> probe_map;
        map<string, int> probe_cp_map;
        map<string, int> probe_cpb_map;
        long mapsize=0;
        char buf[MAX_LINE_SIZE];
        flptr.getline(buf, MAX_LINE_SIZE); //header
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",syllabusName);
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: %s should have headers that start with \"Chr\".\n", syllabusName);
            exit(EXIT_FAILURE);
        }
        int lineNum(0);
        while(!flptr.eof())
        {
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0') lineNum++;
        }
        flptr.clear(ios::goodbit);
        flptr.seekg (0, ios::beg);
        flptr.getline(buf, MAX_LINE_SIZE); //header
        int fcount=0;
        while(!flptr.eof())
        {
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=7 && col_num!=8) {
                    printf("ERROR: column number is not correct in row %d with probe name %s.\n", fcount+2, vs_buf[1].c_str());
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na")
                {
                    printf("ERROR: probe name is missing in row %d.\n", fcount+2);
                    exit(EXIT_FAILURE);
                }
                to_upper(vs_buf[0]);
                to_upper(vs_buf[3]);
                if(vs_buf[0]=="NA" || vs_buf[3]=="NA")
                {
                    printf("ERROR: probe chromosome information or probe position information of probe %s is missing in row %d.\n", vs_buf[1].c_str(),fcount+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[6]=="NA" || vs_buf[6]=="na")
                {
                    printf("ERROR: ESD file path of probe %s is missing in row %d.\n", vs_buf[1].c_str(),fcount+2);
                    exit(EXIT_FAILURE);
                }
                string cpstr=vs_buf[0]+":"+vs_buf[1];
                string cpbstr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[3];
                probe_cp_map.insert(pair<string,int>(cpstr,fcount));
                probe_cpb_map.insert(pair<string,int>(cpbstr,fcount));
                probe_map.insert(pair<string,int>(vs_buf[6],fcount));
                if(probe_cp_map.size()!=probe_cpb_map.size())
                {
                    printf("ERROR: probe %s has multiple positions.\n",vs_buf[1].c_str() );
                    exit(EXIT_FAILURE);
                }
                if(mapsize==probe_map.size())
                {
                    printf("ERROR: Duplicated ESD file name found for %s.\n",vs_buf[6].c_str() );
                    exit(EXIT_FAILURE);
                } else {
                    mapsize=probe_map.size();
                }
                probeinfolst tmpinfo;
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                tmpinfo.probechr=tmpchr;
                strcpy2(&tmpinfo.probeId, vs_buf[1]);
                tmpinfo.gd=atoi(vs_buf[2].c_str());
                tmpinfo.bp=atoi(vs_buf[3].c_str());
                strcpy2(&tmpinfo.genename, vs_buf[4]);
                tmpinfo.orien=vs_buf[5].c_str()[0];
                strcpy2(&tmpinfo.esdpath, vs_buf[6]);
                if(col_num==8){
                    strcpy2(&tmpinfo.bfilepath, vs_buf[7]);
                } else tmpinfo.bfilepath=NULL;
                prbiflst.push_back(tmpinfo);
                fcount++;
            }
        }
        flptr.close();
        printf("%d ESD files of total %ld probes to be included from \"%s\".\n",lineNum, probe_cp_map.size(), syllabusName);
        return lineNum;
    }
    void read_smr_sa(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<string> &a1, vector<string> &a2, vector<float> &freq,vector<float> &beta,vector<float> &se,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Error: can not open the file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL summary data from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"Chr\"");
            exit(EXIT_FAILURE);
        }
        map<string, int> currs_map;
        long currssize=0;
        while(!flptr.eof() && !gzeof(gzfile))
        {
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=9) {
                    printf("ERROR: column number is not correct in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                string tmpstr=vs_buf[0]+":"+vs_buf[1];
                currs_map.insert(pair<string,int>(tmpstr,lineNum));
                if(currssize<currs_map.size())
                {
                    if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                        printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    int tmpchr;
                    if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                    else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                    else tmpchr=atoi(vs_buf[0].c_str());
                    chr.push_back(tmpchr);
                    rs.push_back(vs_buf[1]);
                    if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                        printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }

                    bp.push_back(atoi(vs_buf[2].c_str()));
                    to_upper(vs_buf[3]);
                    a1.push_back(vs_buf[3]);
                    to_upper(vs_buf[4]);
                    a2.push_back(vs_buf[4]);
                    if(vs_buf[5].compare("NA") && vs_buf[5].compare("na"))
                    {
                        float fq=atof(vs_buf[5].c_str());
                        if(fq<0 || fq>1)
                        {
                            printf("ERROR: allele frequency should be between 0 and 1 in row %d.\n",lineNum+2);
                            exit(EXIT_FAILURE);
                        }
                        freq.push_back(fq);
                    } else {
                        freq.push_back(-9);
                    }
                    to_upper(vs_buf[6]);
                    to_upper(vs_buf[7]);
                    to_upper(vs_buf[8]);
                    if(vs_buf[6].compare("NA") && (vs_buf[7].compare("NA") || vs_buf[8].compare("NA"))){
                        beta.push_back(atof(vs_buf[6].c_str()));
                        if(vs_buf[8]=="NA") se.push_back(atof(vs_buf[7].c_str()));
                        else {
                            double betatmp=atof(vs_buf[6].c_str());
                            double ptmp=atof(vs_buf[8].c_str());
                            if(ptmp<0)
                            {
                                printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                                printf("%s\n",buf);
                                exit(EXIT_FAILURE);
                            }
                            if(ptmp==0) {
                                ptmp=__DBL_MIN__;
                                printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                                printf("%s\n",buf);
                            }
                            if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[7]!="NA") || ptmp==1) se.push_back(atof(vs_buf[7].c_str()));
                            else se.push_back(adjSE(betatmp, ptmp));
                        }
                    } else {
                        printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                        printf("%s\n",buf);
                        beta.push_back(-9);
                        se.push_back(-9);
                    }
                    lineNum++;
                    currssize=currs_map.size();
                } else {
                    printf("WARNING: duplicated SNP ID %s in file %s. The duplicated entry is omitted:\n", vs_buf[1].c_str(),esdpath.c_str());
                    printf("%s.\n",buf);
                }
            }
        }
        
        if(gzflag) gzclose(gzfile);
        else flptr.close();
        cout << lineNum << " SNPs to be included from [" + string(esdpath) + "]." << endl;
    }
    void read_smr_sa(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"Chr\"");
            exit(EXIT_FAILURE);
        }
        map<string, int> currs_map;
        long currssize=0;
        while(!flptr.eof() && !gzeof(gzfile))
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=9) {
                    printf("ERROR: the number of columns is incorrect in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                string tmpstr=vs_buf[0]+":"+vs_buf[1];
                currs_map.insert(pair<string,int>(tmpstr,lineNum));
                if(currssize<currs_map.size())
                {
                    if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                        printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpesd.snpchr=23;
                    else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpesd.snpchr=24;
                    else tmpesd.snpchr=atoi(vs_buf[0].c_str());
                    strcpy2(&tmpesd.snprs,vs_buf[1]);
                    if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                        printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    tmpesd.bp=atoi(vs_buf[2].c_str());
                     to_upper(vs_buf[3]);
                     strcpy2(&tmpesd.a1,vs_buf[3]);
                     to_upper(vs_buf[4]);
                    strcpy2(&tmpesd.a2,vs_buf[4]);
                    if(vs_buf[5].compare("NA") && vs_buf[5].compare("na"))
                    {
                        float fq=atof(vs_buf[5].c_str());
                        if(fq<=0 || fq>=1)
                        {
                            printf("ERROR: freq should be between 0 and 1 in row %d.\n", lineNum+2);
                            exit(EXIT_FAILURE);
                        }
                        tmpesd.freq=fq;
                    } else {
                        tmpesd.freq=-9;
                    }
                     to_upper(vs_buf[6]);
                     to_upper(vs_buf[7]);
                    to_upper(vs_buf[8]);
                    if(vs_buf[6].compare("NA") && (vs_buf[7].compare("NA") || vs_buf[8].compare("NA"))){
                        tmpesd.beta=atof(vs_buf[6].c_str());
                        if(vs_buf[8]=="NA") tmpesd.se=atof(vs_buf[7].c_str());
                        else {
                            double betatmp=atof(vs_buf[6].c_str());
                            double ptmp=atof(vs_buf[8].c_str());
                            if(ptmp<0)
                            {
                                printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                                printf("%s\n",buf);
                                exit(EXIT_FAILURE);
                            }
                            if(ptmp==0) {
                                ptmp=__DBL_MIN__;
                                printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                                printf("%s\n",buf);
                            }
                            if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[7]!="NA") || ptmp==1) tmpesd.se=atof(vs_buf[7].c_str());
                            else tmpesd.se=adjSE(betatmp, ptmp);
                        }
                        snpinfo.push_back(tmpesd);
                    } else {
                        printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                        printf("%s\n",buf);
                        tmpesd.beta=-9;
                        tmpesd.se=-9;
                    }
                    lineNum++;
                    currssize=currs_map.size();
                } else {
                    printf("WARNING: duplicated SNP ID %s in file %s. The duplicated entry is omitted:\n", vs_buf[1].c_str(),esdpath.c_str());
                    printf("%s.\n",buf);
                }
            }
        }
        
        if(gzflag) gzclose(gzfile);
        else flptr.close();
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
    }
    void read_plink_qassoc_gz(vector<string> &rs,vector<int> &chr,vector<int> &bp, vector<float> &beta,vector<float> &se,string esdpath){
        
        cout << "Reading summary information from [" + esdpath + "]." << endl;
        char tbuf[MAX_LINE_SIZE];
        gzFile gzfile = gzopen(esdpath.c_str(), "rb");
        if (!(gzfile)) {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     esdpath.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        int lineNum(0);
        gzgets(gzfile, tbuf, MAX_LINE_SIZE); //head
        if(tbuf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(tbuf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        
        while(!gzeof(gzfile))
        {
            gzgets(gzfile, tbuf, MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(tbuf, vs_buf, ", \t\n");
                if(col_num!=9) {
                    printf("ERROR: the number of columns is incorrect in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                rs.push_back(vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }

                bp.push_back(atoi(vs_buf[2].c_str()));
                to_upper(vs_buf[4]);
                to_upper(vs_buf[5]);
                to_upper(vs_buf[8]);
                if(vs_buf[4].compare("NA") && (vs_buf[5].compare("NA") || vs_buf[8].compare("NA")))
                {
                    beta.push_back(atof(vs_buf[4].c_str()));
                    if(vs_buf[8]=="NA") se.push_back(atof(vs_buf[5].c_str()));
                    else {
                        double betatmp=atof(vs_buf[4].c_str());
                        double ptmp=atof(vs_buf[8].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",tbuf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",tbuf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[5]!="NA") || ptmp==1) se.push_back(atof(vs_buf[5].c_str()));
                        else se.push_back(adjSE(betatmp, ptmp));
                    }
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",tbuf);
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
            }
        }
        gzclose(gzfile);
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }
    
    void read_plink_qassoc(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<float> &beta,vector<float> &se,string esdpath){
        
        ifstream flptr(esdpath.c_str());
        if (!flptr) throw ("Error: can not open the file [" + string(esdpath) + "] to read.");
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        flptr.getline(buf,MAX_LINE_SIZE);
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        while(!flptr.eof())
        {
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=9) {
                    printf("ERROR: the number of columns is incorrect in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                rs.push_back(vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                bp.push_back(atoi(vs_buf[2].c_str()));
                 to_upper(vs_buf[4]);
                 to_upper(vs_buf[5]);
                to_upper(vs_buf[8]);
                if(vs_buf[4].compare("NA") && (vs_buf[5].compare("NA") || vs_buf[8].compare("NA")))
                {
                    beta.push_back(atof(vs_buf[4].c_str()));
                    if(vs_buf[8]=="NA") se.push_back(atof(vs_buf[5].c_str()));
                    else {
                        double betatmp=atof(vs_buf[4].c_str());
                        double ptmp=atof(vs_buf[8].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",buf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",buf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[5]!="NA") || ptmp==1) se.push_back(atof(vs_buf[5].c_str()));
                        else se.push_back(adjSE(betatmp, ptmp));
                    }
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",buf);
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
            }
        }
        flptr.close();
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }

    void read_plink_qassoc(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        while(!flptr.eof() && !gzeof(gzfile))
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=9) {
                    printf("ERROR: the number of columns is incorrect in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpesd.snpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpesd.snpchr=24;
                else tmpesd.snpchr=atoi(vs_buf[0].c_str());
                strcpy2(&tmpesd.snprs,vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                tmpesd.bp=atoi(vs_buf[2].c_str());
                 to_upper(vs_buf[4]);
                 to_upper(vs_buf[5]);
                to_upper(vs_buf[8]);
                if(vs_buf[4].compare("NA") && (vs_buf[5].compare("NA") || vs_buf[8].compare("NA")))
                {
                    tmpesd.beta=atof(vs_buf[4].c_str());
                    if(vs_buf[8]=="NA") tmpesd.se=atof(vs_buf[5].c_str());
                    else {
                        double betatmp=atof(vs_buf[4].c_str());
                        double ptmp=atof(vs_buf[8].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",buf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",buf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[5]!="NA") || ptmp==1) tmpesd.se=atof(vs_buf[5].c_str());
                        else tmpesd.se=adjSE(betatmp, ptmp);
                    }
                    snpinfo.push_back(tmpesd);
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",buf);
                   tmpesd.beta=-9;
                   tmpesd.se=-9;
                }
                
                lineNum++;
            }
        }
        if(gzflag) gzclose(gzfile);
        else flptr.close();
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }
    void read_merlin_qassoc_gz(vector<string> &rs,vector<int> &chr,vector<string> &a1, vector<string> &a2, vector<float> freq, vector<float> &beta,vector<float> &se,string esdpath){
        
        cout << "Reading summary information from [" + esdpath + "]." << endl;
        char tbuf[MAX_LINE_SIZE];
        gzFile gzfile = gzopen(esdpath.c_str(), "rb");
        if (!(gzfile)) {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     esdpath.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        int lineNum(0);
        gzgets(gzfile, tbuf, MAX_LINE_SIZE); //head
        if(tbuf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(tbuf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        
        while(!gzeof(gzfile))
        {
            gzgets(gzfile, tbuf, MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(tbuf, vs_buf, ", \t\n");
                if(col_num!=11) {
                    printf("ERROR: the number of columns is incorrect in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                rs.push_back(vs_buf[1]);
                
                to_upper(vs_buf[2]);
                a1.push_back(vs_buf[2]);
                to_upper(vs_buf[3]);
                a2.push_back(vs_buf[3]);
                if(vs_buf[4].compare("NA") || vs_buf[4].compare("na"))
                {
                    float fq=atof(vs_buf[4].c_str());
                    if(fq<=1e-8 || fq>=1)
                    {
                        printf("ERROR: Allele frequency should be between 0 and 1 in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    freq.push_back(fq);
                } else {
                    freq.push_back(-9);
                }

                to_upper(vs_buf[6]);
                to_upper(vs_buf[7]);
                to_upper(vs_buf[10]);
                if(vs_buf[6].compare("NA") && (vs_buf[7].compare("NA") || vs_buf[10].compare("NA")))
                {
                    beta.push_back(atof(vs_buf[6].c_str()));
                    if(vs_buf[10]=="NA") se.push_back(atof(vs_buf[7].c_str()));
                    else {
                        double betatmp=atof(vs_buf[6].c_str());
                        double ptmp=atof(vs_buf[10].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",tbuf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",tbuf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[7]!="NA")  || ptmp==1) se.push_back(atof(vs_buf[7].c_str()));
                        else se.push_back(adjSE(betatmp, ptmp));

                    }
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",tbuf);
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
            }
        }
        gzclose(gzfile);
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }
    void read_gemma_qassoc(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<string> &a1, vector<string> &a2, vector<float> freq, vector<float> &beta,vector<float> &se,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"chr\"");
            exit(EXIT_FAILURE);
        }
        
        while(!flptr.eof() && !gzeof(gzfile))
        {
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=11) {
                    printf("ERROR: the number of columns is incorrect in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                
                rs.push_back(vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                bp.push_back(atoi(vs_buf[2].c_str()));
                 to_upper(vs_buf[4]);
                a1.push_back(vs_buf[4]);
                 to_upper(vs_buf[5]);
                a2.push_back(vs_buf[5]);
                if(vs_buf[6].compare("NA") && vs_buf[6].compare("na"))
                {
                    float fq=atof(vs_buf[6].c_str());
                    if(fq<0 || fq>1)
                    {
                        printf("ERROR: freq should be between 0 and 1 in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    freq.push_back(fq);
                    
                } else {
                    freq.push_back(-9);
                }
                
                 to_upper(vs_buf[7]);
                 to_upper(vs_buf[8]);
                 to_upper(vs_buf[10]);
                if(vs_buf[7].compare("NA") && (vs_buf[8].compare("NA") || vs_buf[10].compare("NA"))) {
                    beta.push_back(atof(vs_buf[7].c_str()));
                    if(vs_buf[10]=="NA") se.push_back(atof(vs_buf[8].c_str()));
                    else {
                        double betatmp=atof(vs_buf[7].c_str());
                        double ptmp=atof(vs_buf[10].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",buf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",buf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[8]!="NA") || ptmp==1) se.push_back(atof(vs_buf[8].c_str()));
                        else se.push_back(adjSE(betatmp, ptmp));
                    }
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",buf);
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
                
            }
        }
        
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
        if(gzflag) gzclose(gzfile);
        else flptr.close();
    }
    void read_bolt_qassoc(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<string> &a1, vector<string> &a2, vector<float> freq, vector<float> &beta,vector<float> &se,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: the headers should start with \"SNP\"");
            exit(EXIT_FAILURE);
        }
        
        while(!flptr.eof() && !gzeof(gzfile))
        {
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=12) {
                    printf("ERROR: the number of columns is incorrect in row %d.\n", lineNum+2);
                    printf("\t%s\n",buf);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[1]=="X" || vs_buf[1]=="x") tmpchr=23;
                else if(vs_buf[1]=="Y" || vs_buf[1]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[1].c_str());
                chr.push_back(tmpchr);
                
                rs.push_back(vs_buf[0]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                bp.push_back(atoi(vs_buf[2].c_str()));
                to_upper(vs_buf[4]);
                a1.push_back(vs_buf[4]);
                to_upper(vs_buf[5]);
                a2.push_back(vs_buf[5]);
                if(vs_buf[6].compare("NA") && vs_buf[6].compare("na"))
                {
                    float fq=atof(vs_buf[6].c_str());
                    if(fq<=1e-8 || fq>=1)
                    {
                        printf("ERROR: Allele frequency should be between 0 and 1 in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    freq.push_back(fq);
                } else {
                    freq.push_back(-9);
                }
                to_upper(vs_buf[8]);
                to_upper(vs_buf[9]);
                to_upper(vs_buf[11]);
                if(vs_buf[8].compare("NA") && (vs_buf[9].compare("NA") || vs_buf[11].compare("NA"))) {
                    beta.push_back(atof(vs_buf[8].c_str()));
                    if(vs_buf[11]=="NA") se.push_back(atof(vs_buf[9].c_str()));
                    else {
                        double betatmp=atof(vs_buf[8].c_str());
                        double ptmp=atof(vs_buf[11].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",buf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",buf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[9]!="NA") || ptmp==1) se.push_back(atof(vs_buf[9].c_str()));
                        else se.push_back(adjSE(betatmp, ptmp));
                    }
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",buf);
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
                
            }
        }
        
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
        if(gzflag) gzclose(gzfile);
        else flptr.close();
    }
    void read_gemma_qassoc(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"chr\"");
            exit(EXIT_FAILURE);
        }
        
        while(!flptr.eof() && !gzeof(gzfile))
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=11) {
                    printf("ERROR: column number is not correct in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpesd.snpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpesd.snpchr=24;
                else tmpesd.snpchr=atoi(vs_buf[0].c_str());
                strcpy2(&tmpesd.snprs,vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                tmpesd.bp=atoi(vs_buf[2].c_str());
                 to_upper(vs_buf[4]);
                strcpy2(&tmpesd.a1,vs_buf[4]);
                 to_upper(vs_buf[5]);
                strcpy2(&tmpesd.a2,vs_buf[5]);
                if(vs_buf[6].compare("NA") && vs_buf[6].compare("na"))
                {
                    float fq=atof(vs_buf[6].c_str());
                    if(fq<=1e-8 || fq>1)
                    {
                        printf("ERROR: allele frequency should be between 0 and 1 in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    
                    tmpesd.freq=fq;
                } else {
                    tmpesd.freq=-9;
                }
                 to_upper(vs_buf[7]);
                 to_upper(vs_buf[8]);
                to_upper(vs_buf[10]);
                if(vs_buf[7].compare("NA") && (vs_buf[8].compare("NA") || vs_buf[10].compare("NA"))) {
                    tmpesd.beta=atof(vs_buf[7].c_str());
                    if(vs_buf[10]=="NA") tmpesd.se=atof(vs_buf[8].c_str());
                    else {
                        double betatmp=atof(vs_buf[7].c_str());
                        double ptmp=atof(vs_buf[10].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",buf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",buf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[8]!="NA") || ptmp==1) tmpesd.se=atof(vs_buf[8].c_str());
                        else tmpesd.se=adjSE(betatmp, ptmp);
                    }
                    snpinfo.push_back(tmpesd);
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",buf);
                    tmpesd.beta=-9;
                    tmpesd.se=-9;
                }
               
                lineNum++;
            }
        }
        
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
        if(gzflag) gzclose(gzfile);
        else flptr.close();
    }
    void read_bolt_qassoc(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",esdpath.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: the headers should start with \"SNP\"");
            exit(EXIT_FAILURE);
        }
        
        while(!flptr.eof() && !gzeof(gzfile))
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=12) {
                    printf("ERROR: the number of columns is incorrect in row %d.\n", lineNum+2);
                    printf("\t%s\n",buf);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the chromosome is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[1]=="X" || vs_buf[1]=="x") tmpchr=23;
                else if(vs_buf[1]=="Y" || vs_buf[1]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[1].c_str());
                tmpesd.snpchr=tmpchr;
                strcpy2(&tmpesd.snprs,vs_buf[0]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: the SNP position is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                tmpesd.bp=atoi(vs_buf[2].c_str());
                to_upper(vs_buf[4]);
                strcpy2(&tmpesd.a1,vs_buf[4]);
                to_upper(vs_buf[5]);
                strcpy2(&tmpesd.a2,vs_buf[5]);
                if(vs_buf[6].compare("NA") && vs_buf[6].compare("na"))
                {
                    float fq=atof(vs_buf[6].c_str());
                    if(fq<=1e-8 || fq>1)
                    {
                        printf("ERROR: allele frequency should be between 0 and 1 in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    
                    tmpesd.freq=fq;
                } else {
                    tmpesd.freq=-9;
                }
                to_upper(vs_buf[8]);
                to_upper(vs_buf[9]);
                to_upper(vs_buf[11]);
                if(vs_buf[8].compare("NA") && (vs_buf[9].compare("NA") || vs_buf[11].compare("NA"))) {
                    tmpesd.beta=atof(vs_buf[8].c_str());
                    if(vs_buf[11]=="NA") tmpesd.se=atof(vs_buf[9].c_str());
                    else {
                        double betatmp=atof(vs_buf[8].c_str());
                        double ptmp=atof(vs_buf[11].c_str());
                        if(ptmp<0)
                        {
                            printf("ERROR: p-value should be positive in row %d.\n",lineNum+2);
                            printf("%s\n",buf);
                            exit(EXIT_FAILURE);
                        }
                        if(ptmp==0) {
                            ptmp=__DBL_MIN__;
                            printf("WARNING: p-value of 0 found in row %d and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                            printf("%s\n",buf);
                        }
                        if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[9]!="NA") || ptmp==1) tmpesd.se=atof(vs_buf[9].c_str());
                        else tmpesd.se=adjSE(betatmp, ptmp);
                    }
                    snpinfo.push_back(tmpesd);
                } else {
                    printf("WARNING: this row is omitted because beta is missing or both SE and P are missing (\"NA\").\n");
                    printf("%s\n",buf);
                    tmpesd.beta=-9;
                    tmpesd.se=-9;
                }
                lineNum++;
                
            }
        }
        
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
        if(gzflag) gzclose(gzfile);
        else flptr.close();
    }
    
    void read_mapfile(bInfo* bdata,string bimfile) {
        // Read mapfile in the format as: chr snp bp
        int ibuf = 0;
        string cbuf = "0";
        double dbuf = 0.0;
        string str_buf;
        ifstream Bim(bimfile.c_str());
        if (!Bim) throw ("Error: can not open the file [" + bimfile + "] to read.");
        cout << "Reading PLINK BIM file from [" + bimfile + "]." << endl;
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        while (Bim) {
            Bim >> ibuf;
            if (Bim.eof()) break;
            bdata->_chr.push_back(ibuf);
            Bim >> str_buf;
            bdata->_snp_name.push_back(str_buf);
            Bim >> ibuf;
            bdata->_bp.push_back(ibuf);
        }
        Bim.close();
        bdata->_snp_num = (int)bdata->_chr.size();
        cout << bdata->_snp_num << " SNPs to be included from [" + bimfile + "]." << endl;
        
        // Initialize _include
        bdata->_include.clear();
        bdata->_include.resize( bdata->_snp_num);
        bdata->_snp_name_map.clear();
        
        for (int i = 0; i <  bdata->_snp_num; i++) {
            
            bdata->_include[i] = i;
            if( bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()){
                cout << "Warning: Duplicated SNP ID \"" + bdata->_snp_name[i] + "\" ";
                stringstream ss;
                ss << bdata->_snp_name[i] << "_" << i + 1;
                bdata->_snp_name[i] = ss.str();
                cout<<"has been changed to \"" + bdata->_snp_name[i] + "\".\n";
            }
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }

    uint64_t get_esi_info(vector<snpinfolst> &snpinfo, probeinfolst** prbiflst,int prb_num,int fformat)
    {
        probeinfolst* locinfolst=*prbiflst;
        map<string, int> rs_map;
        //map<string, int> rsbp_map;
        //map<string, int> rsaa_map;
        long size = 0;
        uint64_t ttl_value_num =0;
        char buf[MAX_LINE_SIZE];
        double disp=0;
        for(int j=0;j<prb_num;j++)
        {
            progress(j, disp, (int)prb_num);
       
            //map<string, int> currs_map;
            //long currssize=0;
            string esdfilename=(locinfolst+j)->esdpath;
            gzFile gzfile = NULL;
            ifstream flptr;
            bool gzflag=has_suffix(esdfilename, "gz");
            if(gzflag)
            {
                gzfile = gzopen(esdfilename.c_str(), "rb");
                if (!(gzfile)) {
                    fprintf (stderr, "%s: Couldn't open file %s\n",
                             esdfilename.c_str(), strerror (errno));
                    exit (EXIT_FAILURE);
                }
            }
            else
            {
                flptr.open((locinfolst+j)->esdpath);
                if (!flptr) throw ("Error: can not open the file [" + string((locinfolst+j)->esdpath) + "] to read.");
            }
            
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
            else flptr.getline(buf,MAX_LINE_SIZE);// the header
            vector<string> vs_buf;
            if(buf[0]=='\0')
            {
                printf("ERROR: the first row of the file %s is empty.\n",esdfilename.c_str());
                exit(EXIT_FAILURE);
            }
            int col_num = split_string(buf, vs_buf, ", \t\n");
            to_upper(vs_buf[0]);
            if(fformat!=3 && vs_buf[0]!="CHR") {
                printf("ERROR: the headers of file %s should start with \"Chr\".\n",esdfilename.c_str());
                exit(EXIT_FAILURE);
            }
            if(fformat==3 && vs_buf[0]!="SNP") {
                printf("ERROR: the headers of file %s should start with \"SNP\".\n",esdfilename.c_str());
                exit(EXIT_FAILURE);
            }
            
            if(fformat==0)
            {
                int lineNum=0;
                while(!flptr.eof() && !gzeof(gzfile))
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        ttl_value_num++;
                        vs_buf.clear();
                        col_num = split_string(buf, vs_buf, ", \t\n");
                        if(col_num!=9) {
                            printf("ERROR: the number of columns is incorrect of row %d in esd file \"%s\"!\n", lineNum+2, esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        
                        if(vs_buf[0]=="NA" ||vs_buf[0]=="na" || vs_buf[1]=="NA" ||vs_buf[1]=="na" ||vs_buf[2]=="NA" ||vs_buf[2]=="na")
                        {
                            printf("ERROR: chromosome or SNP name or SNP position of row %d can't be \"NA\" in esd file \"%s\"!\n", lineNum+2, esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        /*
                        currs_map.insert(pair<string, int>(vs_buf[1].c_str(), lineNum));
                        if(currs_map.size()==currssize)
                        {
                            printf("WARING: duplicate SNP %s found in file %s.\n", vs_buf[1].c_str(),esdfilename.c_str());
                        } else {
                            currssize=currs_map.size();
                        }
                         */
                         to_upper(vs_buf[3]);
                         to_upper(vs_buf[4]);
                        if(vs_buf[3]=="NA" ||vs_buf[4]=="NA")
                        {
                            printf("WARNING: Allele infomation of row %d is \"NA\" in esd file \"%s\"! This may cause incorrect positive/negative sign of effect size when conducting allele check.\n", lineNum+2, esdfilename.c_str());
                           
                        }
                        string crsstr=vs_buf[0]+":"+vs_buf[1];
                        //string crsbpstr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2];
                        //string crsbpaastr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2]+":"+vs_buf[3]+":"+vs_buf[4];
                        //string crsbpaastr_iv=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2]+":"+vs_buf[4]+":"+vs_buf[3];
                        rs_map.insert(pair<string, int>(crsstr.c_str(), lineNum));
                        //rsbp_map.insert(pair<string, int>(crsbpstr.c_str(), lineNum));
                        //rsaa_map.insert(pair<string, int>(crsbpaastr.c_str(), lineNum));
                        //if(rs_map.size() != rsbp_map.size())
                        //{
                        //    printf("ERROR: SNP %s on Chromosome %s has multiple BPs, please check.\n", vs_buf[1].c_str(), vs_buf[0].c_str());
                        //    exit(EXIT_FAILURE);
                        //}
                        bool newsnp=false;
                        if (size < rs_map.size())
                        {
                            newsnp=true;
                            size = rs_map.size();
                        }/*else {
                            //check multi-allelic SNPs
                            if(rs_map.size() != rsaa_map.size())
                            {
                                long tmpsize=rsaa_map.size();
                                rsaa_map.insert(pair<string, int>(crsbpaastr_iv.c_str(), lineNum));
                                if(tmpsize==rsaa_map.size())
                                {
                                    printf("WARING: switched ref allele with alt allele of SNP %s found.\n", vs_buf[1].c_str());
                                } else {
                                    newsnp=true;
                                    printf("WARING: multi-allelic SNPs with duplicate SNP ID %s found.\n", vs_buf[1].c_str());
                                    rs_map.insert(pair<string, int>((crsbpaastr+"_m").c_str(), lineNum));
                                    rsbp_map.insert(pair<string, int>((crsbpaastr+"_m").c_str(), lineNum));
                                }
                                rs_map.insert(pair<string, int>((crsbpaastr+"_").c_str(), lineNum));
                                rsbp_map.insert(pair<string, int>((crsbpaastr+"_").c_str(), lineNum));
                                size = rs_map.size();
                            }
                        }*/
                        if(newsnp)
                        {
                            snpinfolst snptmp;
                            int tmpchr;
                            if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                            else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                            else tmpchr=atoi(vs_buf[0].c_str());
                            snptmp.snpchr=tmpchr;
                            strcpy2(&snptmp.snprs, vs_buf[1]);
                            snptmp.bp=atoi(vs_buf[2].c_str());
                            strcpy2(&snptmp.a1, vs_buf[3]);
                            strcpy2(&snptmp.a2, vs_buf[4]);
                            if(vs_buf[5]=="NA" || vs_buf[5]=="na" )
                            {
                                snptmp.freq=-9;
                            } else {
                                snptmp.freq=atof(vs_buf[5].c_str());
                            }
                            snptmp.gd=0;
                            snpinfo.push_back(snptmp);
                        }
                        lineNum++;
                    }
                }
            }
            else if(fformat==1)
            {
                bInfo binfo;
                string bimfname=string((locinfolst+j)->bfilepath)+".bim";
                read_bimfile(&binfo, bimfname);
                
                while(!flptr.eof() && !gzeof(gzfile))
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        string tmpStr;
                        int tmpchr;
                        snpinfolst snptmp;
                        ttl_value_num++;
                        istringstream iss(buf);
                        iss>>tmpStr; //chr
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: chromosome  can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                        else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                        else tmpchr=atoi(tmpStr.c_str());
                        iss>>tmpStr; //rs
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: SNP name  can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }

                        if(tmpStr[0]!='\0') rs_map.insert(pair<string, int>(tmpStr.c_str(), tmpchr));
                        if (size < rs_map.size())
                        {
                            snptmp.snpchr=tmpchr;
                            strcpy2(&snptmp.snprs, tmpStr);
                            map<string, int>::iterator iter;
                            iter=binfo._snp_name_map.find(tmpStr);
                            if( iter != binfo._snp_name_map.end())
                            {
                                strcpy2(&snptmp.a1, binfo._allele1[iter->second]);
                                strcpy2(&snptmp.a2, binfo._allele2[iter->second]);
                            }else {
                                printf("ERROR: can not find SNP %s in the bim file %s.\n", tmpStr.c_str(), bimfname.c_str());
                                exit(EXIT_FAILURE);
                            }
                            iss>>tmpStr; //BP
                            if(tmpStr=="NA" ||tmpStr=="na" )
                            {
                                printf("ERROR: SNP BP is \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                                exit(EXIT_FAILURE);
                            }

                            snptmp.bp=atoi(tmpStr.c_str());
                            snptmp.gd=0;
                            snptmp.freq=-9;
                            size = rs_map.size();
                            snpinfo.push_back(snptmp);
                        }
                    }
                }
            }
            else if(fformat==2)
            {
                while(!flptr.eof() && !gzeof(gzfile) )
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        string tmpStr;
                        int tmpchr;
                        snpinfolst snptmp;
                        ttl_value_num++;
                        istringstream iss(buf);
                        iss>>tmpStr; //chr
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: chromosome can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                        else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                        else tmpchr=atoi(tmpStr.c_str());
                        iss>>tmpStr; //rs
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: SNP name can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr[0]!='\0') rs_map.insert(pair<string, int>(tmpStr.c_str(), tmpchr));
                        if (size < rs_map.size())
                        {
                            snptmp.snpchr=tmpchr;
                            strcpy2(&snptmp.snprs, tmpStr);
                            iss>>tmpStr; //BP
                            if(tmpStr=="NA" ||tmpStr=="na" )
                            {
                                printf("ERROR: SNP position can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                                exit(EXIT_FAILURE);
                            }
                            snptmp.bp=atoi(tmpStr.c_str());
                            iss>>tmpStr;//n_miss
                            iss>>tmpStr; //A1 minor allele
                             to_upper(tmpStr);
                            strcpy2(&snptmp.a1, tmpStr);
                            iss>>tmpStr; //A0 major allele
                             to_upper(tmpStr);
                            strcpy2(&snptmp.a2, tmpStr);
                            iss>>tmpStr; //freq
                            if(tmpStr=="NA" || tmpStr=="na" )
                            {
                                snptmp.freq=-9;
                            } else {
                                snptmp.freq=atof(tmpStr.c_str());
                            }
                            snptmp.gd=0;
                            size = rs_map.size();
                            snpinfo.push_back(snptmp);
                        }
                    }
                }
            }
            else if(fformat==3)
            {
                while(!flptr.eof() && !gzeof(gzfile) )
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        string tmpStr;
                        int tmpchr;
                        snpinfolst snptmp;
                        ttl_value_num++;
                        istringstream iss(buf);
                        iss>>tmpStr; //rs
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: SNP name can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr[0]!='\0') rs_map.insert(pair<string, int>(tmpStr.c_str(), tmpchr));
                        if (size < rs_map.size())
                        {
                            strcpy2(&snptmp.snprs, tmpStr);
                            iss>>tmpStr; //chr
                            if(tmpStr=="NA" ||tmpStr=="na" )
                            {
                                printf("ERROR: chromosome can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                                exit(EXIT_FAILURE);
                            }
                            if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                            else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                            else tmpchr=atoi(tmpStr.c_str());
                            snptmp.snpchr=tmpchr;
                            iss>>tmpStr; //BP
                            if(tmpStr=="NA" ||tmpStr=="na" )
                            {
                                printf("ERROR: SNP position is \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                                exit(EXIT_FAILURE);
                            }
                            snptmp.bp=atoi(tmpStr.c_str());
                            iss>>tmpStr;//genpos
                            iss>>tmpStr; //A0 minor allele
                            to_upper(tmpStr);
                            strcpy2(&snptmp.a1, tmpStr);
                            iss>>tmpStr; //A1 major allele
                            to_upper(tmpStr);
                            strcpy2(&snptmp.a2, tmpStr);
                            iss>>tmpStr; //freq
                            if(tmpStr=="NA" || tmpStr=="na" )
                            {
                                snptmp.freq=-9;
                            } else {
                                snptmp.freq=atof(tmpStr.c_str());
                            }
                            snptmp.gd=0;
                            size = rs_map.size();
                            snpinfo.push_back(snptmp);
                        }
                    }
                }
            }
            else if(fformat==4)
            {
                //NOTE: only for personal use with BSGS data.
                bInfo binfo;
                string bimfname=string((locinfolst+j)->bfilepath);
                read_mapfile(&binfo, bimfname);
                while(!flptr.eof() && !gzeof(gzfile))
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        string tmpStr;
                        int tmpchr;
                        string a1;
                        string a2;
                        string frqstr;
                        float frq;
                        
                        snpinfolst snptmp;
                        ttl_value_num++;
                        istringstream iss(buf);
                        iss>>tmpStr; //chr
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: chromosome  can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                        else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                        else tmpchr=atoi(tmpStr.c_str());
                        iss>>tmpStr; //rs
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: SNP name  can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        iss >> a1;
                        to_upper(a1);
                        iss >> a2;
                        to_upper(a2);
                        iss >> frqstr;
                        if(frqstr=="NA" ||frqstr=="na" )
                        {
                            frq=atof(frqstr.c_str());
                            if(frq<=1e-8 || frq>1)
                            {
                                printf("ERROR: allele frequency should be between 0 and 1 in file %s.\n", esdfilename.c_str());
                                exit(EXIT_FAILURE);
                            }
                        } else frq=-9;

                        if(tmpStr[0]!='\0') rs_map.insert(pair<string, int>(tmpStr.c_str(), tmpchr));
                        if (size < rs_map.size())
                        {
                            snptmp.snpchr=tmpchr;
                            strcpy2(&snptmp.snprs, tmpStr);
                            strcpy2(&snptmp.a1, a1);
                            strcpy2(&snptmp.a2, a2);
                            map<string, int>::iterator iter;
                            iter=binfo._snp_name_map.find(tmpStr);
                            if( iter != binfo._snp_name_map.end())
                            {
                                snptmp.bp=binfo._bp[iter->second];
                            }else {
                                printf("ERROR: can not find SNP %s in the bim file %s.\n", tmpStr.c_str(), bimfname.c_str());
                                exit(EXIT_FAILURE);
                            }
                            snptmp.gd=0;
                            snptmp.freq=frq;
                            size = rs_map.size();
                            snpinfo.push_back(snptmp);
                        }
                    }
                }
            }
        }
        return ttl_value_num;
    }
    void save_txts_dbesd(char* outFileName, long esiNum, long epiNum,vector<int> &epi2esd,probeinfolst* prbiflst, int fformat,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2, int addn)
    {
        // get esd info
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        uint32_t filetype=DENSE_FILE_TYPE_3;
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        if(addn!=-9)
            printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
        ten_ints[1]=addn;
        ten_ints[2]=(int)esiNum;
        ten_ints[3]=(int)epiNum;
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
        
        uint64_t bsize=(uint64_t)esiNum<<1;
        float* buffer=(float*)malloc (sizeof(float)*bsize);
        if (NULL == buffer) {
            printf("ERROR: failed to allocate write buffer for file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        map<string, int>::iterator iter;
        double disp=0;
        for(int j=0;j<epiNum;j++)
        {
              progress(j, disp, (int)epiNum);
            int esdnumcurprb=epi2esd[j+1]-epi2esd[j];
            printf("Reading the probe %s from %d text file(s).\n",prbiflst[epi2esd[j]].probeId, esdnumcurprb);
            for(int k=0;k<bsize;k++) buffer[k]=-9; //init
            for(int k=epi2esd[j];k<epi2esd[j+1];k++)
            {
                vector<string> _rs;
                vector<float> _beta;
                vector<float> _se;
                vector<int> _chr;
                vector<string> _a1;
                vector<string> _a2;
                vector<int> _bp;
                vector<float> _freq;
                string zname=prbiflst[k].esdpath;
                bool gzflag=has_suffix(zname, "gz");
                switch (fformat)
                {
                    case 0:
                        read_smr_sa(_rs,_chr,_bp, _a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 1:
                        if(gzflag) read_plink_qassoc_gz(_rs, _chr,_bp, _beta, _se, prbiflst[k].esdpath);
                        else read_plink_qassoc(_rs,_chr, _bp,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 2:
                        read_gemma_qassoc(_rs,_chr, _bp,_a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 3:
                        read_bolt_qassoc(_rs,_chr, _bp,_a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 4:
                        read_merlin_qassoc_gz(_rs,_chr,_a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                }
               // when users using --geno-uni, even the SNPs in each files are the same, but can't guarantee in the same order. so this step is mandatory
                vector<int> rsid(_rs.size());
                for (int l = 0; l<_rs.size(); l++){
                    iter = esi_map.find(_rs[l]);
                    if (iter != esi_map.end()) rsid[l]=iter->second;
                    else {
                        printf("ERROR: SNP %s is not in output SNP set. If you are using --geno-uni, please disable it then try again. Otherwise please report this bug.\n",_rs[l].c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                for(int l=0;l<rsid.size();l++)
                {
                    if(fformat!=1)
                    {
                        if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l]){
                            buffer[rsid[l]]=_beta[l];
                            buffer[esiNum+rsid[l]]=_se[l];
                        } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l]){
                            buffer[rsid[l]]=-1.0*_beta[l];
                            buffer[esiNum+rsid[l]]=_se[l];
                        } else {
                            printf("ERROR: inconsistent allele pairs of SNP %s found.\n",_rs[l].c_str());
                            printf("Discrepant Allele pairs: (%s,%s) with (%s,%s).\n",esi_a1[rsid[l]].c_str(), esi_a2[rsid[l]].c_str(),_a1[l].c_str(), _a2[l].c_str());
                            exit(EXIT_FAILURE);
                            
                            //this part is for multi-allelic SNPs. since we don't save multi-allelic SNPs anymore, so we should disable it.
                            /*
                            int did=-9;
                            float sig=1.0;
                            for(int m=0;m<esi_rs.size();m++)
                            {
                                if(esi_rs[m]==_rs[l])
                                {
                                    if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                    {
                                        did=m;
                                        break;
                                    }
                                    if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                    {
                                        did=m;
                                        sig=-1.0;
                                        break;
                                    }
                                }
                            }
                            if(did==-9)
                            {
                                printf("ERROR: This would not go to happen. Please report this bug.");
                                exit(EXIT_FAILURE);
                            }
                            buffer[did]=sig*_beta[l];
                            buffer[esiNum+did]=_se[l];
                             */
                        }
                    } else {
                        buffer[rsid[l]]=_beta[l];
                        buffer[esiNum+rsid[l]]=_se[l];
                    }
                }
            }
            fwrite (buffer,sizeof(float), bsize, smr1);
        }
        fclose (smr1);
        free(buffer);
        cout<<"Effect sizes (beta) and SE for "<<epiNum<<" probes and "<<esiNum<<" SNPs have been saved in a binary file [" + esdfile + "]." <<endl;

    }
    
     void save_full_txts_sbesd(char* outFileName, long esiNum, long epiNum,vector<int> &epi2esd,probeinfolst* prbiflst, int fformat,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2, int addn)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        // get esd info
        
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        
        uint32_t filetype=SPARSE_FILE_TYPE_3;
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        if(addn!=-9)
            printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
        ten_ints[1]=addn;
        ten_ints[2]=(int)esiNum;
        ten_ints[3]=(int)epiNum;
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
        
        vector<uint64_t> cols((epiNum<<1)+1);;
        vector<uint32_t> rowids;
        vector<float> val;
        cols[0]=0;
        
        map<string, int>::iterator iter;
        double disp=0;
        for(int j=0;j<epiNum;j++)
        {
            progress(j, disp, (int)epiNum);
            int esdnumcurprb=epi2esd[j+1]-epi2esd[j];
            printf("\nReading %d text file(s) of probe %s.\n", esdnumcurprb, prbiflst[epi2esd[j]].probeId);
            
            map<string, int> rsa_map;
            long rsNum=0;
            
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            for(int k=epi2esd[j];k<epi2esd[j+1];k++)
            {
                vector<string> _rs;
                vector<float> _beta;
                vector<float> _se;
                vector<int> _chr;
                vector<string> _a1;
                vector<string> _a2;
                vector<int> _bp;
                vector<float> _freq;
                string zname=prbiflst[k].esdpath;
                bool gzflag=has_suffix(zname, "gz");
                switch (fformat)
                {
                    case 0:
                        read_smr_sa(_rs,_chr,_bp, _a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 1:
                        if(gzflag) read_plink_qassoc_gz(_rs, _chr,_bp, _beta, _se, prbiflst[k].esdpath);
                        else read_plink_qassoc(_rs,_chr, _bp,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 2:
                        read_gemma_qassoc(_rs, _chr, _bp, _a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 3:
                        read_bolt_qassoc(_rs,_chr, _bp,_a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 4:
                        read_merlin_qassoc_gz(_rs,_chr,_a1, _a2, _freq,_beta, _se, prbiflst[k].esdpath);
                        break;
                }
                
                vector<int> rsid(_rs.size());
                for (int l = 0; l<_rs.size(); l++){
                    iter = esi_map.find(_rs[l]);
                    if (iter != esi_map.end()) rsid[l]=iter->second;
                    else {
                        printf("ERROR: SNP %s is not in output SNP set. If you are using --geno-uni, please disable it then try again. Otherwise please report this bug.\n",_rs[l].c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                for(int l=0;l<rsid.size();l++)
                {
                    if(fabs(_se[l]+9)>1e-6)
                    {
                        string chckstr=(fformat==1?_rs[l]:(_rs[l]+":"+_a1[l]+":"+_a2[l]));
                        rsa_map.insert(pair<string,int>(chckstr,l));
                        if(rsNum<rsa_map.size())
                        {
                            if(fformat!=1)
                            {
                                if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l] ){
                                    val.push_back(_beta[l]);
                                    rowids.push_back(rsid[l]);
                                    tmpse.push_back(_se[l]);
                                    tmprid.push_back(rsid[l]);
                                } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l] ){
                                    val.push_back(-1.0*_beta[l]);
                                    rowids.push_back(rsid[l]);
                                    tmpse.push_back(_se[l]);
                                    tmprid.push_back(rsid[l]);
                                } else {
                                    printf("ERROR: inconsistent allele pairs of SNP %s found.\n",_rs[l].c_str());
                                    printf("Discrepant Allele pairs: (%s,%s) with (%s,%s).\n",esi_a1[rsid[l]].c_str(), esi_a2[rsid[l]].c_str(),_a1[l].c_str(), _a2[l].c_str());
                                    exit(EXIT_FAILURE);
                                    
                                    //this part is for multi-allelic SNPs. since we don't save multi-allelic SNPs anymore, so we should disable it.
                                    /*
                                        int did=-9;
                                        float sig=1.0;
                                        for(int m=0;m<esi_rs.size();m++)
                                        {
                                            if(esi_rs[m]==_rs[l])
                                            {
                                                if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                                {
                                                    did=m;
                                                    break;
                                                }
                                                if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                                {
                                                    did=m;
                                                    sig=-1.0;
                                                    break;
                                                }
                                            }
                                        }
                                        if(did==-9)
                                        {
                                            printf("ERROR: This would not go to happen. Please report this bug.");
                                            exit(EXIT_FAILURE);
                                        }
                                        val.push_back(sig*_beta[l]);
                                        rowids.push_back(did);
                                        tmpse.push_back(_se[l]);
                                        tmprid.push_back(did);
                                     */
                                }
                            } else {
                                    val.push_back(_beta[l]);
                                    rowids.push_back(rsid[l]);
                                    tmpse.push_back(_se[l]);
                                    tmprid.push_back(rsid[l]);
                            }
                            rsNum=rsa_map.size();
                        } else {
                            printf("WARNING: duplicated SNP ID %s for the probe %s. The duplicated entry is omitted.\n",_rs[l].c_str(),prbiflst[k].probeId);
                        }
                    } else {
                        printf("WARNING: SNP %s  with \"NA\" value is found and skipped.\n",_rs[l].c_str());
                    }
                   
                    
                }
            }
            for(int k=0;k<tmpse.size();k++)
            {
                val.push_back(tmpse[k]);
                rowids.push_back(tmprid[k]);
            }
            uint64_t real_num=tmpse.size();
            cols[(j<<1)+1]=real_num+cols[j<<1];
            cols[j+1<<1]=(real_num<<1)+cols[j<<1];
        }
        uint64_t valNum=val.size();
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
        fwrite (&val[0],sizeof(float), val.size(), smr1);
        fclose (smr1);
        
        cout<<"Effect sizes (beta) and SE for "<<epiNum<<" probes and "<<esiNum<<" SNPs have been saved in a binary file [" + esdfile + "]." <<endl;
    }

    void slct_sparse_per_prb(vector<int> &slct_idx, probeinfolst* prbifo, vector<snpinfolst> &snpinfo, long cis_itvl, long trans_itvl,double transThres,double restThres,FILE* logfile, bool extract_cis_only, bool techHit)
    {
        // no null se in snpinfo has been guaranteed before here.
        // for log file
        map<string, int> rsa_map;
        map<string, int>::iterator iter;
        long rsNum=0;

        vector<int> tran_chr;
        vector<uint64_t> lowerBp;
        vector<uint64_t> upperBp;
        vector<uint64_t> nsnp;
        vector<bool> extend;
        vector<bool> merge;
        vector<int> cis_idx;
        map<string, int> other_idx;
        
        //extract info
        string probid=prbifo->probeId;
        long probbp=prbifo->bp;
        long probchr=prbifo->probechr;
        long cisuperBounder= probbp+cis_itvl;
        long cislowerBounder=((probbp-cis_itvl>0)?(probbp-cis_itvl):0);
       
        //printf("Extracting information of the probe %s...\n", probid.c_str());
        for(int l=0;l<snpinfo.size();l++)
        {
            if(fabs(snpinfo[l].se+9)>1e-6)
            {
                double zsxz=snpinfo[l].beta/snpinfo[l].se;
                double pxz=pchisq(zsxz*zsxz, 1);
                
                if(snpinfo[l].snpchr == probchr && snpinfo[l].bp<=cisuperBounder && snpinfo[l].bp>=cislowerBounder)
                {
                    if(techHit)
                    {
                        //printf("The following SNP in the cis-region is excluded due to the technical eQTL.\n");
                        double z=(snpinfo[l].beta/snpinfo[l].se);
                        double p=pchisq(z*z, 1);
                        string tmp=atos(snpinfo[l].snprs)+"\t"+ atos(snpinfo[l].snpchr)+"\t"+ atos(snpinfo[l].bp)+"\t"+ atos(snpinfo[l].a1)+"\t"+ atos(snpinfo[l].a2)+"\t"+ atos(snpinfo[l].freq)+"\t"+ atos(prbifo->probeId)+"\t"+ atos(prbifo->probechr)+"\t"+ atos(prbifo->bp)+"\t" + atos(prbifo->genename)+"\t"+ atos(prbifo->orien)+"\t"+ atos(snpinfo[l].beta)+"\t"+ atos(snpinfo[l].se)+"\t"+ dtos(p)+"\n";
                        if(techeQTLfile) {
                            fputs(tmp.c_str(),techeQTLfile);
                            fflush(techeQTLfile);
                        }
                        
                        printf("%s\n",tmp.c_str());
                        
                    } else {
                        string chckstr=string(snpinfo[l].snprs);
                        //if(snpinfo[l].a1 && snpinfo[l].a2) chckstr += ":"+string(snpinfo[l].a1)+":"+string(snpinfo[l].a2);
                        rsa_map.insert(pair<string,int>(chckstr,rsNum));
                        if(rsNum<rsa_map.size()){
                            slct_idx.push_back(l);
                            rsNum=rsa_map.size();
                            cis_idx.push_back(l);
                        }

                    }
                    
                }
                else if(!extract_cis_only && pxz<=transThres)
                {
                    uint64_t transNum=0;
                    int curChr=snpinfo[l].snpchr;
                    string chckstr=string(snpinfo[l].snprs);
                    //if(snpinfo[l].a1 && snpinfo[l].a2) chckstr += ":"+string(snpinfo[l].a1)+":"+string(snpinfo[l].a2);
                    rsa_map.insert(pair<string, int>(chckstr,rsNum));
                    if(rsNum<rsa_map.size()){
                        slct_idx.push_back(l);
                        rsNum=rsa_map.size();
                        transNum++;
                    }
                    
                    long transbp=snpinfo[l].bp;
                    long translowerBounder=((transbp-trans_itvl>0)?(transbp-trans_itvl):0);
                    long transuperBounder=transbp+trans_itvl;
                    bool extended=false;
                    bool merged=false;
                    
                    int startptr=l-1;
                    while(startptr>=0 && curChr == snpinfo[startptr].snpchr && transbp-snpinfo[startptr].bp<=trans_itvl)
                    {
                        if(fabs(snpinfo[startptr].se+9)>1e-6)
                        {
                            translowerBounder=snpinfo[startptr].bp;
                            if(upperBp.size()>0 && translowerBounder<=upperBp[upperBp.size()-1] && tran_chr[upperBp.size()-1]==snpinfo[startptr].snpchr) //trans region merges
                            {
                                merged=true;
                                transNum=transNum+nsnp[nsnp.size()-1];
                                break;
                            }
                            if(snpinfo[startptr].snpchr == probchr && translowerBounder<=cisuperBounder && translowerBounder>=cislowerBounder) // trans touches cis region
                            {
                                break;
                            }
                            string chckstr=string(snpinfo[startptr].snprs);
                            //if(snpinfo[startptr].a1 && snpinfo[startptr].a2) chckstr += ":"+string(snpinfo[startptr].a1)+":"+string(snpinfo[startptr].a2);
                            iter=other_idx.find(chckstr);
                            if(iter!=other_idx.end()) // trans merges other eqtl
                            {
                                other_idx.erase(iter->first);
                                transNum++;
                            } else {
                                rsa_map.insert(pair<string, int>(chckstr,rsNum));
                                if(rsNum<rsa_map.size()){
                                    slct_idx.push_back(startptr);
                                    rsNum=rsa_map.size();
                                    transNum++;
                                }
                            }
                        }
                        
                        startptr--;
                    }
                    startptr=l+1;
                    while (startptr<snpinfo.size() && curChr == snpinfo[startptr].snpchr && snpinfo[startptr].bp - transbp <= trans_itvl)
                    {
                        if(fabs(snpinfo[startptr].se+9)>1e-6)
                        {
                            transuperBounder=snpinfo[startptr].bp;
                            if(snpinfo[startptr].snpchr == probchr && transuperBounder>=cislowerBounder && transuperBounder<=cisuperBounder) // trans touches cis region
                            {
                                break;
                                
                            }  else {
                                
                                string chckstr=string(snpinfo[startptr].snprs);
                                //if(snpinfo[startptr].a1 && snpinfo[startptr].a2) chckstr += ":"+string(snpinfo[startptr].a1)+":"+string(snpinfo[startptr].a2);
                                rsa_map.insert(pair<string, int>(chckstr,rsNum));
                                if(rsNum<rsa_map.size()){
                                    slct_idx.push_back(startptr);
                                    rsNum=rsa_map.size();
                                    transNum++;
                                    double zsxz_tmp=snpinfo[startptr].beta/snpinfo[startptr].se;
                                    double pxz_tmp=pchisq(zsxz_tmp*zsxz_tmp, 1);
                                    if(pxz_tmp<transThres) // trans region extends
                                    {
                                        extended=true;
                                        transbp=snpinfo[startptr].bp;
                                    }
                                    l=startptr;
                                }
                            }
                        }
                        startptr++;
                    }
                    if(merged)
                    {
                        upperBp[upperBp.size()-1]=transuperBounder;
                        nsnp[nsnp.size()-1]=transNum;
                        extend[nsnp.size()-1]=extended;
                        merge[nsnp.size()-1]=merged;
                    }
                    else
                    {
                        tran_chr.push_back(curChr);
                        upperBp.push_back(transuperBounder);
                        lowerBp.push_back(translowerBounder);
                        nsnp.push_back(transNum);
                        extend.push_back(extended);
                        merge.push_back(merged);
                    }
                    
                }
                else if(!extract_cis_only && pxz<restThres)
                {
                    string chckstr=string(snpinfo[l].snprs);
                    //+":"+string(snpinfo[l].a1)+":"+string(snpinfo[l].a2);
                    rsa_map.insert(pair<string,int>(chckstr,rsNum));
                    if(rsNum<rsa_map.size()){
                        slct_idx.push_back(l);
                        rsNum=rsa_map.size();
                        other_idx.insert(pair<string,int>(chckstr,l));
                    }
                }
            } else {
                // I don't think this can happen. because when reading NA value is not saved.
                //printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",snpinfo[l].snprs);
            }
        }
        long cissnpnum=cis_idx.size();
        long transnum=nsnp.size();
        long transsnpnum=0;
        long othersnpnum=other_idx.size();
        for(int l=0;l<nsnp.size();l++) transsnpnum+=nsnp[l];
        if(slct_idx.size()-cissnpnum-transsnpnum-othersnpnum != 0)
        {
            printf("Something is wrong with this selection methold. Please report this bug.\n");
            exit(EXIT_FAILURE);
        }
       // printf("%ld SNPs in the cis-region, %ld SNPs in total %ld trans-region(s), %ld other SNPs have been extracted.\n",cissnpnum,transsnpnum,transnum,othersnpnum);
        //log
       
        if(cis_idx.size() || nsnp.size() || other_idx.size())
        {
            string logstr="{"+probid+","+atos(probchr)+","+atos(probbp)+"}\t";
            if(cis_idx.size()>0)
            {                
                logstr+= "["+atos(probchr)+","+ atos(snpinfo[cis_idx[0]].bp)+","+atos(snpinfo[cis_idx[cis_idx.size()-1]].bp)+","+atos(cis_idx.size())+"]\t";
            }else logstr+="[]\t";
            if(nsnp.size()>0)
            {
                for(int h=0;h<nsnp.size();h++)
                {
                    logstr+= "<"+atos(tran_chr[h])+","+ atos(lowerBp[h])+","+atos(upperBp[h])+","+atos(nsnp[h])+">\t";
                }
            }else logstr+="<>\t";
            
            if(other_idx.size()>0)
            {
                //for(iter=other_idx.begin();iter!=other_idx.end();iter++) logstr+="("+atos(snpinfo[iter->second].bp)+"), ";
                //logstr+="\n";
                logstr+="("+atos(other_idx.size())+")\n";
            }
            else logstr+="()\n";
            
            fputs(logstr.c_str(),logfile);
            fflush(logfile);
        } else {
            string logstr="{"+probid+","+atos(probchr)+","+atos(probbp)+"}\t";
            logstr+="[]\t";
            logstr+="<>\t";
            logstr+="()\n";
            fputs(logstr.c_str(),logfile);
            fflush(logfile);
        }
    }
    void save_slct_txts_sbesd(char* outFileName, long esiNum, long epiNum,vector<int> &epi2esd,probeinfolst* prbiflst, int fformat,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2,int cis_itvl, int trans_itvl, float transThres, float restThres, int addn)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        // get esd info
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        uint32_t filetype=SPARSE_FILE_TYPE_3;
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        if(addn!=-9)
            printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
        ten_ints[1]=addn;
        ten_ints[2]=(int)esiNum;
        ten_ints[3]=(int)epiNum;
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
        
        vector<uint64_t> cols((epiNum<<1)+1);;
        vector<uint32_t> rowids;
        vector<float> val;
        cols[0]=0;
        
        map<string, int>::iterator iter;
        
        
        // log file
        FILE* logfile=NULL;
        string logfname = string(outFileName)+".summary";
        logfile=fopen(logfname.c_str(), "w");
        if (!(logfile)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        string logstr="cis-window:\t"+itos(cis_itvl)+"Kb\ntrans-window:\t"+itos(trans_itvl)+"Kb\np-value threshold of trans:\t"+dtos(transThres)+"\np-value threshold of others:\t"+dtos(restThres)+"\n";
        logstr+="\ncis region is indicated by [Chr, Start bp, End bp, nsnp];\ntrans region is indicated by <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is indicated by (NumSNPs beyond cis and trans).\n";
        
        logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
        fputs(logstr.c_str(),logfile);
        fflush(logfile);
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        vector<snpinfolst> snpinfo;
        double disp=0;
        for(int j=0;j<epiNum;j++)
        {
           progress(j, disp, (int)epiNum);
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            int esdnumcurprb=epi2esd[j+1]-epi2esd[j];
            printf("\nReading %d text file(s) of probe %s.\n", esdnumcurprb, prbiflst[epi2esd[j]].probeId);
            snpinfo.clear();
            for(int k=epi2esd[j];k<epi2esd[j+1];k++)
            {
                string zname=prbiflst[k].esdpath;
                bool gzflag=has_suffix(zname, "gz");
                switch (fformat)
                {
                    case 0:
                        read_smr_sa(snpinfo, prbiflst[k].esdpath);
                        break;
                    case 1:
                        read_plink_qassoc(snpinfo, prbiflst[k].esdpath);
                        break;
                    case 2:
                        read_gemma_qassoc(snpinfo, prbiflst[k].esdpath);
                        break;
                    case 3:
                        read_bolt_qassoc(snpinfo, prbiflst[k].esdpath);
                        break;
                }
            }
            snpinfolst* sortptr=&snpinfo[0];
            qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
            
            probeinfolst prbifo=prbiflst[epi2esd[j]];
            vector<int> slct_idx;
            slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl, trans_itvl, transThres, restThres,logfile,false); //slct_idx with no order if there are trans-rgeions
            stable_sort(slct_idx.begin(),slct_idx.end());
            vector<string> _rs(slct_idx.size()), _a1(slct_idx.size()),_a2(slct_idx.size());
            vector<float> _beta(slct_idx.size()), _se(slct_idx.size());
            
            for(int l=0;l<slct_idx.size();l++) {
                _rs[l]=snpinfo[slct_idx[l]].snprs;
                if(fformat!=1) _a1[l]=snpinfo[slct_idx[l]].a1;
                if(fformat!=1) _a2[l]=snpinfo[slct_idx[l]].a2;
                _beta[l]=snpinfo[slct_idx[l]].beta;
                _se[l]=snpinfo[slct_idx[l]].se;
            }
            
            //when using --geno-uni, rsid is not the same as slct_idx, if there are NA values in the raw files.
            vector<int> rsid(_rs.size());
            for (int l = 0; l<_rs.size(); l++){
                iter = esi_map.find(_rs[l]);
                if (iter != esi_map.end()) rsid[l]=iter->second;
                else {
                    printf("ERROR: SNP %s is not in the output SNP set. If you are using --geno-uni, please disable it then try again. Otherwise please report this bug.\n",_rs[l].c_str());
                    exit(EXIT_FAILURE);
                }
            }
          //  map<string, int> rsa_map;
          //  long rsNum=0;
            for(int l=0;l<rsid.size();l++)
            {
                //if(fabs(_se[l]+9)>1e-6) // can move this. the NA is controled in slct_sparse_per_prb
                //{
                  //  string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                  //  rsa_map.insert(pair<string,int>(chckstr,l)); // in slct_sparse_per_prb, ras_map can privent selecting duplicate SNPs and double-slelecting SNPs. so we can move rsa_map here.
                   // if(rsNum<rsa_map.size())
                   // {
                        if(fformat!=1) //not for plink format
                        {
                            if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l] ){
                                val.push_back(_beta[l]);
                                rowids.push_back(rsid[l]);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(rsid[l]);
                            } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l] ){
                                val.push_back(-1.0*_beta[l]);
                                rowids.push_back(rsid[l]);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(rsid[l]);
                            } else {
                                printf("ERROR: inconsistent allele pairs of SNP %s found.\n",_rs[l].c_str());
                                printf("Discrepant Allele pairs: (%s,%s) with (%s,%s).\n",esi_a1[rsid[l]].c_str(), esi_a2[rsid[l]].c_str(),_a1[l].c_str(), _a2[l].c_str());
                                exit(EXIT_FAILURE);
                                
                                //this part is for multi-allelic SNPs. since we don't save multi-allelic SNPs anymore, so we should disable it.
                                /*
                                int did=-9;
                                float sig=1.0;
                                for(int m=0;m<esi_rs.size();m++)
                                {
                                    if(esi_rs[m]==_rs[l])
                                    {
                                        if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                        {
                                            did=m;
                                            break;
                                        }
                                        if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                        {
                                            did=m;
                                            sig=-1.0;
                                            break;
                                        }
                                    }
                                }
                                if(did==-9)
                                {
                                    printf("ERROR: This would not go to happen. Please report this bug.");
                                    exit(EXIT_FAILURE);
                                }
                                val.push_back(sig*_beta[l]);
                                rowids.push_back(did);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(did);
                                 */
                            }
                        } else {
                            val.push_back(_beta[l]);
                            rowids.push_back(rsid[l]);
                            tmpse.push_back(_se[l]);
                            tmprid.push_back(rsid[l]);
                        }
                 //       rsNum=rsa_map.size();
                //    } else {
                 //       printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s. \n",_rs[l].c_str(),prbiflst[epi2esd[j]].probeId.c_str());
                 //   }
                //} else {
                //    printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",_rs[l].c_str());
               // }
            }
            for(int k=0;k<tmpse.size();k++)
            {
                val.push_back(tmpse[k]);
                rowids.push_back(tmprid[k]);
            }
            uint64_t real_num=tmpse.size();
            cols[(j<<1)+1]=real_num+cols[j<<1];
            cols[j+1<<1]=(real_num<<1)+cols[j<<1];
            
            for(int k=0;k<snpinfo.size();k++)
            {
                if(fformat!=1) delete(snpinfo[k].a1);
                if(fformat!=1) delete(snpinfo[k].a2);
                delete(snpinfo[k].snprs);
            }

        }
        uint64_t valNum=val.size();
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
        fwrite (&val[0],sizeof(float), val.size(), smr1);
        fclose (smr1);
        
        printf("Summary data of the specified SNPs and probes has been saved in %s.\n", logfname.c_str());
        cout<<"\nEffect sizes (beta) and SE for "<<epiNum<<" Probes have been saved in a binary file [" + esdfile + "]." <<endl;
        fclose(logfile);
        
    }
    void make_besd(char*outFileName, char* syllabusName, bool gctaflag,bool plinkflag,bool gemmaflag,bool merlinflag,bool boltflag, bool save_dense_flag,int cis_itvl, int trans_itvl, float transThres, float restThres,bool samegeno, int addn)
    {
        if(samegeno) {
            printf("WARNING: --geno-uni is enabled. Please ensure the SNPs and their alleles identical across all the text files.\n");
        }
        if(syllabusName==NULL) throw("Error: please input eQTL file list by the flag --eqtl-flist.");
        int flagcount=0, fformat=0;
        if(gctaflag){ flagcount++; fformat=0;}
        if(plinkflag){ flagcount++;fformat=1;}
        if(gemmaflag){ flagcount++; fformat=2;}
        if(boltflag){ flagcount++; fformat=3;}
        if(merlinflag){ flagcount++; fformat=4; printf("WARNING: THIS OPTION IS ONLY FOR PERSONAL USE. ONLY with --make-besd-dense.\n");}
        if(flagcount==0 || flagcount>1) throw("ERROR: please verify the file format flags. Only one flag can be specified.\n");
        vector<probeinfolst> prbiflstinfo;
        int esdfileNum=read_probeinfolst(prbiflstinfo,syllabusName);
        probeinfolst* prbiflst=&prbiflstinfo[0];
        qsort(prbiflst,esdfileNum,sizeof(probeinfolst),comp);
        
        //get epi
        printf("Generating the .epi file...\n");
        map<string, int> epi_map;
        long epiNum=0;
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi)
        {
            printf("ERROR: can not open the EPI file %s to save!\n", epifile.c_str());
            exit(EXIT_FAILURE);
        }
        vector<int> epi2esd;
        for (int j = 0;j <esdfileNum; j++) {
            epi_map.insert(pair<string,int>(prbiflst[j].probeId,epiNum));
            if(epiNum<epi_map.size())
            {
                epi2esd.push_back(j);
                epi<<prbiflst[j].probechr<<'\t'<<prbiflst[j].probeId<<'\t'<<0<<'\t'<<prbiflst[+j].bp<<'\t'<<prbiflst[j].genename<<'\t'<<prbiflst[j].orien<<'\n';
                epiNum++;
            }
        }
        epi2esd.push_back((int)esdfileNum);
        epi.close();
        printf("%ld probes have been saved in the file \"%s\".\n",epiNum, epifile.c_str());
        
        //get esi: from bim file
        vector<snpinfolst> snpinfo;
        printf("\nScaning %d ESD files to generate ESI file...\n", esdfileNum);
        int filenum2scan=0;
        if(samegeno) filenum2scan=1;
        else filenum2scan=esdfileNum;
        uint64_t ttlv = get_esi_info(snpinfo,&prbiflst,filenum2scan,fformat);
        if(samegeno) ttlv*=esdfileNum;
        long esiNum=snpinfo.size();
        snpinfolst* esiptr=&snpinfo[0];
        qsort(esiptr,esiNum,sizeof(snpinfolst),comp_esi);
        vector<string> esi_rs(esiNum);
        vector<string> esi_a1(esiNum);
        vector<string> esi_a2(esiNum);
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("ERROR: can not open the ESI file to save!");
        for (int j = 0;j <esiNum; j++) {
            esi<<snpinfo[j].snpchr<<'\t'<<snpinfo[j].snprs<<'\t'<<snpinfo[j].gd<<'\t'<<snpinfo[j].bp<<'\t'<<snpinfo[j].a1<<'\t'<<snpinfo[j].a2<<'\t'<<(fabs(snpinfo[j].freq+9)>1e-6?atos(snpinfo[j].freq):"NA")<<'\n';
            esi_rs[j]=snpinfo[j].snprs;
            esi_a1[j]=snpinfo[j].a1;
            esi_a2[j]=snpinfo[j].a2;
        }
        esi.close();
        cout<<esiNum<<" SNPs have been saved in the file [" + esifile + "]."<<endl;
        printf("\nScaning %d ESD files to generate BESD file...\n", esdfileNum);
        if(save_dense_flag)
        {
            double sparsity=1.0*ttlv/(esiNum*epiNum);
            if(sparsity>=0.4)
            {
                //printf("The density of your data is %f. The data will be saved in dense format.\n", sparsity);
                save_txts_dbesd(outFileName, esiNum, epiNum,epi2esd, prbiflst,fformat, esi_rs, esi_a1,esi_a2,addn);
            } else {
                //printf("The density of your data is %f. The data will be saved in sparse format.\n", sparsity);
                save_full_txts_sbesd( outFileName,  esiNum,  epiNum,epi2esd,prbiflst,fformat,esi_rs,esi_a1,esi_a2,addn);
            }
            
        } else {
            save_slct_txts_sbesd(outFileName, esiNum, epiNum,epi2esd,prbiflst,fformat,esi_rs,esi_a1,esi_a2, cis_itvl,  trans_itvl,  transThres,  restThres,addn);
        }
        free_probelist(prbiflstinfo);
        free_snplist(snpinfo);
        
    }
    void make_besd_fmat(char* fmatfileName, char* outFileName,bool mateqtlflag, bool fastnflag, bool fastpflag, bool qtltoolsnflag, bool qtltoolspflag,int addn)
    {
        int tcount=0;
        if(mateqtlflag) tcount++;
        if(fastnflag) tcount++;
        if(fastpflag) tcount++;
        if(qtltoolsnflag)  tcount++;
        if(qtltoolspflag)  tcount++;
        if(tcount!=1) {
            printf("ERROR: Please specify only one file type.\n ");
            exit(EXIT_FAILURE);
        }
        int rspos=-9, prbpos=-9, betapos=-9, ppos=-9, cnum=-9,prbchrpos=-9, prbbppos=-9, prbstrandpos=-9, rschrpos=-9, rsbppos=-9;
        if(mateqtlflag) {
            cnum=5; // general 6; but with noFDRsaveMemory = TRUE, the output column is 5.
            rspos=0;
            prbpos=1;
            betapos=2;
            ppos=4;
        }
        else if(fastnflag) {
            cnum=5;
            rspos=1;
            prbpos=0;
            betapos=4;
            ppos=3;
        }
        else if(fastpflag) {
            cnum=11;
            rspos=5;
            prbpos=0;
            betapos=8;
            ppos=7;
        }
        else if(qtltoolsnflag) {
            cnum=14;
            rspos=7;
            prbpos=0;
            betapos=12;
            ppos=11;
            prbchrpos=1;
            prbbppos=2;
            prbstrandpos=4;
            rschrpos=8;
            rsbppos=9;
 
        }
        else if(qtltoolspflag) {
            cnum=19;
            rspos=7;
            prbpos=0;
            betapos=16;
            ppos=15;
            prbchrpos=1;
            prbbppos=2;
            prbstrandpos=4;
            rschrpos=8;
            rsbppos=9;
        }
        gzFile gzfile=NULL;
        FILE* qfile=NULL;
        bool gzflag=has_suffix(fmatfileName, "gz");
        if(gzflag)
        {
            gzfile = gzopen(fmatfileName, "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         fmatfileName, strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            qfile=fopen(fmatfileName, "r");
            if (!qfile) {
                printf("ERROR: Can't open file %s.\n ", fmatfileName);
                exit(EXIT_FAILURE);
            }
        }
        
        printf("Reading eQTL summary data from %s ...\n", fmatfileName);
        vector<string> prbs;
        vector<string> snps;
        vector<int> prbchr;
        vector<int> prbbp;
        vector<string> prbstrand;
        vector<int> rschr;
        vector<int> rsbp;
        vector< vector<uint32_t> > _ttl_rsid;
        vector< vector<float> > _ttl_beta;
        vector< vector<float> > _ttl_se;
        
        
        char buf[MAX_LINE_SIZE];
        vector<string> vs_buf;
        if(mateqtlflag)
        {
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
            else fgets(buf, MAX_LINE_SIZE, qfile); //header
            if(buf[0]=='\0')
            {
                printf("ERROR: the first row of the file %s is empty.\n",fmatfileName);
                exit(EXIT_FAILURE);
            }
            split_string(buf, vs_buf, ", \t\n");
            to_upper(vs_buf[0]);
            if(vs_buf[0]!="SNP") {
                printf("ERROR: the headers of query file should start with \"SNP\".\n");
                exit(EXIT_FAILURE);
            }
            if(vs_buf.size() < 5 || vs_buf.size() > 7) {
                printf("ERROR: the input file is not in Matrix eQTL output format.\n");
                exit(EXIT_FAILURE);
            }
            if(vs_buf.size() == 5) {
                cnum=5;
                printf("Total 5 columns in the Matrix eQTL output.\n");
            } else {
                cnum=6;
            }
        }
        
        long lineNum(0);
        long esiNum(0);
        long epiNum(0);
        map<string, int> esi_map;
        map<string,int>::iterator iter;
        map<string, int> epi_map;
        map<string, int> rs_prb_map;
        long rsprbmapsize(0);
        while((qfile!=NULL && !feof(qfile)) ||(gzfile!=NULL && !gzeof(gzfile)))
        {
            buf[0]='\0';
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else fgets(buf, MAX_LINE_SIZE, qfile);
            if(buf[0]!='\0')
            {
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=cnum) {
                    printf("ERROR: the number of columns is incorrect in row %ld.\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if( (vs_buf[rspos]=="NA" || vs_buf[rspos]=="na"))
                {
                    printf("ERROR: the SNP name is \"NA\".\n");
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[prbpos]=="NA" || vs_buf[prbpos]=="na")
                {
                    printf("ERROR: the probe name is \"NA\".\n");
                    exit(EXIT_FAILURE);
                }
                
                if(vs_buf[betapos]=="NA" || vs_buf[betapos]=="na") {
                    
                    //printf("WARNING: this row is omitted because the effect size of the SNP is missing (\"NA\").\n");
                    //printf("%s\n",buf);
                    continue;
                }
                if(vs_buf[ppos]=="NA" || vs_buf[ppos]=="na") {
                    
                    //printf("WARNING: the p-value of the SNP is missing (\"NA\"), this row is omitted.\n");
                    //printf("%s\n",buf);
                    continue;
                }
                if(atof(vs_buf[ppos].c_str())<0) {
                    printf("ERROR: p-value should be positive in row %ld.\n", lineNum+2);
                    printf("%s\n",buf);
                    exit(EXIT_FAILURE);
                }
                string rsprb=vs_buf[rspos]+":"+vs_buf[prbpos];
                rs_prb_map.insert(pair<string, int>(rsprb.c_str(), lineNum));
                if(rsprbmapsize<rs_prb_map.size())
                {
                    rsprbmapsize=rs_prb_map.size();
                } else {
                    printf("ERROR: duplicated records found for the SNP %s and the probe %s.\n", vs_buf[0].c_str(),vs_buf[1].c_str());
                    exit(EXIT_FAILURE);
                }
                int rsid=-9;
                iter=esi_map.find(vs_buf[rspos]);
                if(iter==esi_map.end()){
                    esi_map.insert(pair<string, int>(vs_buf[rspos].c_str(), esiNum));
                    snps.push_back(vs_buf[rspos].c_str());
                    if(rschrpos>=0)
                    {
                        string chrstrtmp=vs_buf[rschrpos];
                        if(has_prefix(vs_buf[rschrpos], "chr") || has_prefix(vs_buf[rschrpos], "CHR")) chrstrtmp.erase(0,3);
                        rschr.push_back(atoi(chrstrtmp.c_str()));
                    } else rschr.push_back(-9);
                    if(rsbppos>=0 && (qtltoolsnflag || qtltoolspflag)) rsbp.push_back((atoi(vs_buf[rsbppos].c_str())+atoi(vs_buf[rsbppos+1].c_str()))/2);
                    else rsbp.push_back(-9);
                    rsid=(int)esiNum++;
                }
                else
                {
                    rsid=iter->second;
                }
                iter=epi_map.find(vs_buf[prbpos]);
                if(iter!=epi_map.end()){
                    int idx=iter->second;
                    double betatmp=atof(vs_buf[betapos].c_str());
                    double ptmp=atof(vs_buf[ppos].c_str());
                    _ttl_rsid[idx].push_back(rsid);
                    _ttl_beta[idx].push_back(betatmp);
                    _ttl_se[idx].push_back(adjSE(betatmp, ptmp));
                    
                } else {
                    
                    epi_map.insert(pair<string,int>(vs_buf[prbpos],epiNum));
                    prbs.push_back(vs_buf[prbpos].c_str());
                    if(prbchrpos>=0)
                    {
                        string chrstrtmp=vs_buf[prbchrpos];
                        if(has_prefix(vs_buf[prbchrpos], "chr") || has_prefix(vs_buf[prbchrpos], "CHR")) chrstrtmp.erase(0,3);
                        prbchr.push_back(atoi(chrstrtmp.c_str()));
                    } else prbchr.push_back(-9);
                    if(prbbppos>=0 && (qtltoolsnflag || qtltoolspflag)) prbbp.push_back((atoi(vs_buf[prbbppos].c_str())+atoi(vs_buf[prbbppos+1].c_str()))/2);
                    else prbbp.push_back(-9);
                    if(prbstrandpos>=0) prbstrand.push_back(vs_buf[prbstrandpos]);
                    else prbstrand.push_back("NA");
                    vector<uint32_t> rsidtmp;
                    rsidtmp.push_back(rsid);
                    _ttl_rsid.push_back(rsidtmp);
                    double betatmp=atof(vs_buf[betapos].c_str());
                    double ptmp=atof(vs_buf[ppos].c_str());
                    vector<float> btmp;
                    btmp.push_back(betatmp);
                    _ttl_beta.push_back(btmp);
                    vector<float> setmp;
                    setmp.push_back(adjSE(betatmp, ptmp));
                    _ttl_se.push_back(setmp);
                    epiNum++;
                }
                lineNum++;
            }
        }
        printf("%ld rows to be included from %s.\n", lineNum, fmatfileName);
        if(gzflag) gzclose(gzfile);
        else fclose(qfile);
        
        printf("\nGenerating the .epi file...\n");
        
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the EPI file " + epifile + " to save!");
        for (int j = 0;j <epiNum; j++) {
            epi<<(prbchr[j]>0?atos(prbchr[j]):"NA")<<'\t'<<prbs[j]<<'\t'<<0<<'\t'<<(prbchr[j]>0?atos(prbbp[j]):"NA")<<'\t'<<"NA"<<'\t'<<prbstrand[j]<<'\n';
        }
        epi.close();
        printf("%ld probes have been saved in the file %s.\n",epiNum,epifile.c_str());
        
        
        printf("\nGenerating the .esi file...\n");
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the ESI file to save!");
        esi_map.clear();
        epi_map.clear();
        for (int j = 0;j <esiNum; j++) {
            esi<<(rschr[j]>0?atos(rschr[j]):"NA")<<'\t'<<snps[j]<<'\t'<<0<<'\t'<<(rsbp[j]>0?atos(rsbp[j]):"NA")<<'\t'<<"NA"<<'\t'<<"NA"<<'\t'<<"NA"<<'\n';
            esi_map.insert(pair<string,int>(snps[j],j));
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",esiNum,esifile.c_str());
        
        double sparsity=1.0*lineNum/(esiNum*epiNum);
        printf("\nGenerating the .besd file...\n");
        
        // quanlity control and correctness are done above.
            if(sparsity>0.4)
            {
                //printf("The density of your data is %f. The data will be saved in dense format.\n", sparsity);
                
                string esdfile=string(outFileName)+string(".besd");
                FILE * smr1;
                smr1 = fopen (esdfile.c_str(), "wb");
                if (!(smr1)) {
                    printf("ERROR: failed to open file %s.\n",esdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                uint32_t filetype=DENSE_FILE_TYPE_3;
                vector<int> ten_ints(RESERVEDUNITS);
                ten_ints[0]=filetype;
                if(addn!=-9)
                    printf("Saving the sample size %d to the file %s.\n", addn, esdfile.c_str());
                ten_ints[1]=addn;
                ten_ints[2]=(int)esiNum;
                ten_ints[3]=(int)epiNum;
                for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
                
                uint64_t bsize=(uint64_t)esiNum<<1;
                float* buffer=(float*)malloc (sizeof(float)*bsize);
                if (NULL == buffer) {
                    fprintf(stderr, "Malloc failed\n");
                    exit(-1);
                }
                double disp=0;
                for(int j=0;j<epiNum;j++)
                {
                    progress(j, disp, (int)epiNum);
                  
                    for(int k=0;k<bsize;k++) buffer[k]=-9; //init
                    for(int l=0;l<_ttl_rsid[j].size();l++)
                    {
                        buffer[_ttl_rsid[j][l]] = _ttl_beta[j][l];
                        buffer[esiNum + _ttl_rsid[j][l]] = _ttl_se[j][l];
                    }
                    fwrite (buffer,sizeof(float), bsize, smr1);
                }
                fclose (smr1);
                free(buffer);
                cout<<"Effect sizes (beta) and SE for "<<epiNum<<" probes and "<<esiNum<<" SNPs have been saved in a binary file [" + esdfile + "]." <<endl;
                
            }
            else
            {
                //printf("The density of your data is %f. The data will be saved in sparse format.\n", sparsity);
                
                string esdfile=string(outFileName)+string(".besd");
                FILE * smr1;
                smr1 = fopen (esdfile.c_str(), "wb");
                if (!(smr1)) {
                    printf("ERROR: failed to open file %s.\n",esdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                uint32_t filetype=SPARSE_FILE_TYPE_3;
                vector<int> ten_ints(RESERVEDUNITS);
                ten_ints[0]=filetype;
                if(addn!=-9)
                    printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
                ten_ints[1]=addn;
                ten_ints[2]=(int)esiNum;
                ten_ints[3]=(int)epiNum;
                for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
                
                vector<uint64_t> cols((epiNum<<1)+1);;
                uint64_t valNum=0;
                cols[0]=0;
                
                map<string, int>::iterator iter;
                for(int j=0;j<epiNum;j++)
                {
                    uint64_t real_num = _ttl_beta[j].size();
                    cols[(j<<1)+1]=real_num+cols[j<<1];
                    cols[j+1<<1]=(real_num<<1)+cols[j<<1];
                    valNum+=real_num*2;
                }
                fwrite (&valNum,sizeof(uint64_t), 1, smr1);
                fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
                double disp=0;
                for(int j=0;j<epiNum;j++)
                {
                    progress(j, disp, (int)epiNum);
               
                    fwrite (&_ttl_rsid[j][0],sizeof(uint32_t), _ttl_rsid[j].size(), smr1);
                    fwrite (&_ttl_rsid[j][0],sizeof(uint32_t), _ttl_rsid[j].size(), smr1);
                }
                disp=0;
                for(int j=0;j<epiNum;j++)
                {
                    progress(j, disp, (int)epiNum);
                    fwrite(&_ttl_beta[j][0], sizeof(float), _ttl_beta[j].size(), smr1);
                    fwrite(&_ttl_se[j][0], sizeof(float), _ttl_se[j].size(), smr1);
                }
                fclose (smr1);
                cout<<"Effect sizes (beta) and SE for "<<epiNum<<" probes and "<<esiNum<<" SNPs have been saved in a binary file [" + esdfile + "]." <<endl;
            }
        
        
    }
    void make_besd_byQfile(char* qfileName, char* outFileName,bool save_dense_flag,int cis_itvl, int trans_itvl, float transThres, float restThres,int addn)
    {
        FILE* qfile=fopen(qfileName, "r");
        if (!qfile) {
            printf("ERROR: Can't open file %s.\n ", qfileName);
            exit(EXIT_FAILURE);
        }
        
        printf("Reading eQTL summary data from %s ...\n", qfileName);
        vector<probeinfolst> prbiflst;
        vector<snpinfolst> snpiflst;
        vector< vector<string> > _ttl_rs;
        vector< vector<float> > _ttl_beta;
        vector< vector<float> > _ttl_se;
        
        vector< vector<string> > _ttl_a1;
        vector< vector<string> > _ttl_a2;
        vector< vector<int> > _ttl_chr;
        vector< vector<int> > _ttl_bp;
        
        bool warningdone=false;
        bool warningdone1=false;
        char buf[MAX_LINE_SIZE];
        fgets(buf, MAX_LINE_SIZE, qfile); //header
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",qfileName);
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: the headers of query file should start with \"SNP\".\n");
            exit(EXIT_FAILURE);
        }
        if(vs_buf.size()==13) {
            printf("ERROR: the input file could be in the old Query format. Please update it by adding one column of the effect allele frequency.\n");
            printf("http://cnsgenomics.com/software/smr/query.html\n");
            exit(EXIT_FAILURE);
        }
        long lineNum(0);
        long esiNum(0);
        long epiNum(0);
        map<string, int> esi_map;
        map<string,int>::iterator iter;
        map<string, int> epi_map;
        map<string, int> rs_prb_map;
        long rsprbmapsize(0);
        while(fgets(buf, MAX_LINE_SIZE, qfile))
        {
            vs_buf.clear();
            int col_num = split_string(buf, vs_buf, ", \t\n");
            if(col_num!=14) {
                printf("ERROR: the number of columns is incorrect in row %ld.\n", lineNum+2);
                exit(EXIT_FAILURE);
            }
            if( (vs_buf[0]=="NA" || vs_buf[0]=="na"))
            {
                printf("ERROR: SNP name is \"NA\".\n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[1]=="NA" || vs_buf[1]=="na"))
            {
                printf("ERROR: chromosome infomation of SNP is \"NA\".\n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[6]=="NA" || vs_buf[6]=="na"))
            {
                printf("ERROR: Probe name is missing.\n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[7]=="NA" || vs_buf[7]=="na"))
            {
                printf("ERROR: chromosome infomation of probe is \"NA\". \n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[2]=="NA" || vs_buf[2]=="na"))
            {
                printf("ERROR: SNP position is \"NA\". \n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[8]=="NA" || vs_buf[8]=="na"))
            {
                printf("ERROR: probe position is \"NA\". \n");
                exit(EXIT_FAILURE);
            }
            if((vs_buf[3]=="NA" || vs_buf[4]=="NA" || vs_buf[3]=="na" || vs_buf[4]=="na") && !warningdone)
            {
               printf("WARING: Allele information of one or more SNPs is \"NA\". This may cause incorrect positive/negative sign of effect size when conducting allele check.\n");
                warningdone=true;
               
            }
            if((vs_buf[5]=="NA" || vs_buf[5]=="na" ) )
            {
                if(!warningdone1)
                {
                    printf("WARING: allele frequency of one or more SNPs are missing.\n");
                    warningdone1=true;
                }
            } else {
                float fq=atof(vs_buf[5].c_str());
                if(fq<=1e-8 || fq>=1)
                {
                    printf("ERROR: Allele frequency should be between 0 and 1 in row %ld.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
            }
            int tmpchr;
            to_upper(vs_buf[1]);
            //if(vs_buf[1]=="X" ) tmpchr=23;
            //else if(vs_buf[1]=="Y") tmpchr=24;
            //else tmpchr=atoi(vs_buf[1].c_str());
            if(!vs_buf[1].compare("X")) tmpchr=23;
            else if(!vs_buf[1].compare("Y")) tmpchr=24;
            else tmpchr=atoi(vs_buf[1].c_str());
            to_upper(vs_buf[3]);
            to_upper(vs_buf[4]);
            to_upper(vs_buf[5]);
            to_upper(vs_buf[11]);
            to_upper(vs_buf[12]);
            to_upper(vs_buf[13]);
            if(!vs_buf[11].compare("NA")) {
                
                printf("WARNING: this row is omitted because the effect size of the SNP is missing (\"NA\").\n");
                printf("%s\n",buf);
                continue;
            }
            if(!vs_buf[12].compare("NA") && vs_buf[13]=="NA") {
                
                printf("WARNING: both of SE and P are missing (\"NA\"), this row is omitted.\n");
                printf("%s\n",buf);
                continue;
            }
            if(vs_buf[13]!="NA" && atof(vs_buf[13].c_str())<0) {
                printf("ERROR: p-value should be positive in row %ld.\n", lineNum+2);
                printf("%s\n",buf);
                exit(EXIT_FAILURE);
            }
            string rsprb=vs_buf[0]+":"+vs_buf[6];
            rs_prb_map.insert(pair<string, int>(rsprb.c_str(), lineNum));
            if(rsprbmapsize<rs_prb_map.size())
            {
                rsprbmapsize=rs_prb_map.size();
            } else {
                printf("ERROR: duplicated records found for the SNP %s and the probe %s.\n", vs_buf[0].c_str(),vs_buf[6].c_str());
                exit(EXIT_FAILURE);
            }
            iter=esi_map.find(vs_buf[0]);
            if(iter!=esi_map.end()){
                int idx=iter->second;
                string esi_a1=string(snpiflst[idx].a1);
                string esi_a2=string(snpiflst[idx].a2);
                int esi_bp=snpiflst[idx].bp;
                int esi_chr=snpiflst[idx].snpchr;
                if(tmpchr!=esi_chr) {
                    printf("ERROR: the same SNP ID %s on different chromosomes, please check.\n", vs_buf[0].c_str());
                    exit(EXIT_FAILURE);
                }
                if(atoi(vs_buf[2].c_str())!=esi_bp)
                {
                    printf("ERROR: multiple positions found for the same SNP %s on Chromosome %s, please check.\n", vs_buf[0].c_str(), vs_buf[1].c_str());
                    exit(EXIT_FAILURE);
                }
                if(esi_a1==vs_buf[3] && esi_a2==vs_buf[4]) {
                    
                } else if(esi_a1==vs_buf[4] && esi_a2==vs_buf[3]) {
                    printf("WARING: switched ref allele with alt allele of SNP %s found. Then the sign of beta value is changed.\n", vs_buf[0].c_str());
                    vs_buf[11]=atos(-1.0*atof(vs_buf[11].c_str()));
                    string a0=vs_buf[3];
                    vs_buf[3]=vs_buf[4];
                    vs_buf[4]=a0;
                } else {
                    printf("ERROR: inconsistent alleles for the SNP %s with multiple alleles <%s,%s> and <%s,%s>.\n", vs_buf[0].c_str(),esi_a1.c_str(),esi_a2.c_str(),vs_buf[3].c_str(),vs_buf[4].c_str());
                    exit(EXIT_FAILURE);
                }
                
            } else {
                esi_map.insert(pair<string, int>(vs_buf[0].c_str(), esiNum));
                snpinfolst snptmp;
                snptmp.snpchr=tmpchr;
                strcpy2(&snptmp.snprs, vs_buf[0]);
                snptmp.bp=atoi(vs_buf[2].c_str());
                strcpy2(&snptmp.a1,vs_buf[3]);
                strcpy2(&snptmp.a2,vs_buf[4]);
                snptmp.gd=0;
                snptmp.freq=vs_buf[5]=="NA"?-9:atof(vs_buf[5].c_str());
                snpiflst.push_back(snptmp);
                esiNum++;
            }
            int prbchr=0;
            if(vs_buf[7]=="X" || vs_buf[7]=="x") prbchr=23;
            else if(vs_buf[7]=="Y" || vs_buf[7]=="y") prbchr=24;
            else prbchr=atoi(vs_buf[7].c_str());
            iter=epi_map.find(vs_buf[6]);
            if(iter!=epi_map.end()){
                int idx=iter->second;
                int epi_bp=prbiflst[idx].bp;
                int epi_chr=prbiflst[idx].probechr;
                if(epi_chr != prbchr)
                {
                    printf("ERROR: the same probe ID %s found on different chromosomes.\n",vs_buf[5].c_str() );
                    exit(EXIT_FAILURE);
                }
                if(epi_bp != atoi(vs_buf[8].c_str()))
                {
                    printf("ERROR: multiple posisions found for the same probe %s.\n",vs_buf[5].c_str() );
                    exit(EXIT_FAILURE);
                }
                _ttl_rs[idx].push_back(vs_buf[0].c_str());
                _ttl_beta[idx].push_back(atof(vs_buf[11].c_str()));
                if(vs_buf[13]=="NA") _ttl_se[idx].push_back(atof(vs_buf[12].c_str()));
                else {
                    double betatmp=atof(vs_buf[11].c_str());
                    double ptmp=atof(vs_buf[13].c_str());
                    if(ptmp<0)
                    {
                        printf("ERROR: p-value should be positive in row %ld.\n",lineNum+2);
                        printf("%s\n",buf);
                        exit(EXIT_FAILURE);
                    }
                    if(ptmp==0) {
                        ptmp=__DBL_MIN__;
                        printf("WARNING: p-value of 0 found in row %ld and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                        printf("%s\n",buf);
                    }
                    if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[12]!="NA") || ptmp==1) _ttl_se[idx].push_back(atof(vs_buf[12].c_str()));
                    else _ttl_se[idx].push_back(adjSE(betatmp, ptmp));
                }
                if(!save_dense_flag) {
                    _ttl_a1[idx].push_back(vs_buf[3].c_str());
                    _ttl_a2[idx].push_back(vs_buf[4].c_str());
                    _ttl_chr[idx].push_back(tmpchr);
                    _ttl_bp[idx].push_back(atoi(vs_buf[2].c_str()));
                }
                
            } else {
                probeinfolst tmpprbinfo;
                epi_map.insert(pair<string,int>(vs_buf[6],epiNum));
                strcpy2(&tmpprbinfo.probeId,vs_buf[6]);
                vs_buf[6].c_str();
                tmpprbinfo.probechr=prbchr;
                tmpprbinfo.bp=atoi(vs_buf[8].c_str());
                strcpy2(&tmpprbinfo.genename,vs_buf[9]);
                tmpprbinfo.orien=vs_buf[10][0];
                tmpprbinfo.gd=0;
                tmpprbinfo.bfilepath=NULL;
                tmpprbinfo.esdpath=NULL;
                prbiflst.push_back(tmpprbinfo);
                vector<string> rstmp;
                rstmp.push_back(vs_buf[0].c_str());
                _ttl_rs.push_back(rstmp);
                vector<float> betatmp;
                betatmp.push_back(atof(vs_buf[11].c_str()));
                _ttl_beta.push_back(betatmp);
                vector<float> setmp;
                if(vs_buf[13]=="NA") setmp.push_back(atof(vs_buf[12].c_str()));
                else {
                    double betatmp=atof(vs_buf[11].c_str());
                    double ptmp=atof(vs_buf[13].c_str());
                    if(ptmp<0)
                    {
                        printf("ERROR: p-value should be positive in row %ld.\n",lineNum+2);
                        printf("%s\n",buf);
                        exit(EXIT_FAILURE);
                    }
                    if(ptmp==0) {
                        ptmp=__DBL_MIN__;
                        printf("WARNING: p-value of 0 found in row %ld and changed to min double value %e.\n",lineNum+2,__DBL_MIN__);
                        printf("%s\n",buf);
                    }
                    if((ptmp<MIN_PVAL_ADJUSTED && vs_buf[12]!="NA") || ptmp==1) setmp.push_back(atof(vs_buf[12].c_str()));
                    else setmp.push_back(adjSE(betatmp, ptmp));
                }
                _ttl_se.push_back(setmp);
                if(!save_dense_flag) {
                    vector<string> a1tmp;
                    a1tmp.push_back(vs_buf[3].c_str());
                    _ttl_a1.push_back(a1tmp);
                    vector<string> a2tmp;
                    a2tmp.push_back(vs_buf[4].c_str());
                    _ttl_a2.push_back(a2tmp);
                    vector<int> chrtmp;
                    chrtmp.push_back(tmpchr);
                    _ttl_chr.push_back(chrtmp);
                    vector<int> bptmp;
                    bptmp.push_back(atoi(vs_buf[2].c_str()));
                    _ttl_bp.push_back(bptmp);
                }
                
                epiNum++;
            }
            lineNum++;
        }
        printf("%ld SNPs to be included from %s.\n", lineNum, qfileName);
        fclose(qfile);
		for (int j = 0; j < epiNum; j++) prbiflst[j].gd = j;
        printf("\nGenerating the .epi file...\n");
        probeinfolst* epiptr=&prbiflst[0];
        qsort(epiptr,epiNum,sizeof(probeinfolst),comp);
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the EPI file " + epifile + " to save!");
        for (int j = 0;j <epiNum; j++) {
            epi<<(prbiflst[j].probechr==0?"NA":atos(prbiflst[j].probechr))<<'\t'<<prbiflst[j].probeId<<'\t'<<0<<'\t'<<(prbiflst[j].bp==0?"NA":atos(prbiflst[j].bp))<<'\t'<<prbiflst[j].genename<<'\t'<<prbiflst[j].orien<<'\n';
        }
        epi.close();
        printf("%ld probes have been saved in the file %s.\n",epiNum,epifile.c_str());
        
        
        printf("\nGenerating the .esi file...\n");
        snpinfolst* esiptr=&snpiflst[0];
        qsort(esiptr,esiNum,sizeof(snpinfolst),comp_esi);
        epi_map.clear();
        esi_map.clear();
        
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the ESI file to save!");
        for (int j = 0;j <esiNum; j++) {
            esi<<(snpiflst[j].snpchr==0?"NA":atos(snpiflst[j].snpchr))<<'\t'<<snpiflst[j].snprs<<'\t'<<snpiflst[j].gd<<'\t'<<(snpiflst[j].bp==0?"NA":atos(snpiflst[j].bp))<<'\t'<<snpiflst[j].a1<<'\t'<<snpiflst[j].a2<<'\t'<<(snpiflst[j].freq+9>1e-6?atos(snpiflst[j].freq):"NA")<<'\n';
            esi_map.insert(pair<string,int>(snpiflst[j].snprs,j));
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",esiNum,esifile.c_str());
        
        double sparsity=1.0*lineNum/(esiNum*epiNum);
        printf("\nGenerating the .besd file...\n");
       
        // quanlity control and correctness are done above.
        if(save_dense_flag)
        {
            if(sparsity>0.4)
            {
                //printf("The density of your data is %f. The data will be saved in dense format.\n", sparsity);
                
                string esdfile=string(outFileName)+string(".besd");
                FILE * smr1;
                smr1 = fopen (esdfile.c_str(), "wb");
                if (!(smr1)) {
                    printf("ERROR: failed to open file %s.\n",esdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                uint32_t filetype=DENSE_FILE_TYPE_3;
                vector<int> ten_ints(RESERVEDUNITS);
                ten_ints[0]=filetype;
                if(addn!=-9)
                    printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
                ten_ints[1]=addn;
                ten_ints[2]=(int)esiNum;
                ten_ints[3]=(int)epiNum;
                for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
                
                uint64_t bsize=(uint64_t)esiNum<<1;
                float* buffer=(float*)malloc (sizeof(float)*bsize);
                if (NULL == buffer) {
                    fprintf(stderr, "Malloc failed\n");
                    exit(-1);
                }
                map<string, int>::iterator iter;
                double disp=0;
                for(int j=0;j<epiNum;j++)
                {
                    progress(j, disp, (int)epiNum);
            
					int pkey = prbiflst[j].gd;
                    for(int k=0;k<bsize;k++) buffer[k]=-9; //init
					vector<int> rsid(_ttl_rs[pkey].size());
					for (int l = 0; l<_ttl_rs[pkey].size(); l++){
						iter = esi_map.find(_ttl_rs[pkey][l]);
                        if (iter != esi_map.end()) rsid[l]=iter->second;
                        else {
                            printf("ERROR: SNP is not in SNP map. Please report this bug.\n");
                            exit(EXIT_FAILURE);
                        }
                    }
                    for(int l=0;l<rsid.size();l++)
                    {
						buffer[rsid[l]] = _ttl_beta[pkey][l];
						buffer[esiNum + rsid[l]] = _ttl_se[pkey][l];
                    }
                    fwrite (buffer,sizeof(float), bsize, smr1);
                }
                fclose (smr1);
                free(buffer);
                free_probelist(prbiflst);
                free_snplist(snpiflst);
                cout<<"Effect sizes (beta) and SE for "<<epiNum<<" probes and "<<esiNum<<" SNPs have been saved in a binary file [" + esdfile + "]." <<endl;
                
            } else {
                
                //printf("The density of your data is %f. The data will be saved in sparse format.\n", sparsity);
                
                string esdfile=string(outFileName)+string(".besd");
                FILE * smr1;
                smr1 = fopen (esdfile.c_str(), "wb");
                if (!(smr1)) {
                    printf("ERROR: failed to open file %s.\n",esdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                uint32_t filetype=SPARSE_FILE_TYPE_3;
                vector<int> ten_ints(RESERVEDUNITS);
                ten_ints[0]=filetype;
                if(addn!=-9)
                    printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
                ten_ints[1]=addn;
                ten_ints[2]=(int)esiNum;
                ten_ints[3]=(int)epiNum;
                for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
                
                vector<uint64_t> cols((epiNum<<1)+1);;
                uint64_t valNum=0;
                cols[0]=0;
                
                map<string, int>::iterator iter;
                for(int j=0;j<epiNum;j++)
                {
					int pkey = prbiflst[j].gd;
					uint64_t real_num = _ttl_beta[pkey].size();
                    cols[(j<<1)+1]=real_num+cols[j<<1];
                    cols[j+1<<1]=(real_num<<1)+cols[j<<1];
                    valNum+=real_num*2;
                }
                fwrite (&valNum,sizeof(uint64_t), 1, smr1);
                fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
                double disp=0;
                for(int j=0;j<epiNum;j++)
                {
                    progress(j, disp, (int)epiNum);
					int pkey = prbiflst[j].gd;
					vector<uint32_t> rowids(_ttl_rs[pkey].size());
					for (int l = 0; l<_ttl_rs[pkey].size(); l++){
						iter = esi_map.find(_ttl_rs[pkey][l]);
                        if (iter != esi_map.end()) rowids[l]=iter->second;
                        else {
                            printf("ERROR: SNP is not in SNP map. Please report this bug.\n");
                            exit(EXIT_FAILURE);
                        }
                    }
                    fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
                    fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
                }
                disp=0;
                for(int j=0;j<epiNum;j++)
                {
                    progress(j, disp, (int)epiNum);
					int pkey = prbiflst[j].gd;
					fwrite(&_ttl_beta[pkey][0], sizeof(float), _ttl_beta[pkey].size(), smr1);
					fwrite(&_ttl_se[pkey][0], sizeof(float), _ttl_se[pkey].size(), smr1);
                }
                fclose (smr1);
                free_probelist(prbiflst);
                free_snplist(snpiflst);
                cout<<"Effect sizes (beta) and SE for "<<epiNum<<" probes and "<<esiNum<<" SNPs have been saved in a binary file [" + esdfile + "]." <<endl;
            }
        } else {
            // get esd info
            string esdfile=string(outFileName)+string(".besd");
            FILE * smr1;
            smr1 = fopen (esdfile.c_str(), "wb");
            if (!(smr1)) {
                printf("ERROR: failed to open file %s.\n",esdfile.c_str());
                exit(EXIT_FAILURE);
            }
            uint32_t filetype=SPARSE_FILE_TYPE_3;
            vector<int> ten_ints(RESERVEDUNITS);
            ten_ints[0]=filetype;
            if(addn!=-9)
                printf("Saving sample size %d to the file %s.\n", addn, esdfile.c_str());
            ten_ints[1]=addn;
            ten_ints[2]=(int)esiNum;
            ten_ints[3]=(int)epiNum;
            for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
            fwrite (&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
            
            vector<uint64_t> cols((epiNum<<1)+1);;
            vector<uint32_t> rowids;
            vector<float> val;
            cols[0]=0;
            
            map<string, int>::iterator iter;
            
            
            // log file
            FILE* logfile=NULL;
            string logfname = string(outFileName)+".summary";
            logfile=fopen(logfname.c_str(), "w");
            if (!(logfile)) {
                printf("Error: Failed to open log file.\n");
                exit(1);
            }
            string logstr="cis-window:\t"+itos(cis_itvl)+"Kb\ntrans-window:\t"+itos(trans_itvl)+"Kb\np-value threshold of trans:\t"+dtos(transThres)+"\np-value threshold of others:\t"+dtos(restThres)+"\n";
            logstr+="\ncis region is indicated by [Chr, Start bp, End bp, nsnp];\ntrans region is indicated by <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is indicated by (NumSNPs beyond cis and trans).\n";
            
            logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
            fputs(logstr.c_str(),logfile);
            fflush(logfile);
            cis_itvl=cis_itvl*1000;
            trans_itvl=trans_itvl*1000;
            vector<snpinfolst> snpinfoperprb;
            double disp=0;
            for(int j=0;j<epiNum;j++)
            {
               progress(j, disp, (int)epiNum);
				int pkey = prbiflst[j].gd;
                vector<uint32_t> tmprid;
                vector<float> tmpse;
                string curprb=prbiflst[j].probeId;
                snpinfoperprb.clear();
                
				for (int k = 0; k<_ttl_beta[pkey].size(); k++)
                {
                    snpinfolst tmpinfo;
                    strcpy2(&tmpinfo.snprs,_ttl_rs[pkey][k]);
                    strcpy2(&tmpinfo.a1,_ttl_a1[pkey][k]);
                    strcpy2(&tmpinfo.a2,_ttl_a2[pkey][k]);
					tmpinfo.beta = _ttl_beta[pkey][k];
					tmpinfo.se = _ttl_se[pkey][k];
					tmpinfo.snpchr = _ttl_chr[pkey][k];
					tmpinfo.bp = _ttl_bp[pkey][k];
                        snpinfoperprb.push_back(tmpinfo);
                }
                
                snpinfolst* sortptr=&snpinfoperprb[0];
                qsort(sortptr,snpinfoperprb.size(),sizeof(snpinfolst),comp_esi);
                
                probeinfolst prbifo=prbiflst[j];
                vector<int> slct_idx;
                slct_sparse_per_prb(slct_idx, &prbifo, snpinfoperprb,  cis_itvl,  trans_itvl, transThres, restThres,logfile,false); //slct_idx with no order if there are trans-rgeions
                stable_sort(slct_idx.begin(),slct_idx.end());
                vector<string> _rs(slct_idx.size()), _a1(slct_idx.size()),_a2(slct_idx.size());
                vector<float> _beta(slct_idx.size()), _se(slct_idx.size());
                
                for(int l=0;l<slct_idx.size();l++) {
                    _rs[l]=snpinfoperprb[slct_idx[l]].snprs;
                    _a1[l]=snpinfoperprb[slct_idx[l]].a1;
                    _a2[l]=snpinfoperprb[slct_idx[l]].a2;
                    _beta[l]=snpinfoperprb[slct_idx[l]].beta;
                    _se[l]=snpinfoperprb[slct_idx[l]].se;
                }
                
                vector<int> rsid(_rs.size());
                for (int l = 0; l<_rs.size(); l++){
                    iter = esi_map.find(_rs[l]);
                    if (iter != esi_map.end()) rsid[l]=iter->second;
                    else {
                        printf("ERROR: SNP is not in SNP map. Please report this bug.");
                        exit(EXIT_FAILURE);
                    }
                }
                if(rsid.size()!=_rs.size())
                {
                    printf("ERROR: SNP names don't match for probe %s. Please report this bug.", prbiflst[j].probeId);
                    exit(EXIT_FAILURE);
                }
              
                for(int l=0;l<rsid.size();l++)
                {
                    if(fabs(_se[l]+9)>1e-6) // can move this. the NA is controled in slct_sparse_per_prb
                    {
                        val.push_back(_beta[l]);
                        rowids.push_back(rsid[l]);
                        tmpse.push_back(_se[l]);
                        tmprid.push_back(rsid[l]);
                    } else {
                        printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",_rs[l].c_str());
                    }
                }
                for(int k=0;k<tmpse.size();k++)
                {
                    val.push_back(tmpse[k]);
                    rowids.push_back(tmprid[k]);
                }
                uint64_t real_num=tmpse.size();
                cols[(j<<1)+1]=real_num+cols[j<<1];
                cols[j+1<<1]=(real_num<<1)+cols[j<<1];
                
                free_snplist(snpinfoperprb);
            }
            uint64_t valNum=val.size();
            fwrite (&valNum,sizeof(uint64_t), 1, smr1);
            fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
            fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
            fwrite (&val[0],sizeof(float), val.size(), smr1);
            fclose (smr1);
            free_snplist(snpiflst);
            free_probelist(prbiflst);
            
            printf("Summary data of the specified SNPs and probes has been saved in %s.\n", logfname.c_str());
            cout<<"\nEffect sizes (beta) and SE for "<<epiNum<<" Probes have been saved in a binary file [" + esdfile + "]." <<endl;
            fclose(logfile);

        }
        
    }
    
    void read_maf_file(vector<string> &rs,vector<string> &a1,vector<string> &a2,vector< vector<float> > &maf, vector<string> &cohorts, map<string,int> &rs_map, string mafFileName){
        
        cout << "Reading summary information from [" + mafFileName + "]." << endl;
        char tbuf[MAX_LINE_SIZE];
        gzFile gzfile = gzopen(mafFileName.c_str(), "rb");
        if (!(gzfile)) {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     mafFileName.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        int lineNum(0);
        gzgets(gzfile, tbuf, MAX_LINE_SIZE); //head
        if(tbuf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",mafFileName.c_str());
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        int col_num_header = split_string(tbuf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: the headers should start with \"SNP\"");
            exit(EXIT_FAILURE);
        }
        for(int i=3;i<col_num_header;i++)
        {
            cohorts.push_back(vs_buf[i]);
        }
        while(!gzeof(gzfile))
        {
            gzgets(gzfile, tbuf, MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(tbuf, vs_buf, ", \t\n");
                if(col_num!=col_num_header) {
                    printf("ERROR: the number of columns is incorrect in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: SNP is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: AlleleA  is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: AlleleB  is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                rs_map.insert(pair<string,int>(vs_buf[0], lineNum));
                if(rs_map.size()!=lineNum+1)
                {
                    printf("ERROR: duplicate SNP found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                rs.push_back(vs_buf[0]);
                a1.push_back(vs_buf[1]);
                a2.push_back(vs_buf[2]);
                vector<float> tmpmaf;
                for(int i=3;i<col_num;i++)
                {
                    if(vs_buf[i]=="NA" || vs_buf[i]=="na"){
                        tmpmaf.push_back(-9);
                    }
                    else {
                        tmpmaf.push_back(atof(vs_buf[i].c_str()));
                    }
                }
                maf.push_back(tmpmaf);
                lineNum++;
            }
        }
        gzclose(gzfile);
        cout << lineNum << " MAF info to be included from [" + mafFileName + "]." << endl;
    }
    void selective_cpy(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        esdata->_esi_rs.resize(etrait->_esi_include.size());
        esdata->_esi_allele1.resize(etrait->_esi_include.size());
        esdata->_esi_allele2.resize(etrait->_esi_include.size());
        esdata->_esi_gd.resize(etrait->_esi_include.size());
        esdata->_esi_bp.resize(etrait->_esi_include.size());
        esdata->_esi_chr.resize(etrait->_esi_include.size());
        esdata->_esi_freq.resize(etrait->_esi_include.size());
        esdata->_esi_include.resize(etrait->_esi_include.size());
        esdata->_snp_name_map.clear();
        map<int,int> id_map;
        for (int j = 0; j<etrait->_esi_include.size(); j++)
        {
            string rs=etrait->_esi_rs[etrait->_esi_include[j]];
            esdata->_esi_rs[j]=rs;
            esdata->_esi_allele1[j]=etrait->_esi_allele1[etrait->_esi_include[j]];
            esdata->_esi_allele2[j]=etrait->_esi_allele2[etrait->_esi_include[j]];
            esdata->_esi_gd[j]=etrait->_esi_gd[etrait->_esi_include[j]];
            esdata->_esi_bp[j]=etrait->_esi_bp[etrait->_esi_include[j]];
            esdata->_esi_chr[j]=etrait->_esi_chr[etrait->_esi_include[j]];
            esdata->_esi_freq[j]=etrait->_esi_freq[etrait->_esi_include[j]];
            esdata->_esi_include[j]=j;
            esdata->_snp_name_map.insert(pair<string,int>(rs,j));
            id_map.insert(pair<int,int>(etrait->_esi_include[j],j));
        }
        esdata->_snpNum=etrait->_esi_include.size();
        
        // Also do not forget to copy epi information from etrait to esdata.
        esdata->_epi_bp =  etrait->_epi_bp;
        esdata->_epi_gd = etrait->_epi_gd;
        esdata->_epi_chr = etrait->_epi_chr;
        esdata->_epi_gene = etrait->_epi_gene;
        esdata->_epi_orien = etrait->_epi_orien;
        esdata->_epi_prbID = etrait->_epi_prbID;
        esdata->_include = etrait->_include;
        esdata->_probe_name_map = etrait->_probe_name_map;
        esdata->_probNum = etrait->_probNum;
        
        map<int, int>::iterator iter;
        if(etrait->_rowid.empty())
        {
            esdata->_bxz.resize(etrait->_include.size());
            esdata->_sexz.resize(etrait->_include.size());
            for(int i=0;i<etrait->_include.size();i++)
            {
                esdata->_bxz[i].resize(etrait->_esi_include.size());
                esdata->_sexz[i].resize(etrait->_esi_include.size());
            }
            
            for (int j = 0; j<etrait->_esi_include.size(); j++)
            {
                for (int ii = 0; ii<etrait->_include.size(); ii++)
                {
                    // Here the values in etrait->_include should be 0,1,2,....
                    esdata->_bxz[ii][j]=etrait->_bxz[ii][etrait->_esi_include[j]];
                    esdata->_sexz[ii][j]=etrait->_sexz[ii][etrait->_esi_include[j]];
                    
                }
            }
        }
        else
        {
            // Here the values in etrait->_include should be 0,1,2,....
            esdata->_cols.resize((etrait->_include.size()<<1)+1);
            esdata->_cols[0]=0;
            for (int ii = 0; ii<etrait->_include.size(); ii++)
            {
                uint64_t beta_start=etrait->_cols[ii<<1];
                uint64_t se_start=etrait->_cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                int real_num=0;
                for(int j=0;j<numsnps<<1;j++)
                {
                    int ge_rowid=etrait->_rowid[beta_start+j];
                    iter=id_map.find(ge_rowid);
                    if(iter!=id_map.end())
                    {
                        int sid=iter->second;
                        esdata->_rowid.push_back(sid);
                        esdata->_val.push_back(etrait->_val[beta_start+j]);
                        real_num++;
                    }
                }
                esdata->_cols[(ii<<1)+1]=(real_num>>1)+esdata->_cols[ii<<1];
                esdata->_cols[ii+1<<1]=real_num+esdata->_cols[ii<<1];
            }
            esdata->_valNum=esdata->_val.size();
        }
    }
    void write_smr_esi(string outFileName, eqtlInfo* eqtlinfo)
    {
        string epiName=outFileName+".esi";
        FILE* efile=fopen(epiName.c_str(),"w");
        if(efile==NULL) exit(EXIT_FAILURE);
        printf("Saving SNP information ...\n");
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            string chrstr;
            if(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]==23) chrstr="X";
            else if(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]==24) chrstr="Y";
            else chrstr=atosm(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]);
            
            string str=chrstr+'\t'+eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]+'\t'+atos(0)+'\t'+atosm(eqtlinfo->_esi_bp[eqtlinfo->_esi_include[i]])+'\t'+eqtlinfo->_esi_allele1[eqtlinfo->_esi_include[i]]+'\t'+eqtlinfo->_esi_allele2[eqtlinfo->_esi_include[i]]+'\t'+atosm(eqtlinfo->_esi_freq[eqtlinfo->_esi_include[i]])+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                printf("ERROR: in writing file %s .\n", epiName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        fclose(efile);
        printf("%ld SNPs have been saved in the file %s .\n", eqtlinfo->_esi_include.size(), epiName.c_str());
    }
    void write_smr_epi(string outFileName, eqtlInfo* eqtlinfo)
    {
        string epiName=outFileName+".epi";
        FILE* efile=fopen(epiName.c_str(),"w");
        if(efile==NULL) exit(EXIT_FAILURE);
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            string chrstr;
            if(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]==23) chrstr="X";
            else if(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]==24) chrstr="Y";
            else chrstr=atosm(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]);
            
            string str=chrstr+'\t'+eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]+'\t'+atos(0)+'\t'+atosm(eqtlinfo->_epi_bp[eqtlinfo->_include[i]])+'\t'+eqtlinfo->_epi_gene[eqtlinfo->_include[i]]+'\t'+(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]=='*'?"NA":atos(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]))+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                printf("ERROR: in writing file %s .\n", epiName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        fclose(efile);
        printf("%ld probes have been saved in the file %s .\n", eqtlinfo->_include.size(), epiName.c_str());
    }

}
