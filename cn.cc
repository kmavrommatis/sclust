//########################################################################################
//
//Sclust: Mutation Clustering using Smoothing Splines
//Copyright (C) 2017, Peifer Lab, University of Cologne
//
//This program is free software: you can redistribute it and/or modify it under 
//the terms of the GNU General Public License as published by the Free Software 
//Foundation, either version 3 of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful, but 
//WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License along with this
//program. If not, see <http://www.gnu.org/licenses/>.
//########################################################################################


#include <stdio.h>
#include <stdint.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <list>
#include <map>
#include <sstream>
#include <sys/wait.h>

#include "Sclust.h"
#include "numerics/nr.h"
#include "installdir.h"

#define PI 3.141593

using namespace std;

string cn_help = "      SYNOPSIS \n \t Sclust cn  <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -rc       \t read count file (required)\n \
     \t -snp      \t snp file (required)\n  \
     \t -vcf      \t vcf-file of somatic mutation (required for clustering)\n \
     \t -sv       \t structural variants for segmentation (optional)\n \
     \t -o        \t output name (only prefix)\n \
     \t -w        \t window of initial partition in kbp [1]\n\
     \t -min_r    \t minimal number of reads per partition [1000]\n\
     \t -st       \t sex determination threshold [0.2]\n\
     \t -alpha    \t segmentation p-value cutoff, -1 suppress seg. [0.001]\n\
     \t -ms       \t median smoother (0: none) [2]\n\
     \t -ns       \t minimal number of snps for purity estimation [100]\n\
     \t -maxp     \t maximal expected ploidy [4.5]\n\
     \t -minp     \t minimal expected ploidy [1.5]\n\
     \t -minpu    \t minimal purity [0.2]\n \
     \t -f2       \t compute purity from mutations if this flag is set\n \
     \t -phylowgs \t if this flag is set, input for phylowgs is generated\n \
     \t -pyclone  \t if this flag is set, input for pyclone is generated\n \
     \t -as	  \t flag for allele-specific enriched muts_expAF output\n \
     \t -min_seg  \t minimal number of partitions per segment [2]\n\
     \t -max_qp_iter \t max iterations for QP solver [100]\n\n";


bool cmp_profile(seqcn_profile a,seqcn_profile b) {return(a.scale < b.scale);}

double min_sig_test = 0.015;

//prototypes
void compute_seqcn(vector<CN_data> &cn,seqcn_par seqcn,double &a,char &sex,vector<SNP_data> snp_list);
void purity_est(vector<CN_data> &CN, vector<SNP_data> snp_list, vector<vcf_data> vcf,seqcn_par seqcn,double a,char sex);
double theta_corr(double theta,double sig);
double theta_corr_inv(double theta_obs,double sig);
void optimize_purity(vector<p_data> &obs_p,double &purity,double ploidy_exp,double sig,seqcn_par seqcn,double &L_bi,double &L_cn);
void optimize_purity_ploidy(vector<p_data> &obs_p,double sig,seqcn_par seqcn,double &purity,double &ploidy_exp,double &L_bi,double &L_cn);
void est_alleles_seqcn(vector<p_data> &obs_p,double purity,double scale,double sig,seqcn_par &seqcn,double &L_bi,double &L_cn);
void GC_correct(vector<double> n_reads,vector<double> GC,vector<double> &GC_corr);
void segment_seqcn(vector<double> dat,vector<double> sig,vector<long> &segStart,vector<long> &segEnd,double alpha,long minWd);
void GC_variance(vector<double> rss_cn,vector<double> GC,vector<double> &sigCN);
void extract_theta(long mode,vector<CN_data> CN, vector<SNP_data> snp_list,seqcn_par seqcn,vector<p_data> &obs_p,vector<p_data> &out_cn);
void aggregate_cn(double scale,double purity,vector<CN_data> &CN,vector<p_data> out_cn,int mode);
void mut2cn(vector<p_data> obs_p,vector<vcf_data> &vcf);
void compute_af_seg(vector<p_data> &obs_p,vector<vcf_data> &vcf,double purity);
void subclonal_cn(vector<p_data> &cn,double scale,double purity,double sig);
void vcf_expected_AF(string filename,vector<p_data> cn,double p,vector<vcf_data> vcf,seqcn_par seqcn);
void create_seqcn_profile(vector<p_data> obs_p,double sig,seqcn_par &seqcn,vector<seqcn_profile> &profile);
void optimize_fixed_ploidy(vector<p_data> &obs_p,double &purity,double &scale,double sig,seqcn_par seqcn,double ploidy);
void clear_old_seqcn_files(seqcn_par seqcn);
void plot_pdf_output(double purity,double ploidy,double scale,vector<double> chr_size,vector<p_data> obs_p,vector<seqcn_profile> profile,seqcn_par seqcn);
void write_seqcn_output(double purity,double &ploidy,double scale,char sex,vector<CN_data> CN,vector<p_data> out_cn,vector<vcf_data> vcf,vector<double> &chr_size,seqcn_par seqcn);


void cn(int argc, char *argv[])
{
  long i;
  vector<CN_data> CN;
  ifstream in;
  string out_file;
  string fname="";
  string snp_file="";
  string vcf_file="";
  string tmp,sdrop;
  stringstream line;
  seqcn_par seqcn;
  double a;
  char sex;
  string pu_prior="";
  vector<SNP_data> snp_list;
  SNP_data tmp_snp;

  seqcn.min_reads = 1000;
  seqcn.w = 1;
  seqcn.alpha = 0.001;
  seqcn.median_smooth = 2;
  seqcn.min_seg = 2;
  seqcn.out_name = "";
  seqcn.ns = 100;
  seqcn.max_p = 4.5;
  seqcn.min_p = 1.5;
  seqcn.min_pu = 0.2;
  seqcn.st=0.2;
  seqcn.pu_mu=1;
  seqcn.failed=0;
  seqcn.inv_likelihood=0;
  seqcn.mode="";
  seqcn.f2=0;
  seqcn.lambda=1;
  seqcn.sv="";
  seqcn.phylowgs=0;
  seqcn.pyclone=0;
  seqcn.as=0;
  seqcn.max_qp_iter=100;

  int longindex,opt;
  // option definition
  
  static struct option longopts[]={
    {"help"     , 0, 0,  'h'},
    {"help"     , 0, 0,  '?'},
    {"rc"       , 1, 0,    1},
    {"o"        , 1, 0,    2},
    {"w"        , 1, 0,    3},
    {"min_r"    , 1, 0,    4},
    {"alpha"    , 1, 0,    5},
    {"ms"       , 1, 0,    6},
    {"min_seg"  , 1, 0,    7},
    {"snp"      , 1, 0,    8},
    {"sv"       , 1, 0,    9},
    {"ns"       , 1, 0,   10},
    {"maxp"     , 1, 0,   11},
    {"st"       , 1, 0,   12},
    {"minp"     , 1, 0,   13},
    {"vcf"      , 1, 0,   14},
    {"f2"       , 0, 0,   15},
    {"minpu"    , 1, 0,   16},
    {"phylowgs" , 0, 0,   17},
    {"pyclone"  , 0, 0,   18},
    {"as"       , 0, 0,   19},
    {"max_qp_iter" , 1, 0,   20},
    {0, 0, 0, 0}
  }; 

    //parse command line arguments
  while((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
        {
        case 'h':
	  print_header();
	  cerr << cn_help;
	  exit(1);
          break;
        case '?':
	  print_header();
	  cerr << cn_help;
	  exit(1);  
          break;
        case 1:
	  fname=(string)optarg;
	  break;
	case 2:
	  seqcn.out_name = (string)optarg;
	  break;
	case 3:
	  seqcn.w = atof(optarg);
	  break;
	case 4:
	  seqcn.min_reads = atof(optarg);
	  break;
	case 5:
	  seqcn.alpha = atof(optarg);
	  break;
	case 6:
	  seqcn.median_smooth = atol(optarg);
	  break;
	case 7:
	  seqcn.min_seg = atol(optarg);
	  break;
        case 8:
	  snp_file=(string)optarg;
	  break;
	case 9:
	  seqcn.sv=(string)optarg;
	  break;
	case 10:
	  seqcn.ns = atol(optarg);
	  break;
	case 11:
	  seqcn.max_p = atof(optarg);
	  break;
	case 12:
	  seqcn.st = atof(optarg);
	  break;
	case 13:
	  seqcn.min_p = atof(optarg);
	  break;
	case 14:
	  vcf_file = (string)optarg;
	  break;
	case 15:
	  seqcn.f2 = 1;
	  break;
	case 16:
	  seqcn.min_pu = atof(optarg);
	  break;
	case 17:
	  seqcn.phylowgs = 1;
	  break;
	case 18:
	  seqcn.pyclone = 1;
	  break;
        case 19:
          seqcn.as = 1;
          break;
        case 20:
          seqcn.max_qp_iter = atol(optarg);
          break;
	default:
          cerr << "Error: cannot parse arguments.\n";
          exit(1);
	  break;
        }
    }
  
  if((argc-optind)!=1)
    {
      print_header();
      cerr << cn_help;
      exit(1);
    }
  
  if(fname == "" || seqcn.out_name == "" || snp_file=="")
    {
      print_header();
      cerr << cn_help;
      cerr << "Error: please specify all input or output files\n";
      exit(1);
    }

  in.open(fname.c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open " << fname << ".\n";
      exit(1);
    }


  CN_data tmp_cn;
  vcf_data vcf_record;
  vector<vcf_data> vcf_all;

  //check format
  getline(in,tmp,'\n');
  if(tmp.find("#Sclust "+(string)SCLUSTVERSION) < tmp.size())
    {
      //drop header line
      getline(in,tmp,'\n');
      while(in >> tmp_cn.chr >> tmp_cn.start >> tmp_cn.end >> tmp_cn.n_reads_t >> tmp_cn.n_reads_n >> tmp_cn.GC_content)
	{
	  CN.push_back(tmp_cn);
	}
    }
  else if(tmp.find("#PREMUT_VERSION=1.0") < tmp.size())
    {
      while(getline(in,tmp,'\n'))
	{
	  line.str(""); line.clear();
	  line << tmp;
	  line >> sdrop >> tmp_cn.chr >> tmp_cn.start >> tmp_cn.end >>  sdrop >>  sdrop >>  sdrop >> tmp_cn.n_reads_t >> tmp_cn.n_reads_n >> tmp_cn.GC_content;
	  CN.push_back(tmp_cn);
	}
    }
  else
    {
      cerr << "Error: " << fname << "has incorrect format\n";
      exit(1);
    }
  in.close();

  //computing allelic fraction stuff (if file provided)
  if(snp_file!="")
    {
      in.open((snp_file).c_str());
      if(!in.is_open())
	{
	  cerr << "Error: cannot find snp-file.\n";
	  exit(1);
	}
      
      //check format
      getline(in,tmp,'\n');
      if(tmp.find("#Sclust "+(string)SCLUSTVERSION) < tmp.size())
	{
	  //drop header line
	  getline(in,tmp,'\n');
	  while(in >> tmp_snp.target >> tmp_snp.chr >> tmp_snp.pos 
		>> tmp_snp.s_base_t[0] >> tmp_snp.s_base_t[1] >> tmp_snp.s_base_t[2] 
		>> tmp_snp.s_base_t[3] >> tmp_snp.s_base_n[0] >> tmp_snp.s_base_n[1] 
		>> tmp_snp.s_base_n[2] >> tmp_snp.s_base_n[3] >> tmp_snp.allele_A  >> tmp_snp.allele_B)
	    {
	      snp_list.push_back(tmp_snp);
	    }
	}
      else if(tmp.find("#PREMUT_VERSION=1.0") < tmp.size())
	{
	  i=tmp.find(";a=");
	  if(i<tmp.size())
	    a=atof(tmp.substr(i+3).c_str());
	  
	  while(getline(in,tmp,'\n'))
	    {
	      line.str(""); line.clear();
	      line << tmp;
	      line >> tmp_snp.target >> sdrop >> tmp_snp.pos >> tmp_snp.chr; 
	      line >> tmp_snp.s_base_t[0] >> tmp_snp.s_base_t[1] >> tmp_snp.s_base_t[2];
	      line >> tmp_snp.s_base_t[3] >> tmp_snp.s_base_n[0] >> tmp_snp.s_base_n[1];
	      line >> tmp_snp.s_base_n[2] >> tmp_snp.s_base_n[3] >> tmp_snp.allele_A  >> tmp_snp.allele_B;
	      snp_list.push_back(tmp_snp);
	    }
	}
      else
	{
	  cerr << "Error: " << snp_file << "has incorrect format\n";
	  exit(1);
	}
      in.close();
      //do the segnemtation
      compute_seqcn(CN,seqcn,a,sex,snp_list);

      //read vcf-file
      if(vcf_file!="")
	{
	  //read vcf file
	  in.open((vcf_file).c_str());
	  if(!in.is_open())
	    {
	      cerr << "Error: cannot open vcf-file: " << vcf_file << endl;
	      exit(1);
	    }
	  
	  //throw away the header lines
	  getline(in,tmp,'\n');
	  while(tmp[0]=='#')
	    getline(in,tmp,'\n');
	  
	  
	  while(!in.eof())
	    {
	      read_vcf_line(tmp,vcf_record);
	      //filter out only substitions, no snps
	      if(vcf_record.FILTER=="PASS" && vcf_record.DB==0)
		{
		  //only import mutations from autosomes
		  if(vcf_record.chr!="chrX" && vcf_record.chr!="chrY")
		    vcf_all.push_back(vcf_record);
		}
	      getline(in,tmp,'\n');
	    }
	  in.close();
	}
      
      purity_est(CN,snp_list,vcf_all,seqcn,a,sex);

    }  
}


void sv_input(string sv_file,vector<CN_data> cn,vector<long> &seg_start,vector<long> &seg_end,vector<long> seg_map,long end)
{
  ifstream in;
  stringstream line;
  string tmp;
  vector<long> bkp;
  
  string chr;
  long pos,bkp_pre=-1;

  in.open(sv_file.c_str());
  
  if(!in.is_open())
    {
      cerr << "Error: cannot open structural variants.\n";
      exit(1);
    }
  
  for(long i=0;i<seg_start.size()-1;i++)
    bkp.push_back(seg_end[i]);
  
  seg_start.clear(); seg_end.clear();
  
  while(getline(in,tmp,'\n'))
    {
      if(tmp[0]!='#')
	{
	  line.str(""); line.clear();
	  line << tmp;
	  getline(line,tmp,'\t');
	  if(tmp[0]!='c')
	    chr="chr"+tmp;
	  else
	    chr=tmp;
	  getline(line,tmp,'\t');
	  pos=atol(tmp.c_str());
	  
	  for(long i=0;i < cn.size(); i++)
	    {
	      if(cn[i].chr == chr && cn[i].start <= pos && cn[i].end >= pos && seg_map[i]!=end)
		{
		  bkp.push_back(seg_map[i]);
		  break;
		}
	    }
	}
    }
  in.close();
  
  sort(bkp.begin(),bkp.end());
  seg_start.push_back(0);
  for(long i=0;i < bkp.size(); i++)
    {
      if(bkp[i]!=bkp_pre)
	{
	  seg_end.push_back(bkp[i]);
	  seg_start.push_back(bkp[i]+1);
	  bkp_pre=bkp[i];
	}
    }
  seg_end.push_back(end);
}

void compute_seqcn(vector<CN_data> &cn,seqcn_par seqcn,double &a,char &sex,vector<SNP_data> snp_list)
{
  long i,j,win;
  ofstream out;

  //minimal number of reads per partition
  double mean_X=0,mean_Y=0;
  long n_X=0,n_Y=0;
  double min_reads = seqcn.min_reads;
  string out_name = seqcn.out_name;
  vector<double> GC_corr_n,GC_corr_t;
  vector<double> n_reads_n,n_reads_t;
  vector<double> GC;
  vector<double> cn_raw,sig,L;
  vector<long> seg_start,seg_end;
  vector<long> baf_bkp;
  vector<long> part_start,part_end;
  vector<string> part_chr;
  vector<long> seg_map(cn.size(),0);
  //previous chomosome
  string chr_pre = cn[0].chr;
  //temporary read counts per partition
  double tmp_n_reads_n;
  double tmp_n_reads_t;
  double tmp_GC=0;
  long n_GC;
  //mean read density
  double mean_dens=0;
  //segmented copy number (mean)
  vector<double> mean_cn;
  double tmp_mean_cn;
  //temporary partition length
  double tmp_L;
  //total size of genome
  double A=0;
  //total number of reads
  double reads_tot_t=0,reads_tot_n=0;
  //minimal window size in kbp
  double w = seqcn.w;
  vector<double> rss_cn;
  //sex factor
  double sf;
  //partition data
  seg_start.push_back(0);
  part_chr.push_back(cn[0].chr);
  part_start.push_back(cn[0].start);
  win = 0;
  for(i=0;i<cn.size();i++)
    {
      //total size of sequenced genome
      A += (double)(cn[i].end-cn[i].start+1);
      //total aligned reads
      reads_tot_n += (double)cn[i].n_reads_n;
      reads_tot_t += (double)cn[i].n_reads_t;
      
      if(cn[i].chr=="chrX")
	{
	  mean_X+=(double)cn[i].n_reads_t/(double)(cn[i].end-cn[i].start+1);
	  n_X++;
	}

      if(cn[i].chr=="chrY")
	{
	  mean_Y+=(double)cn[i].n_reads_t/(double)(cn[i].end-cn[i].start+1);
	  n_Y++;
	}

      //if chromosome bondary is hit
      if(cn[i].chr != chr_pre)
	{   
	  if(tmp_n_reads_n > 0)
	    {
	      n_reads_n.push_back(tmp_n_reads_n);
	      n_reads_t.push_back(tmp_n_reads_t);
	      GC.push_back(tmp_GC/(double)n_GC);
	      L.push_back(tmp_L);
	      mean_dens += tmp_n_reads_n/tmp_L;
	      seg_map[i] = win;
	      seg_end.push_back(win);
	      win++; 
	      seg_start.push_back(win);
	      tmp_n_reads_n = (double)cn[i].n_reads_n;
	      tmp_n_reads_t = (double)cn[i].n_reads_t;
	      if(cn[i].GC_content >=0)
		{
		  tmp_GC=cn[i].GC_content;
		  n_GC=1;
		}
	      else
		{
		  tmp_GC=0;
		  n_GC=0;
		}
	      tmp_L +=(double)(cn[i].end-cn[i].start+1);
	    }
	  else
	    {
	      //drop window 
	      seg_map[i] = win;
	      seg_end.push_back(win);
	      seg_start.push_back(win+1);
	      tmp_n_reads_n = (double)cn[i].n_reads_n;
	      tmp_n_reads_t = (double)cn[i].n_reads_t;
	      tmp_L +=(double)(cn[i].end-cn[i].start+1);
	    }
	}
      else
	{
	  tmp_n_reads_n += (double)cn[i].n_reads_n;
	  tmp_n_reads_t += (double)cn[i].n_reads_t;
	  if(cn[i].GC_content>=0)
	    {
	      tmp_GC+=cn[i].GC_content;
	      n_GC++;
	    }
	  tmp_L +=(double)(cn[i].end-cn[i].start+1);
	  if((tmp_L >= w*1000 && tmp_n_reads_n >= min_reads) || i == cn.size()-1)
	    {
	      n_reads_n.push_back(tmp_n_reads_n);
	      n_reads_t.push_back(tmp_n_reads_t);
	      L.push_back(tmp_L);
	      GC.push_back(tmp_GC/(double)n_GC);
	      mean_dens += tmp_n_reads_n/tmp_L;
	      seg_map[i] = win;
	      part_end.push_back(cn[i].start);
	      if(i != cn.size()-1)
		{
		  part_chr.push_back(cn[i+1].chr);
		  part_start.push_back(cn[i+1].start);
		}
	      win++; 
	      tmp_L = 0;
	      tmp_GC=0;
	      n_GC=0;
	      tmp_n_reads_n = 0;
	      tmp_n_reads_t = 0;
	    }
	  else
	    {
	      seg_map[i] = win;
	    }
	}
      
      chr_pre = cn[i].chr;
    }
  mean_dens *= 1.0/(double)L.size();
  seg_end.push_back(n_reads_n.size()-1);
  //partitioning done
  
  //read scructural breakpoints into segmentation
  if(seqcn.sv!="")
    sv_input(seqcn.sv,cn,seg_start,seg_end,seg_map,n_reads_n.size()-1);
  
  // X,Y fraction
  mean_X*=1.0/(double)n_X;
  mean_Y*=1.0/(double)n_Y;

  if(mean_Y/mean_X < seqcn.st)
    sex='f';
  else
    sex='m';

  //GC correction
  GC_correct(n_reads_n,GC,GC_corr_n);
  GC_correct(n_reads_t,GC,GC_corr_t);

  //adapt for totally sequenced bases
  a = reads_tot_n/reads_tot_t;
  for(i=0;i<n_reads_n.size();i++)
      n_reads_t[i]  = n_reads_t[i]*a*GC_corr_n[i]/GC_corr_t[i];
  
  //median smooth number of reads in tumor
  if(seqcn.median_smooth > 0)
    median_smooth(n_reads_t,seqcn.median_smooth);
  //median smooth number of reads in normal
  if(seqcn.median_smooth > 0)
    median_smooth(n_reads_n,seqcn.median_smooth);

  //compute data for segmentation
  cn_raw.clear();
  for(i=0;i<n_reads_n.size();i++)
    {
      cn_raw.push_back(n_reads_t[i]/n_reads_n[i]);
    }
  if(seqcn.median_smooth > 0)
    median_smooth(cn_raw,1);
  
  for(i=0;i<n_reads_n.size();i++)
    {
      sig.push_back((1+cn_raw[i])/((L[i]/A)*reads_tot_n*(n_reads_n[i]/L[i])/mean_dens));
    }

  //segment data
  if(seqcn.alpha != -1)
    {
      segment(cn_raw,sig,seg_start,seg_end,seqcn.alpha,seqcn.min_seg);
      seg_end[seg_end.size()-1] = cn_raw.size()-1;
    }

  //compute mean over segments
  mean_cn = vector<double> (n_reads_n.size(),0);
  rss_cn = vector<double> (n_reads_n.size(),0);

  for(i=0;i<seg_start.size();i++)
    {
      tmp_mean_cn = 0;
      for(j=seg_start[i];j <= seg_end[i];j++)
	{
	  tmp_mean_cn += cn_raw[j];
	}
      tmp_mean_cn *= 1.0/((double)(seg_end[i]-seg_start[i]+1));

      for(j=seg_start[i];j <= seg_end[i];j++)
	{
	  mean_cn[j] = 2*tmp_mean_cn;
	  rss_cn[j] = log(cn_raw[j])-log(tmp_mean_cn);
	}
    }

  //write coarse copy numbers
  for(i=0;i<cn.size();i++)
    {
      if(cn[i].chr=="chrX" && sex=='f')
	sf=1;
      else if(cn[i].chr=="chrX" && sex=='m')
	sf=0.5;
      else if(cn[i].chr=="chrY" && sex=='f')
	sf=0;
      else if(cn[i].chr=="chrY" && sex=='m')
	sf=0.5;
      else
	sf=1;

      cn[i].cn_coarse = sf*mean_cn[seg_map[i]];
    }


  vector<double> cn_raw_log;

  for(i=0;i<cn_raw.size();i++)
    cn_raw_log.push_back(log(cn_raw[i]));
  //compute GC dependent variances for refining the segmentation
  GC_variance(rss_cn,GC,sig);


  if(seqcn.alpha != -1)
    {
      segment_seqcn(cn_raw_log,sig,seg_start,seg_end,seqcn.alpha,seqcn.min_seg);
      seg_end[seg_end.size()-1] = cn_raw.size()-1;
    }

 //compute mean over segments
  mean_cn = vector<double> (n_reads_n.size(),0);
  rss_cn = vector<double> (n_reads_n.size(),0);

  for(i=0;i<seg_start.size();i++)
    {
      tmp_mean_cn = 0;
      for(j=seg_start[i];j <= seg_end[i];j++)
	{
	  tmp_mean_cn += cn_raw[j];
	}
      tmp_mean_cn *= 1.0/((double)(seg_end[i]-seg_start[i]+1));

      for(j=seg_start[i];j <= seg_end[i];j++)
	{
	  mean_cn[j] = 2*tmp_mean_cn;
	  rss_cn[j] = cn_raw[j]-tmp_mean_cn;
	}
    }


  for(i=0;i<cn.size();i++)
    {
      if(cn[i].chr=="chrX" && sex=='f')
	sf=1;
      else if(cn[i].chr=="chrX" && sex=='m')
	sf=0.5;
      else if(cn[i].chr=="chrY" && sex=='f')
	sf=0;
      else if(cn[i].chr=="chrY" && sex=='m')
	sf=0.5;
      else
	sf=1;

      cn[i].cn = sf*mean_cn[seg_map[i]];
    }
}

void purity_est(vector<CN_data> &CN, vector<SNP_data> snp_list, vector<vcf_data> vcf,seqcn_par seqcn,double a,char sex)
{
  long i,j,k,iCN;
  long n_snp_list = snp_list.size();
  double min_coverage = 15;
  
  double L_bi,L_cn,ploidy_exp;

  double subcl_cn_p;
  double cov_n;
  double pi_A_t,pi_A_n;
  double tau_max = 0.15;
  vector<double> xi_n;
  double scale;
  vector<p_data> obs_p,out_cn;
  vector<p_cluster_data> p_cluster;
  seqcn_profile tmp_profile;
  vector<double> chr_size;
  ifstream in;
  string tmp;
  stringstream line;
  double ploidy = 0;

  ofstream out;

  //estimate fold back standard deviation from normal
  double sig=0;

  xi_n.clear();
  for(k=0;k<n_snp_list;k++)
    {
      i=snp_list[k].allele_A; j=snp_list[k].allele_B;
      if(snp_list[k].s_base_n[i]+snp_list[k].s_base_n[j] >= min_coverage && 
	 snp_list[k].s_base_t[i]+snp_list[k].s_base_t[j] >= min_coverage)
	{
	  cov_n = (double)(snp_list[k].s_base_n[i]+snp_list[k].s_base_n[j]);
	  pi_A_n=(double)snp_list[k].s_base_n[i]/cov_n;
	  if(pow(0.5-fabs(pi_A_n-0.5),2) >= tau_max*tau_max)
	    {
	      
	      xi_n.push_back(4*(pi_A_n-0.5)*(pi_A_n-0.5));
	    }
	}
    }
  sig = sqrt(median(xi_n)*a*PI/8.0);

  //cerr << "Fold-back standard deviation: " << sig << endl;
  xi_n.clear();
  
  //fishing out the data
  extract_theta(2,CN,snp_list,seqcn,obs_p,out_cn);
  
  long n_p = obs_p.size();
  double purity;

  if(vcf.size()>0)
    mut2cn(obs_p,vcf);

   //create profile for plotting
   vector<seqcn_profile> profile;
   
   optimize_purity_ploidy(obs_p,sig,seqcn,purity,ploidy_exp,L_bi,L_cn);
   scale=ploidy_exp/(2*purity);

   //invariant copy number
   if(purity==0.2)
     seqcn.inv_likelihood=1;

   if(seqcn.inv_likelihood==1 || seqcn.f2 == 1)
     { 
       purity=1; scale=1;
       if(vcf.size()>0)
	 {
	   vector<double> ccf;
	   vector<long> nmut;
	   double max_ccf=0;
	   long nmut_tot=0;
	   ploidy=seqcn.min_p;
	   for(j=0;j<obs_p.size();j++)
	     obs_p[j].is_subclonal=0;

	   for(j=0;j<5;j++)
	     {
	       ccf.clear(); nmut.clear();  max_ccf=0; nmut_tot=0;
	       est_alleles_seqcn(obs_p,purity,scale,sig,seqcn,L_bi,L_cn);
	       mut2cn(obs_p,vcf);
	       vcf_expected_AF(seqcn.out_name,obs_p,purity,vcf,seqcn);
	       int cluster_exit_code = system(((string)"Sclust cluster -i "+seqcn.out_name+" --max_qp_iter "+to_string(seqcn.max_qp_iter)).c_str());
	       if(cluster_exit_code != 0) {
	           // Propagate the exit code from cluster subprocess (e.g., exit code 2 for QP iteration limit)
	           exit(WEXITSTATUS(cluster_exit_code));
	       }
	       in.open((seqcn.out_name+"_mclusters.txt").c_str());
	       if(!in.is_open())
		 {
		   cerr << "Error: cannot open mutation clusters\n";
		   exit(1);
		 }
	       getline(in,tmp,'\n');
	       while(getline(in,tmp,'\n'))
		 {
		   line.str(""); line.clear();
		   line << tmp;
		   getline(line,tmp,'\t');
		   getline(line,tmp,'\t');
		   ccf.push_back(atof(tmp.c_str()));
		   getline(line,tmp,'\t');
		   getline(line,tmp,'\t');
		   nmut.push_back(atol(tmp.c_str()));
		   nmut_tot+=atol(tmp.c_str());
		 }

	       for(i=0;i<ccf.size();i++)
		 {
		   if(ccf[i] > max_ccf && (double)nmut[i]/(double)nmut_tot >= 0.05)
		     max_ccf=ccf[i];
		 }
	       in.close();
	       for(i=0;i<ccf.size();i++)
		 {
		   if(ccf[i]==max_ccf)
		     {
		       purity=purity*ccf[i];
		       break;
		     }
		 }
	       if(purity > 1)
		 purity=1;
	       scale=ploidy/(2*purity);
	       est_alleles_seqcn(obs_p,purity,scale,sig,seqcn,L_bi,L_cn);
	       ploidy=seqcn.ploidy;
	       ploidy_exp=ploidy;
	     }
	 }
       else
	 seqcn.failed=1;
     }

   extract_theta(1,CN,snp_list,seqcn,obs_p,out_cn);
   est_alleles_seqcn(out_cn,purity,scale,sig,seqcn,L_bi,L_cn);
   subclonal_cn(out_cn,scale,purity,sig);
   if(seqcn.alpha != -1)
     aggregate_cn(scale,purity,CN,out_cn,0);
   else
     aggregate_cn(scale,purity,CN,out_cn,1);
   extract_theta(3,CN,snp_list,seqcn,obs_p,out_cn);

   //second optimization round
   optimize_purity_ploidy(obs_p,sig,seqcn,purity,ploidy_exp,L_bi,L_cn);
   scale=ploidy_exp/(2*purity);
   est_alleles_seqcn(out_cn,purity,scale,sig,seqcn,L_bi,L_cn);
   subclonal_cn(out_cn,scale,purity,sig);
       
   est_alleles_seqcn(out_cn,purity,scale,sig,seqcn,L_bi,L_cn);
   est_alleles_seqcn(obs_p,purity,scale,sig,seqcn,L_bi,L_cn);

   //create profile
   create_seqcn_profile(obs_p,sig,seqcn,profile);
   //subclonal copy number
   subclonal_cn(out_cn,scale,purity,sig);
   
   //write allelic fractions etc
   clear_old_seqcn_files(seqcn);
   
   if(seqcn.failed==0)
     {
       //determine the mode of the run
       if(seqcn.inv_likelihood==1 && seqcn.f2 == 0)
	 seqcn.mode="invariant";
       else if(seqcn.f2 == 1)
	 seqcn.mode="forced";
       else
	 seqcn.mode="optimum";
       
       seqcn.L_bi=L_bi;
       seqcn.L_cn=L_cn;
       //genomic imbalances
       for(i=0;i<out_cn.size();i++)
	 out_cn[i].theta_corr=theta_corr_inv(out_cn[i].theta_t,sig);

       write_seqcn_output(purity,ploidy,scale,sex,CN,out_cn,vcf,chr_size,seqcn);
       
       //adding the final minimum to the profile for plotting
       tmp_profile.scale = scale;
       tmp_profile.purity = purity;
       tmp_profile.ploidy = ploidy;
       profile.push_back(tmp_profile);
       sort(profile.begin(),profile.end(),cmp_profile);

       plot_pdf_output(purity,ploidy,ploidy_exp,chr_size,obs_p,profile,seqcn);
     }
   else
     cout << "Sclust failed for sample: " << seqcn.out_name << ".\n";
}

bool cmp_p_data_theta(p_data a,p_data b) {return(a.theta_t > b.theta_t);}

void normal_approx_xi(vector<double> xi,double &theta,double &theta_sig)
{
  long i;
  long n = xi.size();
  double xi_sig=0;

  theta = sqrt(median(xi));
  theta_sig  = 0;
  
  for(i=0;i<n;i++)
    xi_sig += (sqrt(xi[i])-theta)*(sqrt(xi[i])-theta);

  xi_sig *=1.0/(double)(n*n);

  theta_sig = DMAX(sqrt(xi_sig),min_sig_test);
}


void extract_theta(long mode,vector<CN_data> CN, vector<SNP_data> snp_list,seqcn_par seqcn,vector<p_data> &obs_p,vector<p_data> &out_cn)
{
  //modes: 1: using refined copy numbers
  //modes: 2: using coarse copy numbers
  //modes: 3: using aggregated copy numbers

  long i,j,k;
  double cn_pre;
  long n_snp_list = snp_list.size();
  long n_seg =0;
  long start,end;
  string chr,chr_pre;
  double min_coverage = 15;
  long min_snps=20;

  long idum = -time(NULL);
  p_data tmp_p;
  double cov_t,cov_n;
  double pi_A_t,pi_A_n;
  double tau_max = 0.15;
  vector<double> xi_t,xi_n;
  vector<double> cn;

  //copy number vector
  if(mode==1)
    {
      for(i=0;i<CN.size();i++)
	cn.push_back(CN[i].cn);
    }
  else if(mode==2)
    {
      for(i=0;i<CN.size();i++)
	cn.push_back(CN[i].cn_coarse);
    }
  else if(mode==3)
    {
      for(i=0;i<CN.size();i++)
	cn.push_back(CN[i].cn_aggregated);
    }
  else
    {
      cerr << "Error: unknown mode: " << mode << endl;
      exit(1);
    }

  obs_p.clear();
  out_cn.clear();

  k=0;
  cn_pre = cn[snp_list[k].target];
  chr_pre = CN[snp_list[k].target].chr;
  chr = CN[snp_list[k].target].chr;
  start = CN[snp_list[k].target].start;

  for(k=0;k<n_snp_list;k++)
    {
      i=snp_list[k].allele_A; j=snp_list[k].allele_B;
      if((cn_pre!=cn[snp_list[k].target] || k == n_snp_list-1 || chr_pre!=CN[snp_list[k].target].chr) && n_seg >= min_snps)
	{
	  if(n_seg >= seqcn.ns && chr != "chrX" && chr != "chrY" && cn_pre <= 5)
	    {
	      //output stuff
	      tmp_p.chr = chr_pre;
	      tmp_p.start = start;
	      tmp_p.end = end;
	      tmp_p.cn = cn_pre;
	      tmp_p.n_snps = n_seg;
	      normal_approx_xi(xi_t,tmp_p.theta_t,tmp_p.theta_sig_t);
	      normal_approx_xi(xi_n,tmp_p.theta_n,tmp_p.theta_sig_n);
	      obs_p.push_back(tmp_p);
	      out_cn.push_back(tmp_p);
	    }
	  else
	    {
	      //output stuff
	      tmp_p.chr = chr_pre;
	      tmp_p.start = start;
	      tmp_p.end = end;
	      tmp_p.cn = cn_pre;
	      tmp_p.n_snps = n_seg;
	      normal_approx_xi(xi_t,tmp_p.theta_t,tmp_p.theta_sig_t);
	      normal_approx_xi(xi_n,tmp_p.theta_n,tmp_p.theta_sig_n);
	      out_cn.push_back(tmp_p);
	    }

	  n_seg =0;
	  tmp_p.cn_index.clear();
	  xi_t.clear(); xi_n.clear();
	  cn_pre = cn[snp_list[k].target];
	  chr_pre = CN[snp_list[k].target].chr;
	  chr = CN[snp_list[k].target].chr;
	  start = CN[snp_list[k].target].start;
	  
	  if(snp_list[k].s_base_n[i]+snp_list[k].s_base_n[j] >= min_coverage && 
	     snp_list[k].s_base_t[i]+snp_list[k].s_base_t[j] >= min_coverage)
	    {
	      cov_t = (double)(snp_list[k].s_base_t[i]+snp_list[k].s_base_t[j]);
	      cov_n = (double)(snp_list[k].s_base_n[i]+snp_list[k].s_base_n[j]);
	      pi_A_t=(double)snp_list[k].s_base_t[i]/cov_t;
	      pi_A_n=(double)snp_list[k].s_base_n[i]/cov_n;
	      
	      if(pow(0.5-fabs(pi_A_t-0.5),2)+pow(0.5-fabs(pi_A_n-0.5),2) >= tau_max*tau_max)
		{
		  //new snp to add
		  tmp_p.cn_index.push_back(snp_list[k].target);
		  n_seg++;
		  xi_t.push_back(4*(pi_A_t-0.5)*(pi_A_t-0.5));
		  xi_n.push_back(4*(pi_A_n-0.5)*(pi_A_n-0.5));
		}
	    }
	}
      else
	{
	  if(CN[snp_list[k].target].end-start >=0)
	    end = CN[snp_list[k].target].end;
	  cn_pre = cn[snp_list[k].target];
	  if(snp_list[k].s_base_n[i]+snp_list[k].s_base_n[j] >= min_coverage && 
	     snp_list[k].s_base_t[i]+snp_list[k].s_base_t[j] >= min_coverage)
	    {
	      cov_t = (double)(snp_list[k].s_base_t[i]+snp_list[k].s_base_t[j]);
	      cov_n = (double)(snp_list[k].s_base_n[i]+snp_list[k].s_base_n[j]);
	      pi_A_t=(double)snp_list[k].s_base_t[i]/cov_t;
	      pi_A_n=(double)snp_list[k].s_base_n[i]/cov_n;
	      if(pow(0.5-fabs(pi_A_t-0.5),2)+pow(0.5-fabs(pi_A_n-0.5),2) >= tau_max*tau_max)
		{
		  //new snp to add
		  tmp_p.cn_index.push_back(snp_list[k].target);
		  n_seg++;
		  xi_t.push_back(4*(pi_A_t-0.5)*(pi_A_t-0.5));
		  xi_n.push_back(4*(pi_A_n-0.5)*(pi_A_n-0.5));
		}
	    }
	}
    }

}

void aggregate_cn(double scale,double purity,vector<CN_data> &CN,vector<p_data> out_cn,int mode)
{
  long i,j,n;
  string chr_pre;
  long iCN_pre,iCN;
  int sub_pre;
  double theta_exp_pre;
  double cn_mean=0;
  long n_mean=0;
  vector<long> start,end;

  n=CN.size();
  vector<int> is_subclonal(n,-1);
  vector<double> theta_exp(n,-1);
  for(i=0;i<out_cn.size();i++)
    {
      for(j=0;j<out_cn[i].cn_index.size();j++)
	{
	  is_subclonal[out_cn[i].cn_index[j]]=(int)out_cn[i].is_subclonal;
	  theta_exp[out_cn[i].cn_index[j]]=out_cn[i].theta_exp;
	}
    }

  for(i=0;i<n;i++)
    {
      if(is_subclonal[i]!=-1)
	{
	  sub_pre=is_subclonal[i];
	  break;
	}
    }

  for(i=0;i<n;i++)
    {
      if(theta_exp[i]!=-1)
	{
	  theta_exp_pre=theta_exp[i];
	  break;
	}
    }

  //fill up -1 with the previous values
  for(i=0;i<n;i++)
    {
      if(is_subclonal[i]==-1)
	is_subclonal[i]=sub_pre;
      else
	sub_pre=is_subclonal[i];

      if(theta_exp[i]==-1)
	theta_exp[i]=theta_exp_pre;
      else
	theta_exp_pre=theta_exp[i];
	  
    }
  sub_pre=is_subclonal[0];
  theta_exp_pre=theta_exp[0];
  chr_pre=CN[0].chr;
  iCN_pre=(long)((CN[0].cn-2*(1-purity))*scale+0.5);
  start.push_back(0);

  if(mode==0)
    {
      for(i=0;i<n;i++)
	{
	  iCN=(long)((CN[i].cn-2*(1-purity))*scale+0.5);
	  if(chr_pre!=CN[i].chr || theta_exp_pre!=theta_exp[i] ||sub_pre!=is_subclonal[i] || iCN!=iCN_pre)
	    {
	      end.push_back(i-1); 
	      chr_pre=CN[i].chr;
	      iCN_pre=iCN;
	      theta_exp_pre=theta_exp[i];
	      sub_pre=is_subclonal[i];
	      start.push_back(i);
	    }
	}
  
      end.push_back(n-1);
      
      for(i=0;i<start.size();i++)
	{
	  cn_mean=0; n_mean=0;
	  for(j=start[i];j<=end[i];j++)
	    {
	      cn_mean+=CN[j].cn;
	      n_mean++;
	    }
	  cn_mean*=1.0/(double)n_mean;
	  for(j=start[i];j<=end[i];j++)
	    CN[j].cn_aggregated=cn_mean;
	}
    }
  else
    {
      for(i=0;i<n;i++)
	CN[i].cn_aggregated=CN[i].cn;
    }
}


void est_alleles_seqcn(vector<p_data> &obs_p,double purity,double scale,double sig,seqcn_par &seqcn,double &L_bi,double &L_cn)
{
  double min_dist,dist,ploidy,sl;
  long i,j,k,i_min,cn,m;
  long cn_max = 50,min_snps = 5;
  double ls,lls,dL_bi,dL_cn;

  L_bi=0; L_cn=0;

  vector< vector<long> > mapA(cn_max+1,vector<long>(0,0));
  vector< vector<long> > mapB(cn_max+1,vector<long>(0,0));

  //create map
  for(k=1;k<=cn_max;k++)
    {
      for(i=0;i<=(long)(k/2);i++)
	{
	  mapA[k].push_back(k-i);
	  mapB[k].push_back(i);
	}
    }

  //compute dist scores
  for(i=0;i<obs_p.size();i++)
    {
      //iCN type of copy numbers
      cn = (long)((obs_p[i].cn-2*(1-purity))*scale+0.5);
      if(cn > cn_max)
	cn = cn_max-1;
      if(cn<=0)
	cn=0;
      if(obs_p[i].n_snps >= min_snps && cn >=1)
	{
	  i_min = 0; min_dist = 1e10;
	
	  for(j=0;j<mapA[cn].size();j++)
	    {
	      dist = fabs(obs_p[i].theta_t-theta_corr(purity*(mapA[cn][j]-mapB[cn][j])/(2+purity*(cn-2)),sig));
	      if(dist < min_dist)
		{
		  min_dist = dist;
		  i_min = j;
		}
	    }

	  ls=1;
	  if(mapA[cn].size()==1)
	    ls = fabs(theta_corr(purity*(mapA[cn][0]-mapB[cn][0])/(2+purity*(cn-2)),sig)-theta_corr(0,sig));
	  else
	    {
	      for(j=0;j<mapA[cn].size()-1;j++)
		{
		  lls=fabs(obs_p[i].theta_t-theta_corr(purity*(mapA[cn][j]-mapB[cn][j])/(2+purity*(cn-2)),sig)-theta_corr(purity*(mapA[cn][j+1]-mapB[cn][j+1])/(2+purity*(cn-2)),sig));
		  if(lls < ls)
		    ls = lls;
		}
	    }


	  obs_p[i].A=mapA[cn][i_min];
	  obs_p[i].B=mapB[cn][i_min];
	  obs_p[i].ls = 1.0/ls;
	  obs_p[i].theta_exp = theta_corr(purity*(mapA[cn][i_min]-mapB[cn][i_min])/(2+purity*(cn-2)),sig);
	  
	}
      else if(cn==0)
	{
	  obs_p[i].A = 0;
	  obs_p[i].B = 0;
	  obs_p[i].theta_exp = theta_corr(0,sig);
	}
      else
	{
	  obs_p[i].A = -1;
	  obs_p[i].B = -1;
	  obs_p[i].theta_exp = -1;
	}
    }

  //compute likelhood
 
  for(i=0;i<obs_p.size();i++)
    {
      if(obs_p[i].A != -1 && obs_p[i].B != -1)
	{
	  cn = obs_p[i].A+obs_p[i].B;
	  obs_p[i].iCN = cn;
	  dL_bi=pow(obs_p[i].theta_t-theta_corr(purity*(obs_p[i].A-obs_p[i].B)/(2+purity*(cn-2)),sig),2)/(2*obs_p[i].theta_sig_t*obs_p[i].theta_sig_t);
	  dL_cn=pow((obs_p[i].cn-2*(1-purity))*scale-obs_p[i].iCN,2)/(2*obs_p[i].theta_sig_t*obs_p[i].theta_sig_t);
	  L_bi+=dL_bi;
	  L_cn+=dL_cn;
	}
    }

  ploidy=0;sl=0;
  for(i=0;i<obs_p.size();i++)
    {
      cn = (long)((obs_p[i].cn-2*(1-purity))*scale+0.5);
      if(cn<=0)
	cn=1;
      ploidy+=(double)(obs_p[i].end-obs_p[i].start)*cn;
      sl+=(double)(obs_p[i].end-obs_p[i].start);
    }
  ploidy*=1.0/sl;

  seqcn.ploidy=ploidy;
}


double compute_tau(double theta,double p,long A1,long B1,long A2,long B2)
{
  double nm1,np1;
  double nm2,np2;
  double tau;

  nm1=(double)(A1-B1);
  nm2=(double)(A2-B2);
  np1=(double)(A1+B1);
  np2=(double)(A2+B2);

  tau=(nm2-theta*np2-2*theta*(1-p)/p)/(theta*(np1-np2)-nm1+nm2);
  return(tau);
}

void subclonal_cn(vector<p_data> &cn,double scale,double purity,double sig)
{
  long i,j,A,B;
  double raw_A;
  double raw_B;
  double raw_cn;
  double theta;
  double tau_c,max_tau;
  long max_cl;
  double p=purity;
  
  vector<long> clA,clB;
  vector<double> tau;
  bool sel=0;
  double subcl_cn_p;

  for(i=0;i<cn.size();i++)
    {
      cn[i].is_subclonal=0;
      cn[i].p_subclonal=1;
      cn[i].inconsistent_state=0;
      if(cn[i].n_snps >=5 && cn[i].theta_exp >= 0)
	{
	   
	  subcl_cn_p = 1-erf(fabs(cn[i].theta_t-cn[i].theta_exp)/(sqrt(2)*cn[i].theta_sig_t));
	  cn[i].p_subclonal=subcl_cn_p;
	  if( subcl_cn_p <0.001 && cn[i].theta_exp >= 0)
	    {
	      cn[i].is_subclonal=1;
	    }
	  else
	    {
	      cn[i].is_subclonal=0;
	    }
	}
    }
  
  for(i=0;i<cn.size();i++)
    {
      if(cn[i].is_subclonal == 1)
	{
	  if(cn[i].A==cn[i].B && cn[i].theta_t < cn[i].theta_exp)
	    {
	      cn[i].is_subclonal = 0;
	      cn[i].p_subclonal = 1;
	    }
	  if(cn[i].B==0 && cn[i].theta_t > cn[i].theta_exp)
	    {
	      cn[i].is_subclonal = 0;
	      cn[i].p_subclonal = 1;
	    }
	}

      if(cn[i].is_subclonal == 1)
	{
	  raw_cn = (cn[i].cn-2*(1-purity))*scale;
	  if(raw_cn < 0)
	    raw_cn=0;
	  
	  theta=theta_corr_inv(cn[i].theta_t,sig);
	  
	  raw_A = (theta*(2+p*(raw_cn-2))+p*raw_cn)/(2*p);
	  raw_B = raw_cn-raw_A;
	  if(raw_B < 0)
	    cn[i].inconsistent_state=1;
	  else
	    {
	      clA.clear();
	      clB.clear();
	      tau.clear();
	      A=(long)floor(raw_A);
	      B=(long)floor(raw_B);

	      if(A+B==(long)floor(raw_cn))
		{
		  if((A==cn[i].A && B==cn[i].B) ||(A+1==cn[i].A && B==cn[i].B))
		    {
		      if(A==cn[i].A && B==cn[i].B && A+1 >= B)
			{
			  clA.push_back(A+1);
			  clB.push_back(B);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A+1,B);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}
		      if(A+1==cn[i].A && B==cn[i].B && A >= B)
			{
			  clA.push_back(A);
			  clB.push_back(B);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A,B);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}
		    }

		  if((A==cn[i].A && B==cn[i].B) ||(A==cn[i].A && B+1==cn[i].B))
		    {
		      if(A==cn[i].A && B==cn[i].B && A >= B+1)
			{
			  clA.push_back(A);
			  clB.push_back(B+1);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A,B+1);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}
		      if(A==cn[i].A && B+1==cn[i].B && A >= B)
			{
			  clA.push_back(A);
			  clB.push_back(B);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A,B);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}
		    }
		}
	      else
		{
		  if((A+1==cn[i].A && B==cn[i].B) ||(A+1==cn[i].A && B+1==cn[i].B))
		    {
		      if(A+1==cn[i].A && B==cn[i].B && A+1 >= B+1)
			{
			  clA.push_back(A+1);
			  clB.push_back(B+1);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A+1,B+1);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}

		      if(A+1==cn[i].A && B+1==cn[i].B && A+1 >= B)
			{
			  clA.push_back(A+1);
			  clB.push_back(B);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A+1,B);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}
		    }

		  if((A==cn[i].A && B+1==cn[i].B) ||(A+1==cn[i].A && B+1==cn[i].B))
		    {
		      if(A==cn[i].A && B+1==cn[i].B && A+1 >= B+1)
			{
			  clA.push_back(A+1);
			  clB.push_back(B+1);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A+1,B+1);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}

		      if(A+1==cn[i].A && B+1==cn[i].B && A >= B+1)
			{
			  clA.push_back(A);
			  clB.push_back(B+1);
			  tau_c=compute_tau(theta,p,cn[i].A,cn[i].B,A,B+1);
			  if(tau_c <= 0 || tau_c >=1)
			    tau_c=-1;
			  tau.push_back(tau_c);
			}
		    }
		}

	      
	      max_tau=-1; max_cl=0;
	      for(j=0;j<clA.size();j++)
		{
		  if(tau[j] > max_tau)
		    {
		      max_tau=tau[j];
		      max_cl=j;
		    }
		}

	      if(max_tau > 0)
		{
		  j=max_cl;

		  if(cn[i].B==0)
		    sel=0;
		  else if(clB[j]==0)
		    sel=1;
		  else if(clB[j]+clA[j] > cn[i].A+cn[i].B)
		    sel=1;
		  else
		    sel=0;

		  if(sel==0)
		    {
		      cn[i].clone1_A = cn[i].A;
		      cn[i].clone1_B = cn[i].B;
		      cn[i].clone1_tau=tau[j];
		      cn[i].clone2_A = clA[j];
		      cn[i].clone2_B = clB[j];
		      cn[i].clone2_tau=1-tau[j];
		    }
		  else
		    {
		      cn[i].clone2_A = cn[i].A;
		      cn[i].clone2_B = cn[i].B;
		      cn[i].clone2_tau=tau[j];
		      cn[i].clone1_A = clA[j];
		      cn[i].clone1_B = clB[j];
		      cn[i].clone1_tau=1-tau[j];
		    }

		  if(cn[i].clone1_A+cn[i].clone1_B == cn[i].clone2_A+cn[i].clone2_B)
		    cn[i].is_subclonal = 1;
		  else
		    {
		      cn[i].A = cn[i].clone1_A;
		      cn[i].B = cn[i].clone1_B;
		      cn[i].iCN=cn[i].A+cn[i].B;
		    }
		 
		}
	      else
		{
		  cn[i].inconsistent_state=1;
		}
	    }
	}
    }
}


double theta_corr(double theta,double sig)
{
  double tau = 0.5*theta;

  return(2*tau*erf(tau/(sig*sqrt(2)))+sqrt(8/PI)*sig*exp(-tau*tau/(2*sig*sig)));
}

double theta_corr_inv(double theta_obs,double sig)
{
  double eps=1e-6;
  double x_l=0,x_r=1,x_m=0.5;
  double f_l=theta_obs-theta_corr(x_l,sig);
  double f_r=theta_obs-theta_corr(x_r,sig);
  double f_m=theta_obs-theta_corr(x_m,sig);
  long i=0,n_max=100;

  if(f_l < 0)
    return(0);
  
  if(fabs(f_l) <= eps)
    return(x_l);

  if(fabs(f_r) <= eps)
    return(x_r);
    
  while(fabs(f_m) > eps)
    {
      if(f_l*f_m < 0)
	x_r=x_m;
      else
	x_l=x_m;
      
      x_m=0.5*(x_l+x_r);
      
      f_l=theta_obs-theta_corr(x_l,sig);
      f_r=theta_obs-theta_corr(x_r,sig);
      f_m=theta_obs-theta_corr(x_m,sig);

      if(i==n_max)
	break;
      i++;
    }
  
  return(x_m);
}


void create_seqcn_profile(vector<p_data> obs_p,double sig,seqcn_par &seqcn,vector<seqcn_profile> &profile)
{
  long i,j;
  double max_p = DMAX(5,seqcn.max_p);
  double min_p = 1;
  double dx = 0.05;
  double scale;
  double p_exp;
  ofstream out;
  
  double L_bi,L_cn;
  double max_L=-1,min_L=1e10;
  long iCN;
  double purity;
  double ploidy;
  double norm_ploidy;

  seqcn_profile tmp_profile;
  
  profile.clear();
  
  
  for(p_exp=min_p;p_exp<=max_p;p_exp+=dx)
    {   
      
      optimize_purity(obs_p,purity,p_exp,sig,seqcn,L_bi,L_cn);
      scale = p_exp /(2*purity);

      //to determine if the likelihood is invariant
      if(p_exp <= 3)
	{
	  if(L_bi > max_L)
	    max_L = L_bi;
	  if(L_bi < min_L)
	    min_L = L_bi;
	}

      tmp_profile.ploidy_exp = p_exp;
      tmp_profile.scale = scale;
      tmp_profile.purity = purity;
      tmp_profile.L_bi = L_bi;
      tmp_profile.L_cn = L_cn;

      //compute ploidy
      ploidy = 0; norm_ploidy=0;
      
      for(i=0;i<obs_p.size();i++)
  	{
  	  iCN = (long)((obs_p[i].cn-2*(1-purity))*scale+0.5);
  	  if(obs_p[i].chr!="chrX" && obs_p[i].chr!="chrY" && iCN >=0)
  	    {
  	      ploidy += (double)iCN*((double)(obs_p[i].end - obs_p[i].start));
  	      norm_ploidy += (double)(obs_p[i].end - obs_p[i].start);
  	    }
  	}
      ploidy *= 1.0/norm_ploidy;
      tmp_profile.ploidy = ploidy;
      
      profile.push_back(tmp_profile);
    }

  //output profile
  out.open((seqcn.out_name+"_cn_profile.txt").c_str());
  out << "Ploidy_Exp\tPurity\tPloidy\tScale\tLog-Likelihood-Bi-Allele\tLog-Likelihood-CN\n";
  for(i=0;i<profile.size();i++)
    {
      out << profile[i].ploidy_exp << "\t" << profile[i].purity << "\t" << profile[i].ploidy << "\t" << profile[i].scale;
      out << "\t" << profile[i].L_bi << "\t" << profile[i].L_cn << endl;
    }
  out.close();
}

void optimize_purity(vector<p_data> &obs_p,double &purity,double ploidy_exp,double sig,seqcn_par seqcn,double &L_bi,double &L_cn)
{
 
  double p_min=seqcn.min_pu,L_min=1e10,p;
  
  for(p=seqcn.min_pu;p<=1;p=p+0.01)
    {
      est_alleles_seqcn(obs_p,p,ploidy_exp/(2*p),sig,seqcn,L_bi,L_cn);
      if(L_bi < L_min)
	{
	  p_min=p;
	  L_min=L_bi;
	}
    }
  //output results
  est_alleles_seqcn(obs_p,p_min,ploidy_exp/(2*p_min),sig,seqcn,L_bi,L_cn);
  purity=p_min;
}

void optimize_purity_ploidy(vector<p_data> &obs_p,double sig,seqcn_par seqcn,double &purity,double &ploidy_exp,double &L_bi,double &L_cn)
{

  double p_min=seqcn.min_p,L_min=1e10,p;
  
  for(p=seqcn.min_p;p<=seqcn.max_p;p=p+0.01)
    {
      optimize_purity(obs_p,purity,p,sig,seqcn,L_bi,L_cn);
      if(L_cn < L_min)
	{
	  p_min=p;
	  L_min=L_cn;
	}
    }

  //output results
  optimize_purity(obs_p,purity,p_min,sig,seqcn,L_bi,L_cn);
  ploidy_exp=p_min;
}

void GC_correct(vector<double> n_reads,vector<double> GC,vector<double> &GC_corr)
{
  long nbin=1000;
  long i,j,n;
  double read_median;
  double alpha=1;

  GC_corr.clear();

  vector< vector<double> > GCbin(nbin+1,vector<double>(0,0));
  vector<double> readsBin(nbin+1,0);

  if(n_reads.size()!=GC.size())
    {
      cerr << "Error: read depth data and GC content does not match.\n";
      exit(1);
    }
 
  //do binning for the GC data
  for(i=0;i<GC.size();i++)
    {
      GCbin[(long)(nbin*GC[i])].push_back(n_reads[i]);
    }

  read_median = median(n_reads);
  n=0;
  for(i=0;i<=nbin;i++)
    {
      if(GCbin[i].size()!=0)
	{
	  readsBin[i]=median(GCbin[i]);
	  n++;
	}
    }

  //variables for splines
  double *x=dvector(1,n);
  double *y=dvector(1,n);
  double *w=dvector(1,n);
  double *g=dvector(1,n);
  double *gam=dvector(1,n);

  j=1;
  for(i=0;i<=nbin;i++)
    {
      if(GCbin[i].size()!=0)
	{
	  x[j]=(double)i/(double)nbin;
	  y[j]=readsBin[i];
	  w[j]=sqrt((double)GCbin[i].size());
	  j++;
	}
    }
  
  alpha=1;
  splines_gcv(x,y,w,n,&alpha,g,gam);

  gam[n]=0;
  for(i=n-1;i>=2;i--)
    gam[i]=gam[i-1];
  gam[1]=0;

  for(i=0;i<GC.size();i++)
    {
      GC_corr.push_back(plot_splines(x,g,gam,n,GC[i])/read_median);
    }

  free_dvector(x,1,n);
  free_dvector(y,1,n);
  free_dvector(w,1,n);
  free_dvector(g,1,n);
  free_dvector(gam,1,n);
}

void GC_variance(vector<double> rss_cn,vector<double> GC,vector<double> &sigCN)
{
  long nbin=100;
  long i,j;
  double var;
  double max_var=0;

  sigCN.clear();

  vector< vector<double> > GCbin(nbin+1,vector<double>(0,0));
  vector<double> sigBin(nbin+1,0);

  if(rss_cn.size()!=GC.size())
    {
      cerr << "Error: read depth data and GC content does not match.\n";
      exit(1);
    }
 
  //do binning for the GC data
  for(i=0;i<GC.size();i++)
    {
      GCbin[(long)(nbin*GC[i])].push_back(rss_cn[i]);
    }

  for(i=0;i<=nbin;i++)
    {
      if(GCbin[i].size()>5)
  	{
  	  sigBin[i]=0;
	  for(j=0;j<GCbin[i].size();j++)
	    sigBin[i]+=GCbin[i][j]*GCbin[i][j];
	  sigBin[i]*=1.0/(double)GCbin[i].size();
	  if(sigBin[i] > max_var)
	    max_var = sigBin[i];
  	}
      else
  	sigBin[i] = -1;
    }

   for(i=0;i<GC.size();i++)
    {
      if(sigBin[(long)(nbin*GC[i])]==-1)
	sigCN.push_back(max_var);
      else
	sigCN.push_back(sigBin[(long)(nbin*GC[i])]);
    }
}


void segment_seqcn(vector<double> dat,vector<double> sig,vector<long> &segStart,vector<long> &segEnd,double alpha,long minWd)
{
  long i,n,bkp;
  vector<long> tmpStart,tmpEnd;
  list<long> finalBKP;
  double score,maxScore,S1,S2,n1,n2,SQ1,SQ2,v1,v2;
  n=dat.size();
  
  for(i=0;i<segStart.size();i++)
    {
      tmpStart.push_back(segStart[i]);
      tmpEnd.push_back(segEnd[i]);
      finalBKP.push_back(segStart[i]);
      finalBKP.push_back(segEnd[i]);
    }
   
  //clear seg-vectors
  segStart.clear(); segEnd.clear();
  //main loop
  while(tmpStart.size()!=0)
    {
      S1=0; S2=0;
      SQ1=0; SQ2=0;
      n1=0; n2=0;
      v1=0; v2=0;
      maxScore=-1;
      //initialize sub segment
      for(i=tmpStart[0];i<=tmpStart[0]+2*minWd;i++)
	{
	  if(i-tmpStart[0]<=minWd-1)
	    {
	      S1+=dat[i];
	      SQ1+=dat[i]*dat[i];
	      v1 += sig[i];
	      n1+=1;
	    }
	  else if(i-tmpStart[0] > minWd)
	    {
	      S2+=dat[i];
	      SQ2+=dat[i]*dat[i];
	      v2 += sig[i];
	      n2+=1;
	    }
	}
      for(i=1;i<=(tmpEnd[0]-tmpStart[0])-2*minWd;i++)
	{
	  //compute t-statistic for all sub-segments
	  S1+=dat[tmpStart[0]+minWd+i-1];
	  SQ1+=dat[tmpStart[0]+minWd+i-1]*dat[tmpStart[0]+minWd+i-1];
	  v1+=sig[tmpStart[0]+minWd+i-1];
	  S1-=dat[tmpStart[0]+i-1];
	  SQ1-=dat[tmpStart[0]+i-1]*dat[tmpStart[0]+i-1];
	  v1-=sig[tmpStart[0]+i-1];
	  S2+=dat[tmpStart[0]+2*minWd+i];
	  SQ2+=dat[tmpStart[0]+2*minWd+i]*dat[tmpStart[0]+2*minWd+i];
	  v2+=sig[tmpStart[0]+2*minWd+i];
	  S2-=dat[tmpStart[0]+minWd+i];
	  SQ2-=dat[tmpStart[0]+minWd+i]*dat[tmpStart[0]+minWd+i];
	  v2-=sig[tmpStart[0]+minWd+i];
	
	  score=fabs(S1/n1-S2/n2)/sqrt(v1/(n1*n1)+v2/(n2*n2));
	  if(score>maxScore)
	    {
	      maxScore=score;
	      bkp=tmpStart[0]+minWd+i;
	    }
	}
      //is maxScore significantly high?
      if(erfc(maxScore/sqrt(2)) <=alpha && tmpEnd[0]-tmpStart[0]+1 >=minWd)
	{
	  finalBKP.push_back(bkp);
	  //subdivide segment 
	  tmpStart.push_back(tmpStart[0]);
	  tmpEnd.push_back(bkp);
	  tmpStart.push_back(bkp+1);
	  tmpEnd.push_back(tmpEnd[0]);
	  maxScore = 0;
	}
      tmpStart.erase(tmpStart.begin());
      tmpEnd.erase(tmpEnd.begin());
    }
  //sort break points
  finalBKP.sort();
  //get segments from break points
  list<long>::iterator it;
  for(it=finalBKP.begin();it!=finalBKP.end();it++)
    {
      if(*it!=n-1)
	{
	  if(*it==0)
	    {
	      segStart.push_back(0);
	      it++;
	      segEnd.push_back(*it);
	      it--;
	    }
	  else
	    {
	      segStart.push_back(*it+1);
	      it++;
	      segEnd.push_back(*it);
	      it--;
	    }
	}

    }
}


void mut2cn(vector<p_data> obs_p,vector<vcf_data> &vcf)
{
  long i,j;

  for(i=0;i<vcf.size();i++)
    {
      vcf[i].index=-1;
      //quick and dirty can be optimized
      for(j=0;j<obs_p.size();j++)
	{
	  if(vcf[i].chr==obs_p[j].chr && vcf[i].pos >= obs_p[j].start && vcf[i].pos <= obs_p[j].end)
	    {
	      vcf[i].index=j;
	      break;
	    }
	}
    }
  
}


bool cmp_vcf_index(vcf_data a,vcf_data b){return(a.index < b.index);}

void compute_af_seg(vector<p_data> &obs_p,vector<vcf_data> &vcf,double p)
{
  long i,j,m,k;
  long cn,n_max;
  vector< vector<double> > m_af;
  vector< vector<double> > m_sig2;
  vector<long> n_af;
  vector<long> vcf_index;
  long lower,upper;
  double af;
  double sig2;

  sort(vcf.begin(),vcf.end(),cmp_vcf_index);

  for(i=0;i<vcf.size();i++)
    vcf_index.push_back(vcf[i].index);
  
  for(i=0;i<obs_p.size();i++)
    {
      obs_p[i].seg_af.clear();
      obs_p[i].n_seg_af.clear();
      lower= lower_bound (vcf_index.begin(), vcf_index.end(),i)-vcf_index.begin();
      upper= upper_bound (vcf_index.begin(), vcf_index.end(),i)-vcf_index.begin();

      if(lower<upper)
	{
	  cn=(long)(obs_p[i].A+obs_p[i].B);
	  m_af.clear(); n_af.clear();
	  n_af = vector<long>(cn+2,0); 
	  m_af = vector< vector<double> > (cn+2,vector<double>(0,0));
	  m_sig2 = vector< vector<double> > (cn+2,vector<double>(0,0));
	  for(j=lower;j<upper;j++)
	    {
	      af=vcf[j].AF;
	      sig2=af*(1-af)/((double)vcf[j].DP);
	      if(af > p*(double)(cn+0.5)/(2*(1-p)+p*(double)cn))
		{
		  m_af[cn+1].push_back(af);
		  m_sig2[cn+1].push_back(sig2);
		  n_af[cn+1]++;
		}
	      else if(af < (p*0.5)/(2*(1-p)+p*(double)cn))
		{
		  m_af[0].push_back(af);
		  m_sig2[0].push_back(sig2);
		  n_af[0]++;
		}
	      else
		{
		  m=(long)(af*(2*(1-p)+p*(double)cn)/p+0.5);
		  m_af[m].push_back(af);
		  m_sig2[m].push_back(sig2);
		  n_af[m]++;
		}
	    }

	  n_max=0;
	  for(j=0;j<cn+2;j++)
	    {
	      if(n_af[j] > n_max)
		n_max = n_af[j]; 
	    }

	  if(n_max==0)
	    n_max=1; //to prevent devision by zero

	  for(j=0;j<cn+2;j++)
	    {
	      if(m_af[j].size() >=20)
		{
		  //very heuristic at the moment
		  if((double)n_af[j]/(double)n_max >= 0.25 || fabs(median(m_af[j])-(double)j*p/(2*(1-p)+p*(double)cn)) <= 0.02)
		    {
		      if(m_af[j].size()>1)
			{
			  obs_p[i].seg_af.push_back(median(m_af[j]));
			  sig2=0;
			  for(k=0;k<m_sig2[j].size();k++)
			    sig2+=m_sig2[j][k];
			  obs_p[i].seg_af_sig.push_back(sqrt(sig2));
			}
		    }
		}
	    }
	}
      else
	obs_p[i].seg_af.clear();
    }
}


void vcf_expected_AF(string filename,vector<p_data> cn,double p,vector<vcf_data> vcf,seqcn_par seqcn)
{
  ofstream out,outp,outpy,outas;
  long i,j,ind,s=0;
  double m,m_raw;
  double icn;
  double af,af_sc;
  double tau;
  double pval,ccf,sig;
  double cnA, cnB;                           // clonal case: copy numbers of alleles A and B
  double cn1A, cn1B, cn2A, cn2B, tau1, tau2; // subclonal case: cn's of alleles A, B of clones 1 and 2 
  //phylowgs input generation
  if(seqcn.phylowgs==1)
    {
      outp.open((filename+".phylowgs.cnv_data.txt").c_str());
      outp.close();
      outp.open((filename+".phylowgs.ssm_data.txt").c_str());
      outp << "id\tgene\ta\td\tmu_r\tmu_v\n";
    }

  //pyclone input generation
  if(seqcn.pyclone==1)
    {
      outpy.open((filename+".pyclone.tsv").c_str());
      outpy << "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\tvariant_case\tvariant_freq\tgenotype\n";
    }

  out.open((filename+"_muts_expAF.txt").c_str());
  out << "Mut_ID\tChr\tPosition\tWt\tMut\tAF_obs\tCoverage\tAF_exp\tMut_Copies\tMut_Copies_Raw\tIs_Subclonal_CN\tiCN\tP_Is_Clonal\n";

  //additional allele-specific output
  if (seqcn.as)
  {   
      outas.open((filename+"_muts_expAF_as.txt").c_str());
      outas << "Mut_ID\tChr\tPosition\tWt\tMut\tAF_obs\tCoverage\tAF_exp\tMut_Copies\tMut_Copies_Raw\tIs_Subclonal_CN\tiCN\tP_Is_Clonal\t";
      outas << "CN_A\tCN_B\ttau1\tCN_Clone1_A\tCN_Clone1_B\ttau2\tCN_Clone2_A\tCN_Clone2_B\n";
  }
  for(i=0;i<vcf.size();i++)
    {
      if(vcf[i].index >= 0 && vcf[i].index < cn.size())
	{
	  ind = vcf[i].index;
	  cnA = (double)cn[ind].A;
	  cnB = (double)cn[ind].B;
	  cn1A = -1; 
	  cn1B = -1;
	  cn2A = -1; 
	  cn2B = -1;
	  tau1 = -1;
	  tau2 = -1;
	  
	  icn = cnA+cnB;
	  if(icn>0)
	    {
	      if(cn[ind].is_subclonal==0 || cn[ind].inconsistent_state==1)
		{
		  m_raw=vcf[i].AF*(2*(1-p)+p*icn)/p;
		  m=floor(m_raw+0.5);
		  if(m<=1)
		    m=1;
		  if(m>=cnA)
		    m=cnA;

		  af=m*p/(2*(1-p)+p*icn);
		}
	      else
		{
		  cnA= -1;
		  cnB= -1;
		  cn1A = (double)cn[ind].clone1_A;
		  cn1B = (double)cn[ind].clone1_B;
		  tau1 = cn[ind].clone1_tau;
		  cn2A = (double)cn[ind].clone2_A;
		  cn2B = (double)cn[ind].clone2_B;
		  tau2 = cn[ind].clone2_tau;
		  icn = (cn1A+cn1B)*tau1;
		  icn+= (cn2A+cn2B)*tau2;
		  m_raw=vcf[i].AF*(2*(1-p)+p*icn)/p;
		  
		  if(cn1A+cn1B > cn2A+cn2B)
		    tau=tau1;
		  else
		    tau=tau2;
		 
		  double dist_m=1e10;
		  m=1;
		  for(long mm=1; mm <= LMIN(cn1A,cn2A); mm++)
		    {
		      if(fabs((double)mm - m_raw) < dist_m)
			{
			  m = (double)mm;
			  dist_m = fabs((double)mm - m_raw);
			}

		      if(fabs((double)mm+tau - m_raw) < dist_m)
			{
			  m = (double)mm+tau;
			  dist_m = fabs((double)mm+tau - m_raw);
			}
		    }
		  
		  af=m*p/(2*(1-p)+p*icn);
		  
		}
	      
	      ccf=vcf[i].AF/af;
	      sig=sqrt(ccf*(1-vcf[i].AF)/(vcf[i].DP*af));
	      if(sig==0)
		pval=1;
	      else
		pval=0.5*(1+erf((ccf-1)/(sqrt(2)*sig)));

	      if(vcf[i].ref_seq.size()==1 && vcf[i].alt_seq.size()==1)
		{
		  if(seqcn.phylowgs==1)
		    {
		      outp << "s" << s << "\t" << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_SNM" << "\t";
		      outp <<  (long)(vcf[i].DP*(1-vcf[i].AF)+0.1) << "\t" << vcf[i].DP << "\t0.999\t" << 1-af << endl;
		      s++;
		    }

		  if(seqcn.pyclone==1)
		    {
		      outpy << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_SNM" << "\t";
		      outpy << (long)(vcf[i].DP*(1-vcf[i].AF)+0.1) << "\t" << vcf[i].DP-(long)(vcf[i].DP*(1-vcf[i].AF)+0.1) << "\t";
		      outpy << "2\t" << cn[ind].B << "\t" << cn[ind].A << "\t" << filename << "\t" << vcf[i].AF << "\t";
		      if(cn[ind].B==0)
			outpy << "BB\n";
		      else
			outpy << "AB\n";	
		    }

		  out << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_SNM" << "\t";
		  out << vcf[i].chr << "\t" << vcf[i].pos << "\t" << vcf[i].ref_seq << "\t";
		  out << vcf[i].alt_seq << "\t" << vcf[i].AF << "\t" << vcf[i].DP << "\t";
		  out << af << "\t" << m  << "\t" << m_raw << "\t" << cn[ind].is_subclonal << "\t";
		  out << icn << "\t" << pval << endl;
		
		  if (seqcn.as)
		  {
			outas << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_SNM" << "\t";
                  	outas << vcf[i].chr << "\t" << vcf[i].pos << "\t" << vcf[i].ref_seq << "\t";
                  	outas << vcf[i].alt_seq << "\t" << vcf[i].AF << "\t" << vcf[i].DP << "\t";
                  	outas << af << "\t" << m  << "\t" << m_raw << "\t" << cn[ind].is_subclonal << "\t";
                  	outas << icn << "\t" << pval << "\t" << cnA << "\t" << cnB << "\t";
                  	outas << tau1 << "\t" << cn1A << "\t" << cn1B << "\t";
                  	outas << tau2 << "\t" << cn2A << "\t" << cn2B << endl;
		  }

		}
	      else if(vcf[i].ref_seq.size()>1 && vcf[i].alt_seq.size()==1)
		{
		  out << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_DEL" << "\t";
		  out << vcf[i].chr << "\t" << vcf[i].pos << "\t" << vcf[i].ref_seq << "\t";
		  out << vcf[i].alt_seq << "\t" << vcf[i].AF << "\t" << vcf[i].DP << "\t";
		  out << af << "\t" << m  << "\t" << m_raw << "\t" << cn[ind].is_subclonal << "\t";
		  out << icn << "\t" << pval << endl;
		
		  if (seqcn.as)
		  {
		    outas << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_DEL" << "\t";
                    outas << vcf[i].chr << "\t" << vcf[i].pos << "\t" << vcf[i].ref_seq << "\t";
                    outas << vcf[i].alt_seq << "\t" << vcf[i].AF << "\t" << vcf[i].DP << "\t";
                    outas << af << "\t" << m  << "\t" << m_raw << "\t" << cn[ind].is_subclonal << "\t";
                    outas << icn << "\t" << pval << "\t" << cnA << "\t" << cnB << "\t";
                    outas << tau1 << "\t" << cn1A << "\t" << cn1B << "\t";
                    outas << tau2 << "\t" << cn2A << "\t" << cn2B << endl;
		  }
		}
	      else if(vcf[i].ref_seq.size()==1 && vcf[i].alt_seq.size()>1)
		{
		  out << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_INS" << "\t";
		  out << vcf[i].chr << "\t" << vcf[i].pos << "\t" << vcf[i].ref_seq << "\t";
		  out << vcf[i].alt_seq << "\t" << vcf[i].AF << "\t" << vcf[i].DP << "\t";
		  out << af << "\t" << m  << "\t" << m_raw << "\t" << cn[ind].is_subclonal << "\t";
		  out << icn << "\t" << pval << endl; 
		
		  if (seqcn.as)
		  {
		    outas << filename << "_" << vcf[i].chr << ":" << vcf[i].pos << "_INS" << "\t"; 
                    outas << vcf[i].chr << "\t" << vcf[i].pos << "\t" << vcf[i].ref_seq << "\t";
                    outas << vcf[i].alt_seq << "\t" << vcf[i].AF << "\t" << vcf[i].DP << "\t";
                    outas << af << "\t" << m  << "\t" << m_raw << "\t" << cn[ind].is_subclonal << "\t";
                    outas << icn << "\t" << pval << "\t" << cnA << "\t" << cnB << "\t";
                    outas << tau1 << "\t" << cn1A << "\t" << cn1B << "\t";
                    outas << tau2 << "\t" << cn2A << "\t" << cn2B << endl;
		  }
		}
	    }
	}
    }

  if(seqcn.phylowgs==1)
    outp.close();
  if(seqcn.pyclone==1)
    outpy.close();
  out.close();
  if (seqcn.as)
    outas.close();
}



void plot_pdf_output(double purity,double ploidy,double p_est,vector<double> chr_size,vector<p_data> obs_p,vector<seqcn_profile> profile,seqcn_par seqcn)
{
   //#### start plotting
   
   //PLOT the whole stuff as pdf
   //file pointer for interactive R
  long i;
  FILE *fp;
  long n_p = obs_p.size();

#ifdef RPLOTTING
  fp=popen("R --slave","w"); // open R
  

  //initialize R stream  
  fprintf(fp,"source(\"%s\")\n",((string)INSTALLDIR+"/R/plot_cn.R").c_str());
  fprintf(fp,"pdf(\"%s_cn_profile.pdf\",width=8,height=5.5)\n",seqcn.out_name.c_str());
  fprintf(fp,"d=numeric(0)\ncluster=numeric(0)\ndata=numeric(0)\n");
  
  fprintf(fp,"plot.profile(\"%s\",%f,\"%s\")\n",seqcn.out_name.c_str(),p_est,(seqcn.out_name+" ("+seqcn.mode+")").c_str());
  
  for(i=0;i<n_p;i++)
    {
      if(obs_p[i].iCN <= 5)
	fprintf(fp,"d = rbind(d,c(%f,%f,%f,%f,%f,%ld,%ld,%ld,%ld,%f,%f))\n",obs_p[i].iCN,obs_p[i].theta_t,obs_p[i].theta_sig_t,obs_p[i].theta_n,obs_p[i].theta_sig_n,obs_p[i].cluster,obs_p[i].n_snps,(long)obs_p[i].A,(long)obs_p[i].B,obs_p[i].theta_exp,(double)(obs_p[i].end - obs_p[i].start)/chr_size[chrom2number(obs_p[i].chr)]);
    }
  
  
  fprintf(fp,"snps.alleles(d,%f,%f,\"%s\")\n",purity,ploidy,seqcn.out_name.c_str());
  fprintf(fp,"plot.cn.profile(\"%s\",%f,%f)\n",seqcn.out_name.c_str(),purity,ploidy);
  
  //close R stream
  fprintf(fp,"a=dev.off()\n");
  pclose(fp);
#endif

   //#### end plotting
  system(("rm -f "+seqcn.out_name+"_cn_profile.txt").c_str());
}

void clear_old_seqcn_files(seqcn_par seqcn)
{
  system(("rm -f "+seqcn.out_name+"_allelic_states.txt").c_str());
  system(("rm -f "+seqcn.out_name+"_cluster_assignments.txt").c_str());
  system(("rm -f "+seqcn.out_name+"_uncorr_cn.seg").c_str());
  system(("rm -f "+seqcn.out_name+"_iCN.seg").c_str());
  system(("rm -f "+seqcn.out_name+"_mclusters.txt").c_str());
  system(("rm -f "+seqcn.out_name+"_muts_expAF.txt").c_str());
  system(("rm -f "+seqcn.out_name+"_subclonal_cn.txt").c_str());
  system(("rm -f "+seqcn.out_name+"_mcluster.pdf").c_str());
}

void write_seqcn_output(double purity,double &ploidy,double scale,char sex,vector<CN_data> CN,vector<p_data> out_cn,vector<vcf_data> vcf,vector<double> &chr_size,seqcn_par seqcn)
{
  long i;
  double het = 0; double norm_het = 0;
  long iCN;
  ofstream out;
  double incon;
  double norm_ploidy = 0;
  double logR;

     //compute fraction of subclonal copy number
   het = 0; norm_het=0;
   for(i=0;i<out_cn.size();i++)
     {
       if(out_cn[i].chr!="chrX" && out_cn[i].chr!="chrY" && out_cn[i].inconsistent_state==0
	  && out_cn[i].n_snps >=5 && out_cn[i].theta_exp >= 0)
	 {
	   norm_het += (double)(out_cn[i].end - out_cn[i].start);
	   if(out_cn[i].is_subclonal==1)
	     {
	       het += (double)(out_cn[i].end - out_cn[i].start);
	     }
	 }
     }
   het *= 1.0/norm_het;
   //compute ploidy
   ploidy = 0; norm_ploidy=0;
   incon = 0;
   chr_size=vector<double>(NCHR,0);
   
   for(i=0;i<out_cn.size();i++)
     {
       iCN = (long)((out_cn[i].cn-2*(1-purity))*scale+0.5);
       if(out_cn[i].chr!="chrX" && out_cn[i].chr!="chrY" && iCN >=0)
	 {
	   if(out_cn[i].inconsistent_state == 1)
	     incon +=  (double)(out_cn[i].end - out_cn[i].start);
	   ploidy += (double)iCN*((double)(out_cn[i].end - out_cn[i].start));
	   norm_ploidy += (double)(out_cn[i].end - out_cn[i].start);
	 }
       chr_size[chrom2number(out_cn[i].chr)]+=(double)(out_cn[i].end - out_cn[i].start);
     }
   //out.close();
   
   ploidy *= 1.0/norm_ploidy;
   incon *= 1.0/norm_ploidy;
   //pupl output to stdout
   out.open((seqcn.out_name+"_cn_summary.txt").c_str());
   out << "sample_name\tpurity\tploidy\tfraction_subclonal_cn\tsex_estimated\tstatus\tfraction_inconsistent_segs\n";
   out << seqcn.out_name << "\t" << purity << "\t" << ploidy << "\t" << het << "\t" << sex  << "\t" << seqcn.mode  << "\t" << incon << endl;
   out.close();

   out.open((seqcn.out_name+"_allelic_states.txt").c_str());
   out << "Sample\tChromosome\tStart\tEnd\tCopy_Nr_Raw\tCopyNr\tA\tB\tLOH\tTheta\tTheta_Exp\tn_SNPs\tIs_Subclonal_CN\tSubclonal_P_value\tIs_Inconsistent_State\n";
   for(i=0;i<out_cn.size();i++)
     {
       if(out_cn[i].n_snps >=1 && out_cn[i].theta_exp >= 0 && out_cn[i].chr != "chrX")
	 {
	   out << seqcn.out_name << "\t" << out_cn[i].chr << "\t" << out_cn[i].start << "\t" << out_cn[i].end << "\t";
	   out << (out_cn[i].cn-2*(1-purity))*scale << "\t";
	   out << out_cn[i].iCN << "\t" << out_cn[i].A << "\t" << out_cn[i].B << "\t";
	   if(out_cn[i].B == 0)
	     out << "1\t";
	   else
	     out << "0\t";
	   out << out_cn[i].theta_t << "\t" << out_cn[i].theta_exp << "\t" << out_cn[i].n_snps << "\t";
	   out << out_cn[i].is_subclonal << "\t" <<  out_cn[i].p_subclonal << "\t";
	   out << out_cn[i].inconsistent_state << endl;
	 }
     }
   out.close();


   //output subclonal copy numbers
    out.open((seqcn.out_name+"_subclonal_cn.txt").c_str());
    out << "Sample\tChromosome\tStart\tEnd\tSubclonal_CN\tClone1_A\tClone1_B\tClone1_Fraction\tClone2_A\tClone2_B\tClone2_Fraction\n";
    for(i=0;i<out_cn.size();i++)
      {
	if(out_cn[i].is_subclonal==1 && out_cn[i].inconsistent_state==0 && (out_cn[i].chr!="chrX" && out_cn[i].chr!="chrY"))
	  {
	    out << seqcn.out_name << "\t" << out_cn[i].chr << "\t" << out_cn[i].start << "\t" << out_cn[i].end << "\t";
	    out << (double)(out_cn[i].clone1_A+out_cn[i].clone1_B)*out_cn[i].clone1_tau+(double)(out_cn[i].clone2_A+out_cn[i].clone2_B)*out_cn[i].clone2_tau << "\t";
	    out << out_cn[i].clone1_A << "\t" << out_cn[i].clone1_B << "\t" << out_cn[i].clone1_tau << "\t";
	    out << out_cn[i].clone2_A << "\t" << out_cn[i].clone2_B << "\t" << out_cn[i].clone2_tau << "\n";
	  }
      }
    out.close();

   ofstream out_iCN;
   out_iCN.open((seqcn.out_name+"_iCN.seg").c_str());
   out.open((seqcn.out_name+"_uncorr_cn.seg").c_str());
   
   string chr_pre;
   chr_pre = CN[0].chr;
   double cn_pre = CN[0].cn_aggregated;
   long win;
   
   out << seqcn.out_name << "\t" << chr_pre << "\t" << CN[0].start;
   out_iCN << seqcn.out_name << "\t" << chr_pre << "\t" << CN[0].start;
   win = 0;
   for(i=0;i<CN.size();i++)
     {
      if(chr_pre != CN[i].chr || cn_pre != CN[i].cn_aggregated)
	{
	  out << "\t" << CN[i-1].end << "\t" << win << "\t" << cn_pre << endl;
	  if(CN[i-1].chr!="chrX" && CN[i-1].chr!="chrY")
	    out_iCN << "\t" << CN[i-1].end << "\t" << win << "\t" << LMAX((long)((cn_pre-2*(1-purity))*scale+0.5),0) << endl;
	  win = 0;
	  out << seqcn.out_name << "\t" << CN[i].chr << "\t" << CN[i].start;
	  if(CN[i].chr!="chrX" && CN[i].chr!="chrY")
	    out_iCN << seqcn.out_name << "\t" << CN[i].chr << "\t" << CN[i].start;
	}
      chr_pre = CN[i].chr; cn_pre = CN[i].cn_aggregated; win++;
     }
   out << "\t" << CN[CN.size()-1].end << "\t" << win << "\t" << cn_pre << endl;
   if(CN[CN.size()-1].chr!="chrX" && CN[CN.size()-1].chr!="chrY")
     out_iCN  << "\t" << CN[CN.size()-1].end << "\t" << win << "\t" <<  LMAX((long)((cn_pre-2*(1-purity))*scale+0.5),0) << endl;
   
   out.close();
   out_iCN.close();   

   //compute expected AF
   if(vcf.size()>0)
     {
       mut2cn(out_cn,vcf);
       vcf_expected_AF(seqcn.out_name,out_cn,purity,vcf,seqcn);
     }
}

//#### further functions 

double median(vector<double> &data)
{
  long i,j,n=data.size();
  vector<double> tmp;
  double median;

  tmp.clear();
  for(i=0;i<n;i++)
    tmp.push_back(data[i]);
  
  sort(tmp.begin(),tmp.end());
  j=tmp.size();
  if(j-2*(long)(j/2)==0)
    {
      median=0.5*(tmp[j/2-1]+tmp[j/2]);
    }
  else
    {
      median=tmp[j/2];
    }
  return(median);
}

//note that the actual bandwidth of the smoother is 2*w+1
void median_smooth(vector<double> &data,long w)
{
  long i,j,k,l;
  long n = data.size();
  vector<double> sd(n,0),med_v(2*w+1,0);

  for(i=0;i<n;i++)
    {
      l=0;
      for(j=-w;j<=w;j++)
	{
	  k = i+j;
	  //boundary treatment: copy end points over
	  if(k < 0)
	    k = (-1)*k;
	  if(k >= n)
	    k = n-2-(k-n);
	  //copy to med_v for median computation
	  med_v[l] = data[k];
	  l++;
	}
      sd[i] = median(med_v);
    }
  
  //replace input data by median smoothed values
  data = sd;
}

void segment(vector<double> dat,vector<double> sig,vector<long> &segStart,vector<long> &segEnd,double alpha,long minWd)
{
  long i,n,bkp;
  vector<long> tmpStart,tmpEnd;
  list<long> finalBKP;
  double score,maxScore,S1,S2,n1,n2,SQ1,SQ2,v1,v2;
  n=dat.size();
  
  for(i=0;i<segStart.size();i++)
    {
      tmpStart.push_back(segStart[i]);
      tmpEnd.push_back(segEnd[i]);
      finalBKP.push_back(segStart[i]);
      finalBKP.push_back(segEnd[i]);
    }
   
  //clear seg-vectors
  segStart.clear(); segEnd.clear();
  //main loop
  while(tmpStart.size()!=0)
    {
      S1=0; S2=0;
      SQ1=0; SQ2=0;
      n1=0; n2=0;
      v1=0; v2=0;
      maxScore=-1;
      //initialize sub segment
      for(i=tmpStart[0];i<=tmpEnd[0];i++)
	{
	  if(i-tmpStart[0]<minWd-1)
	    {
	      S1+=dat[i];
	      SQ1+=dat[i]*dat[i];
	      v1 += sig[i];
	      n1+=1;
	    }
	  else
	    {
	      S2+=dat[i];
	      SQ2+=dat[i]*dat[i];
	      v2 += sig[i];
	      n2+=1;
	    }
	}
      for(i=0;i<=(tmpEnd[0]-tmpStart[0]+1)-2*minWd;i++)
	{
	  //compute t-statistic for all sub-segments
	  S1+=dat[tmpStart[0]+minWd+i];
	  SQ1+=dat[tmpStart[0]+minWd+i]*dat[tmpStart[0]+minWd+i];
	  S2-=dat[tmpStart[0]+minWd+i];
	  SQ2-=dat[tmpStart[0]+minWd+i]*dat[tmpStart[0]+minWd+i];
	  n1+=1;
	  n2-=1;
	  //v1=SQ1/n1-S1*S1/(n1*n1);
	  //v2=SQ2/n2-S2*S2/(n2*n2);
	  v1+=sig[tmpStart[0]+minWd+i];
	  v2-=sig[tmpStart[0]+minWd+i];
	  score=fabs(S1/n1-S2/n2)/sqrt(v1/n1+v2/n2);
	  if(score>maxScore)
	    {
	      maxScore=score;
	      bkp=tmpStart[0]+minWd+i;
	    }
	}
      //is maxScore significantly high?
      if(erfc(maxScore/sqrt(2)) <=alpha && tmpEnd[0]-tmpStart[0]+1 >=minWd)
	{
	  finalBKP.push_back(bkp);
	  //subdivide segment 
	  tmpStart.push_back(tmpStart[0]);
	  tmpEnd.push_back(bkp);
	  tmpStart.push_back(bkp+1);
	  tmpEnd.push_back(tmpEnd[0]);
	  maxScore = 0;
	}
      tmpStart.erase(tmpStart.begin());
      tmpEnd.erase(tmpEnd.begin());
    }
  //sort break points
  finalBKP.sort();
  //get segments from break points
  list<long>::iterator it;
  for(it=finalBKP.begin();it!=finalBKP.end();it++)
    {
      if(*it!=n-1)
	{
	  if(*it==0)
	    {
	      segStart.push_back(0);
	      it++;
	      segEnd.push_back(*it);
	      it--;
	    }
	  else
	    {
	      segStart.push_back(*it+1);
	      it++;
	      segEnd.push_back(*it);
	      it--;
	    }
	}

    }
}

void read_vcf_line(string tmp,vcf_data &vcf_record)
{
  stringstream line1,line2;
  string tmp1,tmp2;
  vcf_record.DB=0;

  line1.str(""); line1.clear();
  line1 << tmp;
  getline(line1,tmp1,'\t');
  vcf_record.chr=tmp1;
  getline(line1,tmp1,'\t');
  vcf_record.pos=atol(tmp1.c_str());
  getline(line1,tmp1,'\t');
  vcf_record.ID=tmp1;
  getline(line1,tmp1,'\t');
  vcf_record.ref_seq=tmp1;
  getline(line1,tmp1,'\t');
  vcf_record.alt_seq=tmp1;
  getline(line1,tmp1,'\t');
  vcf_record.QUAL=atoi(tmp1.c_str());
  getline(line1,tmp1,'\t');
  vcf_record.FILTER=tmp1;
  getline(line1,tmp1,'\t');
  line1.str(""); line1.clear();
  line1 << tmp1;
  while(getline(line1,tmp1,';'))
    {
      line2.str(""); line2.clear();
      line2 << tmp1;
      tmp1=""; tmp2="";
      getline(line2,tmp1,'=');
      getline(line2,tmp2,'=');
      if(tmp1=="DP")
	vcf_record.DP=atol(tmp2.c_str());
      if(tmp1=="AF")
	vcf_record.AF=atof(tmp2.c_str());
     if(tmp1=="DP_N")
	vcf_record.DP_N=atol(tmp2.c_str());
      if(tmp1=="AF_N")
	vcf_record.AF_N=atof(tmp2.c_str());
      if(tmp1=="FR")
	vcf_record.FR=atof(tmp2.c_str());
      if(tmp1=="TG")
	vcf_record.TG=tmp2;
      if(tmp1=="DB")
	vcf_record.DB=1;
    }

  if(vcf_record.AF>1)
    vcf_record.AF=1;
  if(vcf_record.AF_N>1)
    vcf_record.AF_N=1;
}
