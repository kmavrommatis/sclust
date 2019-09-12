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
#include <zlib.h>
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
#include <map>
#include <list>
#include <sstream>

#include "Sclust.h"
#include "installdir.h"
#include "nibtools.h"
#include "samtools/sam.h"

using namespace std;


string presclust_help = "      SYNOPSIS \n \t Sclust bamprocess  <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -t        \t tumor alignment file (sorted, indexed bam-file)\n \
     \t -n        \t normal alignment file (sorted, indexed bam-file)\n \
     \t -i        \t input for merging chromosome files (prefix only)\n \
     \t -o        \t output name\n\
     \t -r        \t chromosome to extract\n \
     \t -build    \t genome build (possible: hg38, hg19, mm10) [hg19]\n\
     \t -part     \t choose genome partition [1]:\n";

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

#define LMAX(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg2) : (lminarg1))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

#define DMAX(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg2) : (dminarg1))


//prototypes
void presclust_core(string chromosome,string part_name,string out_name,string t_file,string n_file,string build,long max_insert, vector<SNP_data> &snp_list, vector<CN_data> &CN);
void pileup_buffer_bam_presclust(string rname,long pos,bam_map &m1,samfile_t *fp,bam1_t *b,vector<pileup> &buffer,bool &is_eof,long &cur_pos,long &n_reads,vector<long> &rs);
void merge(string in_name,string out_name,string build);


void bamprocess(int argc, char *argv[])
{
  int longindex,opt;
  // option definition

  static struct option longopts[]={
    {"help"     , 0, 0,  'h'},
    {"help"     , 0, 0,  '?'},
    {"t"        , 1, 0,    1},
    {"n"        , 1, 0,    2},
    {"o"        , 1, 0,    3},
    {"build"    , 1, 0,    4},
    {"part"     , 1, 0,    5},
    {"r"        , 1, 0,    6},
    {"i"        , 1, 0,    7},
    {0, 0, 0, 0}
  }; 
  
  string t_file="";
  string n_file="";
  string in_file="";
  string out_name="";
  string build="hg19";
  string chromosome="";
  ofstream out;
  ifstream in;
  long i; int list=0;
  string tmp;

  //read partition list file

  vector<string> part_list;
  
  in.open(((string)INSTALLDIR+"/../partitions/partition_list.txt").c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open partition table.\n";
      exit(1);
    }
  while(in >> tmp)
    part_list.push_back(tmp);
  in.close();

  optind=0;
  //parse command line arguments
  while((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
        {
        case 'h':
	  print_header();
	  cerr << presclust_help;
	  for(i=0;i<part_list.size();i++)
	    cerr << "                  \t   (" << i+1 << ") " << part_list[i] << endl;
	  cerr << endl;
	  exit(1);
          break;
        case '?':
	  print_header();
	  cerr << presclust_help;
	  for(i=0;i<part_list.size();i++)
	    cerr << "                  \t   (" << i+1 << ") " << part_list[i] << endl;
	  cerr << endl;
	  exit(1);  
          break;
        case 1:
	  t_file=(string)optarg;
          break;
        case 2:
	  n_file=(string)optarg;
          break;
        case 3:
	  out_name=(string)optarg;
          break;
	case 4:
	  build = (string)optarg;
	  if(build != "hg38" && build != "hg19" && build != "mm10")
	    {
	      cerr << "Error: please use hg38, hg19 or mm10. Falling back to default value: hg19.\n";
	      build = "hg19";	      
	    }
	  break;
	case 5:
	  list=abs(atoi(optarg))-1;
	  if(list < 0 || list >= part_list.size())
	    {
	      cerr << "Error: invalid genome partition.\n";
	      exit(1);
	    }
	  break;
	case 6:
	  chromosome=(string)optarg;
          break;
	case 7:
	  in_file=(string)optarg;
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
      cerr << presclust_help;
      for(i=0;i<part_list.size();i++)
	cerr << "                  \t   (" << i+1 << ") " << part_list[i] << endl;
      cerr << endl;
      exit(1);
    }

  if(in_file=="")
    {
      if(t_file=="" || out_name == "" || chromosome == "")
	{
	  cerr << "Error: please specify alignment files (tumor and normal), output name, or chromosome name.\n";
	  print_header();
	  cerr << presclust_help;
	  for(i=0;i<part_list.size();i++)
	    cerr << "                  \t   (" << i+1 << ") " << part_list[i] << endl;
	  cerr << endl;
	  exit(1);
	}
    }
  else
    {
      if(out_name == "")
	{
	  cerr << "Error: please specify output name for merging.\n";
	  print_header();
	  cerr << presclust_help;
	  for(i=0;i<part_list.size();i++)
	    cerr << "                  \t   (" << i+1 << ") " << part_list[i] << endl;
	  cerr << endl;
	  exit(1);
	}
    }

  //some variables
  long max_insert=1000;
  vector<SNP_data> snp_list;
  vector<CN_data> CN;


  if(in_file=="" && n_file != "")
    {
      presclust_core(chromosome,part_list[list],out_name,t_file,n_file,build,max_insert,snp_list,CN);
      
      //save intermediate data
      out.open((out_name+"_"+chromosome+"_bamprocess_data.txt").c_str());
      out << "<snp_list>\n"; 
      for(i=0;i<snp_list.size();i++)
	{
	  out << snp_list[i].target << "\t" << snp_list[i].chr << "\t" << snp_list[i].pos;
	  for(int ii=0;ii<4;ii++)
	    out << "\t" << snp_list[i].s_base_t[ii];
	  for(int ii=0;ii<4;ii++)
	    out << "\t" << snp_list[i].s_base_n[ii];
	  out << "\t" << snp_list[i].allele_A <<  "\t" << snp_list[i].allele_B << endl;
	}
      out << "</snp_list>\n"; 
      out << "<cn>\n"; 
      for(i=0;i<CN.size();i++)
	{
	  out << CN[i].chr << "\t" << CN[i].start << "\t" << CN[i].end << "\t" << CN[i].n_reads_t << "\t" << CN[i].n_reads_n;
	  out << "\t" << CN[i].GC_content << endl;
	}
      out << "</cn>\n";       
      out.close();
    }
  else
    {
      merge(in_file,out_name,build);
    }
  
}


void presclust_core(string chromosome,string part_name,string out_name,string t_file,string n_file,string build,long max_insert, vector<SNP_data> &snp_list, vector<CN_data> &CN)
{
  long max_buffer_size=10000;
  long min_coverage=15;
  long local_coverage=0;

  snp_list.clear(); CN.clear();
  samfile_t *fp_t,*fp_n;
  bam1_t *b_t=bam_init1();
  bam1_t *b_n=bam_init1();
  
  //open bam file stream
  fp_t = samopen(t_file.c_str(), "rb", 0);
  fp_n = samopen(n_file.c_str(), "rb", 0);
  
  if(fp_t==NULL || fp_n==NULL)
    {
      cerr << "Error: cannot open tumor/normal bam file.\n";
      exit(1);
    }

  //open index stuff
  bam_index_t *idx_t = 0,*idx_n=0;
  idx_t = bam_index_load(t_file.c_str());
  idx_n = bam_index_load(n_file.c_str());

  if(idx_t == 0 || idx_n == 0)
    {
      cerr << "Error: please index bam-files.\n";
      exit(1);
    }

  //parse region
  int tid_t, beg_t, end_t;
  int tid_n, beg_n, end_n;
  bam_parse_region(fp_t->header,chromosome.c_str(), &tid_t, &beg_t, &end_t);
  bam_parse_region(fp_n->header,chromosome.c_str(), &tid_n, &beg_n, &end_n);
  if(tid_t < 0 || tid_n < 0)
    {
      cerr << "Error: chromosome "<< chromosome << " does not exist.\n";
      exit(1);
    }
  //create bam iterator
  bam_iter_t iter_t,iter_n;
  iter_t = bam_iter_query(idx_t, tid_t, beg_t, end_t);
  iter_n = bam_iter_query(idx_n, tid_n, beg_n, end_n);

  // variables for the targetList
  vector< vector<string> > Tname(NCHR,vector<string>(0,""));
  vector< vector<long> >   Tstart(NCHR,vector<long>(0,0));
  vector< vector<long> >   Tend(NCHR,vector<long>(0,0));
  vector<double> avg_cover_t,avg_cover_n;
  double tmp_cover_t=0,tmp_cover_n=0;
  string chr,name; long sta,sto; long curr_target;
  long global_target=0; 
  
  //variables for SNP's
  string SNP_id;
  long pos_SNP;
  char allele_A,allele_B;
  bool SNP_at_pos;
  SNP_data tmp_snp;

  //some variables
  long i,j,k,percentage,n_bases;
  long n_ref,ref_pos;
  string ref_name,tmp;
  string header;
  char wt_base;
  bool is_eof_t=0,is_eof_n=0,hit_last=0;
  ifstream in,inSNP;
  long t_start,t_end,dt,total_t_chr=0,total_t=0;
  CN_data tmp_cn;
  long tmp_n_reads_t,tmp_n_reads_n;
  long reads_target_t,reads_target_n;
  long n_reads_t=0,n_reads_n=0;
  vector<long> rs_t,rs_n;
  long ii,jj,kk;
  double min_e=1;
  stringstream line;
  map<string,string> key_ref;
  map<string,string>::iterator it_ref;
  map<string,long>::iterator it;
  
  //variables bam format
  bam_map m1_t,m1_n,m1_tmp;
  long max_buffer_pos_t, max_buffer_pos_n;
  vector<bam_map> m1_buffer_t,m1_buffer_n;
  long cur_pos_t,cur_pos_n;

  //variables for variant calling
  vector<double> s_t(5,0),s_n(5,0);

  //indel stuff
  long n_indel;
  string indel_pre;

  //nib access
  nib nibObj; int stat;
  
  //output pileup format
  ofstream out_pileup_t;
  ofstream out_pileup_n;
  bool write_out;
  double s_other;

  //create header
  line.str(""); line.clear();
  line << "#Sclust " << SCLUSTVERSION << endl;
  header = line.str();

  //read reference list
  in.open(((string)INSTALLDIR+"/../annotation/"+build+"/ref_names.txt").c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open reference names file.\n";
      exit(1);
    } 

  while(in >> tmp)
    key_ref[tmp]=tmp;
  in.close();

  //read target list
  in.open(((string)INSTALLDIR+"/../partitions/"+build+"_"+part_name+".txt").c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open partition list.\n";
      exit(1);
    }

  while(in >> chr >> sta >> sto >> name)
    {
      i=chrom2number(chr);
      Tstart[i].push_back(sta);
      Tend[i].push_back(sto);
      Tname[i].push_back(name);
    }
  in.close();

  //initialize m1_t and m1_n
  is_eof_t = read_bam_iter(fp_t,iter_t,b_t,m1_t);
  is_eof_n = read_bam_iter(fp_n,iter_n,b_n,m1_n);
  
  int rr=0;
  while(is_eof_t && is_eof_n) //loop over chromosomes
    {
      if(rr!=0)
	cout << "\nTotal time used to process " << ref_name << ":\t" << total_t_chr/60 << " min " << total_t_chr-60*(total_t_chr/60) << " sec  " << endl << endl;
      rr=1;

      if(is_eof_t==0 || is_eof_n==0 || m1_t.rname!=chromosome || m1_n.rname!=chromosome)
	{
	  hit_last=1;
	  break;
	}


      if(m1_t.rname!=m1_n.rname)
	{
	  cerr << "Error: different reference names between tumor and norma\n";
	  exit(1);
	}
      
      ref_name=m1_t.rname;
      k=chrom2number(ref_name);
      //check if chromosome is in traget list
      if(Tend[k].size()==0)
	exit(1);

      cout << "Processing " << ref_name << ":\n";
      total_t_chr=0;
      stat=nibObj.open((string)INSTALLDIR+"/../annotation/"+build+"/"+build+"_"+ref_name+".nib");
      if(stat!=0)
	{
	  cerr << "Error: " << errormsg[stat] << " " << (string)INSTALLDIR+"/"+build+"/"+build+"_"+ref_name+".nib" << endl;
	  exit(1);
	}
      inSNP.open(((string)INSTALLDIR+"/../annotation/SNP_"+build+"/"+build+"_snp_"+ref_name+".txt").c_str());
      if(!inSNP.is_open())
	{
	  cerr << "Error: cannot open " << ((string)INSTALLDIR+"/SNP_"+build+"/"+build+"_snp_"+ref_name+".txt.") << endl;
	  exit(1);
 	}
      //some initializations
      percentage=10; n_bases=nibObj.size();
      t_start=time(NULL);
      pos_SNP=-1;
      max_buffer_pos_t=0; max_buffer_pos_n=0;
      curr_target=0;
      m1_tmp.rname=ref_name; m1_tmp.pos = -10000;
      m1_tmp.size = 0;
      cur_pos_t=-1; cur_pos_n=-1;
      m1_buffer_t =vector<bam_map> (max_buffer_size,m1_tmp);
      m1_buffer_n =vector<bam_map> (max_buffer_size,m1_tmp);
      
      //pileup buffer initialization
      pileup pileup_tmp;
      pileup_tmp.bases_f = vector<double> (5,0);
      pileup_tmp.bases_r = vector<double> (5,0);
      pileup_tmp.bases_spliced = vector<double> (5,0);
      vector<pileup> p_buffer_t(max_buffer_size,pileup_tmp);
      vector<pileup> p_buffer_n(max_buffer_size,pileup_tmp);
      n_reads_t = 0; n_reads_n = 0;
      reads_target_t=0; reads_target_n=0;
      
      vector<char> ref_chromosome; ref_chromosome.push_back('N'); //to occupy the zero position
      cout << "Read reference of " << ref_name << " into memory.\n";
      for(ref_pos=1;ref_pos <= n_bases;ref_pos++)
	{
	  nibObj.nextBase(&wt_base);
	  ref_chromosome.push_back(wt_base);
	}
      cout << "Processing the data: ";
      for(ref_pos=1;ref_pos <= n_bases;ref_pos++) //loop over position
	{
	  if(100*(double)ref_pos/(double)n_bases>=percentage)
	    {
	      t_end=time(NULL); dt=t_end-t_start;
	      total_t+=dt; total_t_chr+=dt;
	      cout << "*";
	      t_start=time(NULL);
	      cout.flush();
	      percentage+=10;
	    }
	  
	  wt_base=ref_chromosome[ref_pos];
	  
	  SNP_at_pos=is_SNP(ref_pos,inSNP,SNP_id,pos_SNP,allele_A,allele_B);
	  if(ref_pos>Tend[k][curr_target] && curr_target < Tend[k].size())
	    {
	      
	      //write cn_data
	      tmp_cn.chr=ref_name;
	      tmp_cn.start=Tstart[k][curr_target];
	      tmp_cn.end=Tend[k][curr_target];
	      tmp_cn.n_reads_t = n_reads_t;
	      tmp_cn.n_reads_n = n_reads_n;
	      //compute GC-data
	      long p_size=0; tmp_cn.GC_content=0;
	      for(long kk=Tstart[k][curr_target];kk<=Tend[k][curr_target];kk++)
		{
		  if(ref_chromosome[kk]!='N')
		    p_size++;
		  if(ref_chromosome[kk]=='G' || ref_chromosome[kk]=='C')
		    tmp_cn.GC_content += 1.0;
		}

	      if(p_size==0)
		tmp_cn.GC_content = -1;
	      else
		tmp_cn.GC_content*= 1.0/(double)p_size;

	      CN.push_back(tmp_cn);

	      n_reads_t = 0; n_reads_n = 0;
	      reads_target_t=0; reads_target_n=0;
	      tmp_cover_t=0; tmp_cover_n=0;
	      curr_target++; global_target++;
	    }

	  pileup_buffer_bam_presclust(ref_name,ref_pos,m1_t,fp_t,b_t,p_buffer_t,is_eof_t,cur_pos_t,tmp_n_reads_t,rs_t);
	  pileup_buffer_bam_presclust(ref_name,ref_pos,m1_n,fp_n,b_n,p_buffer_n,is_eof_n,cur_pos_n,tmp_n_reads_n,rs_n);	   
	  

	   if(ref_pos >= Tstart[k][curr_target] && ref_pos <= Tend[k][curr_target] && wt_base!='N')
	    {
	      n_reads_t += tmp_n_reads_t;
	      n_reads_n += tmp_n_reads_n;
	    }

	  //compute coverage for actual position
	  p_buffer_t[cur_pos_t].coverage = 0;
	  p_buffer_n[cur_pos_t].coverage = 0;
	  ii=base2num(wt_base);
	  s_other = 0;
	  for(kk=0;kk<4;kk++)
	    {
	      p_buffer_t[cur_pos_t].coverage += p_buffer_t[cur_pos_t].bases_f[kk];
	      p_buffer_t[cur_pos_t].coverage += p_buffer_t[cur_pos_t].bases_r[kk];
	      if(kk!=ii)
		{
		  s_other += p_buffer_t[cur_pos_t].bases_f[kk];
		  s_other += p_buffer_t[cur_pos_t].bases_r[kk];
		}
	      p_buffer_n[cur_pos_t].coverage += p_buffer_n[cur_pos_t].bases_f[kk];
	      p_buffer_n[cur_pos_t].coverage += p_buffer_n[cur_pos_t].bases_r[kk];
	    }
	  
	  if(ref_pos >= Tstart[k][curr_target] && ref_pos <= Tend[k][curr_target] && wt_base!='N' && p_buffer_t[cur_pos_t].coverage >= 5 && p_buffer_n[cur_pos_t].coverage >= 5)
	    {
	      tmp_cover_t+=p_buffer_t[cur_pos_t].coverage;
	      tmp_cover_n+=p_buffer_n[cur_pos_t].coverage;
	      if(SNP_at_pos)
	  	{
	  	  i=base2num(allele_A); j=base2num(allele_B);
	  	  tmp_snp.target=global_target;
	  	  tmp_snp.pos=ref_pos;
	  	  tmp_snp.chr=ref_name;
	  	  tmp_snp.allele_A=i;
	  	  tmp_snp.allele_B=j;
	  	  for(ii=0;ii<4;ii++)
	  	    {
	  	      tmp_snp.s_base_t[ii]=p_buffer_t[cur_pos_t].bases_f[ii]+p_buffer_t[cur_pos_t].bases_r[ii];
	  	      tmp_snp.s_base_n[ii]=p_buffer_n[cur_pos_t].bases_f[ii]+p_buffer_n[cur_pos_t].bases_r[ii];
	  	    }
	  	  snp_list.push_back(tmp_snp);
		}
	    }
	} 
      cout << endl;
      //close nib stream
      nibObj.close();
      inSNP.close();
    }

  //close bam stream
  samclose(fp_t); samclose(fp_n);

  cout << "\nProcessing of " << ref_name << " successfully accomplished. Total time used: " << total_t/60 << " min " << total_t-60*(total_t/60) << " sec  " << endl << endl;
}


void pileup_buffer_bam_presclust(string rname,long pos,bam_map &m1,samfile_t *fp,bam1_t *b,
			      vector<pileup> &buffer,bool &is_eof,long &cur_pos,long &n_reads,vector<long> &rs)
{
  long i,j,l,k,pos_buffer_read;
  string tmp;
  double eps;
  long q_cutoff=10;
  bool not_pair;
  long max_buffer_size = 10000;
  long n_N;
  n_reads = 0;
  rs.clear();
  map<string,long>::iterator it;

  //clear last position
  if(cur_pos-1 < 0)
    {
      buffer[max_buffer_size-1].bases_f[0] = 0; buffer[max_buffer_size-1].bases_f[1] = 0;
      buffer[max_buffer_size-1].bases_f[2] = 0; buffer[max_buffer_size-1].bases_f[3] = 0;
      buffer[max_buffer_size-1].bases_r[0] = 0; buffer[max_buffer_size-1].bases_r[1] = 0;
      buffer[max_buffer_size-1].bases_r[2] = 0; buffer[max_buffer_size-1].bases_r[3] = 0;
      buffer[max_buffer_size-1].bases_spliced[0] = 0; buffer[max_buffer_size-1].bases_spliced[1] = 0;
      buffer[max_buffer_size-1].bases_spliced[2] = 0; buffer[max_buffer_size-1].bases_spliced[3] = 0;
      buffer[max_buffer_size-1].insertion.clear();
      buffer[max_buffer_size-1].deletion.clear();
      buffer[max_buffer_size-1].insertion_f.clear();
      buffer[max_buffer_size-1].deletion_f.clear();
      buffer[max_buffer_size-1].ins.clear();
      buffer[max_buffer_size-1].del.clear();
      buffer[max_buffer_size-1].n_ins.clear();
      buffer[max_buffer_size-1].n_del.clear();
    }
  else
    {
      buffer[cur_pos-1].bases_f[0] = 0; buffer[cur_pos-1].bases_f[1] = 0;
      buffer[cur_pos-1].bases_f[2] = 0; buffer[cur_pos-1].bases_f[3] = 0;
      buffer[cur_pos-1].bases_r[0] = 0; buffer[cur_pos-1].bases_r[1] = 0;
      buffer[cur_pos-1].bases_r[2] = 0; buffer[cur_pos-1].bases_r[3] = 0;
      buffer[max_buffer_size-1].bases_spliced[0] = 0; buffer[max_buffer_size-1].bases_spliced[1] = 0;
      buffer[max_buffer_size-1].bases_spliced[2] = 0; buffer[max_buffer_size-1].bases_spliced[3] = 0;
      buffer[cur_pos-1].insertion.clear();
      buffer[cur_pos-1].deletion.clear();
      buffer[cur_pos-1].insertion_f.clear();
      buffer[cur_pos-1].deletion_f.clear();
      buffer[cur_pos-1].ins.clear();
      buffer[cur_pos-1].del.clear();
      buffer[cur_pos-1].n_ins.clear();
      buffer[cur_pos-1].n_del.clear();
    }
  
  //move buffer forward
  cur_pos++;
  if(cur_pos == max_buffer_size)
    cur_pos = 0;

  //center to current position if needed
  while(rname == m1.rname && is_eof == 1 && m1.pos < pos)
    {
      //read reads
      is_eof = read_bam(fp,b,m1);
    }
  //exract all reads from current position and write it to buffer
  while(rname == m1.rname && is_eof == 1 && m1.pos == pos)
    {
      //count number of reads starting at actual position
      if(m1.mapq >0)
	{
	  n_reads++;
	  rs.push_back(m1.seq.size());
	}
      
      n_N = 0;
      for(i = 0; i< m1.n_cigar;i++)
	{
	  if(m1.cigar_type[i] == 'N')
	    {
	      n_N += m1.cigar_pos[i];
	    }
	}
      //buffer indels
      if(m1.mapq >= 10 && m1.is_indel == 1 && (m1.flag & 0x400)/0x400==0)
	{
	  j=cur_pos; l=0;
	  for(i = 0; i< m1.n_cigar;i++)
	    {
	      if(m1.cigar_type[i] == 'N')
		{
		  n_N += m1.cigar_pos[i];
		}
	      if(m1.cigar_type[i] == 'M')
		{
		  for(k=0;k<m1.cigar_pos[i];k++)
		    {
		      j++; l++; 
		      if(j == max_buffer_size)
			j=0;
		    }
		}
	      else if(m1.cigar_type[i] == 'S' || m1.cigar_type[i] == 'H')
		{
		  for(k=0;k<m1.cigar_pos[i];k++)
		    {
		      if(i!=0)
			{
			  j++;
			  if(j == max_buffer_size)
			    j=0;
			} 
		      l++; 
		    }
		}
	      else if(m1.cigar_type[i] == 'D')
		{
		  tmp="";
		  for(k=0;k<m1.cigar_pos[i];k++)
		    tmp+="-";
		  buffer[j].deletion.push_back(tmp);

		  it=buffer[j].deletion_f.find(tmp);
		  if(it==buffer[j].deletion_f.end())
		    buffer[j].deletion_f[tmp] = (long)((m1.flag & 0x10)/0x10 ==0);
		  else
		    it->second += (long)((m1.flag & 0x10)/0x10 ==0);

		  for(k=0;k<m1.cigar_pos[i];k++)
		    {
		      j++;
		      if(j == max_buffer_size)
			j=0;
		    }
		}
	      else if(m1.cigar_type[i] == 'I')
		{
		  tmp="";
		  for(k=0;k<m1.cigar_pos[i];k++)
		    tmp+="ACGTN"[m1.seq[l+k]];
		  buffer[j].insertion.push_back(tmp);

		  it=buffer[j].insertion_f.find(tmp);
		  if(it==buffer[j].insertion_f.end())
		    buffer[j].insertion_f[tmp] = (long)((m1.flag & 0x10)/0x10 ==0);
		  else
		    it->second += (long)((m1.flag & 0x10)/0x10 ==0);

		  for(k=0;k<m1.cigar_pos[i];k++)
		    l++;
		}
	      else
		{
		  cerr << "Error: P not yet implemented.\n";
		  m1.mapq=0;
		}
	    } 
	}
      
      //write single bases to buffer
      pos_buffer_read = cur_pos;
      for(i=0;i<LMIN(max_buffer_size-1,m1.mask_size);i++)
	{
	  if(m1.mask_subst[i] < 0)
	    pos_buffer_read++;
	  else
	    {
	      //buffer substitutions
	      if((int)(m1.qual[m1.mask_subst[i]]-33) >= q_cutoff && m1.mapq > 0)
		{
		  eps = 1-pow(10,-0.1*(double)m1.mapq); //mapping quality
		  eps*= 1-pow(10,-0.1*(double)((int)m1.qual[m1.mask_subst[i]]-33)); //base phred score

		  //remove PCR-duplicates
		  if((m1.flag & 0x400)/0x400==1)
		    eps=0;

		  if(eps >= 0.99)
		    {
		      if((m1.flag & 0x10)/0x10 ==0)
			buffer[pos_buffer_read].bases_f[m1.seq[m1.mask_subst[i]]]+= 1;
		      else
			buffer[pos_buffer_read].bases_r[m1.seq[m1.mask_subst[i]]]+= 1;
		      if(n_N >10)
			buffer[pos_buffer_read].bases_spliced[m1.seq[m1.mask_subst[i]]]+= 1;
		    }
		}
	      pos_buffer_read++;
	    }
	  if(pos_buffer_read == max_buffer_size)
	    pos_buffer_read = 0;
	}
      //read next read
      is_eof = read_bam(fp,b,m1);
    }
  
}


void merge(string in_name,string out_name,string build)
{
  string tmp;
  stringstream line;
  ifstream in;
  ofstream out;
  long i,j,ii;
  long n_chr;
  vector<string> chr_names;
  bool is_single = 1;
  string header="#Sclust "+(string)SCLUSTVERSION;


  //read chromosome list
  in.open(((string)INSTALLDIR+"/../annotation/"+build+"/ref_names.txt").c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open chromosome list.\n";
      exit(1);
    }
 
  while(in >> tmp)
    {
      chr_names.push_back(tmp);
    }
  in.close();
  n_chr = chr_names.size();
  
  SNP_data tmp_snp;
  CN_data  tmp_cn;
  long sel=0,sum_cn=0;
  vector< vector<SNP_data> > snp_list(n_chr,vector<SNP_data>(0,tmp_snp));
  vector< vector<CN_data> > cn(n_chr,vector<CN_data>(0,tmp_cn));
  double sum_t=0;
  long sum_n=0;
  
  for(i=0;i<n_chr;i++)
    {
      in.open((in_name+"_"+chr_names[i]+"_bamprocess_data.txt").c_str());
      if(in.is_open())
	{
	  while(getline(in,tmp,'\n'))
	    {
	      if(tmp=="<snp_list>")
		sel = 1;
	      else if(tmp=="<cn>")
		sel = 2;
	      else if(tmp.substr(0,2)=="</")
		sel = 0;
	      else
		{
		  line.str(""); line.clear();
		  line << tmp;
		  //parse snp list
		  if(sel==1)
		    {
		      line >> tmp_snp.target >> tmp_snp.chr >> tmp_snp.pos;
		      line >> tmp_snp.s_base_t[0] >> tmp_snp.s_base_t[1] >> tmp_snp.s_base_t[2];
		      line >> tmp_snp.s_base_t[3] >> tmp_snp.s_base_n[0] >> tmp_snp.s_base_n[1];
		      line >> tmp_snp.s_base_n[2] >> tmp_snp.s_base_n[3] >> tmp_snp.allele_A  >> tmp_snp.allele_B;
		      snp_list[i].push_back(tmp_snp);
		    }
		  if(sel==2)
		    {
		      line >> tmp_cn.chr >> tmp_cn.start >> tmp_cn.end;
		      line >> tmp_cn.n_reads_t >> tmp_cn.n_reads_n >> tmp_cn.GC_content;
		      cn[i].push_back(tmp_cn);
		    }
		}
	    }
	  in.close();		
	}
    }
      
  //write cn data
  out.open((out_name+"_rcount.txt").c_str());
  out << header << endl;
  out << "chromosome\tstart\tend\trcount_t\trcount_n\tGC\n";
  
  for(i=0;i<n_chr;i++)
    {
      for(j=0;j<cn[i].size();j++)
	{ 
	  out << cn[i][j].chr << "\t" << cn[i][j].start << "\t" << cn[i][j].end  << "\t" << cn[i][j].n_reads_t << "\t" << cn[i][j].n_reads_n << "\t" << cn[i][j].GC_content << endl;
	}
    }
  out.close();
  
  //write snp list
  out.open((out_name+"_snps.txt").c_str());
  out << header << endl;
  out << "part_no\tchromosome\tposition\tT_A\tT_C\tT_G\tT_T\tN_A\tN_C\tN_G\tN_T\tallele_A\tallele_B\n";
  
  for(i=0;i<n_chr;i++)
    {
      for(j=0;j<snp_list[i].size();j++)
	{
	  out << snp_list[i][j].target+sum_cn << "\t" << snp_list[i][j].chr << "\t" << snp_list[i][j].pos;
	  for(int ii=0;ii<4;ii++)
	    out << "\t" << snp_list[i][j].s_base_t[ii];
	  for(int ii=0;ii<4;ii++)
	    out << "\t" << snp_list[i][j].s_base_n[ii];
	  out << "\t" << snp_list[i][j].allele_A <<  "\t" << snp_list[i][j].allele_B << endl;
	}
      sum_cn += cn[i].size();
    }
  out.close(); 
}
