#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "samtools/sam.h"

using namespace std;

#define SCLUSTVERSION "1.1"
#define NCHR 24
#define MAXCIGAR 1024
#define MAXRL 1024

inline void print_header()
{
  string header = "\n\t --- Sclust: Mutation Clustering using Smoothing Splines ---\n \
      \t\t    Version "+(string)SCLUSTVERSION+"  21/March/2019\n \
      \t     University of Cologne \n \
      \t     Peifer Lab\n\n";
  cerr << header;
}


typedef struct {
  long target;
  long pos;
  string chr;
  double s_base_t[4];
  double s_base_n[4];
  int allele_A;
  int allele_B;
} SNP_data;


typedef struct {
  string chr;
  long start;
  long end;
  double cn;
  double cn_coarse;
  double cn_aggregated;
  double sig_cn;
  long n_reads_t;
  long n_reads_n;
  double GC_content;
  double loc_bg;
} CN_data;

typedef struct{
  double min_reads;
  double lambda;
  double w;
  double alpha;
  long median_smooth;
  long min_seg;
  string out_name;
  long ns;
  double st;
  double max_p;
  double min_p;
  double min_pu;
  double pu_mu;
  bool failed;
  bool inv_likelihood;
  bool f2;
  string mode;
  double ploidy;
  double L_bi;
  double L_cn;
  string sv;
  bool phylowgs;
  bool pyclone;
  bool as;
} seqcn_par;

typedef struct{
  string chr;
  long start;
  long end;
  double cn;
  long cluster;
  long n_snps;
  double theta_t;
  double theta_sig_t;
  double theta_n;
  double theta_sig_n;
  double A;
  double B;
  double iCN;
  double theta_exp;
  double theta_corr;
  double ls;
  vector<double> seg_af;
  vector<double> seg_af_sig;
  vector<double> seg_af_exp;
  vector<long> n_seg_af;
  vector<long> cn_index;
  bool is_subclonal;
  double p_subclonal;
  long clone1_A;
  long clone1_B;
  double clone1_tau;
  long clone2_A;
  long clone2_B;
  double clone2_tau;
  bool inconsistent_state;
} p_data;

typedef struct{
  double cn_mean;
  double cn_sig;
  double theta_mean;
  double m;
} p_cluster_data;

typedef struct{
  double scale;
  double ploidy_exp;
  double purity;
  double ploidy;
  double L_cn;
  double L_bi;
} seqcn_profile;



typedef struct {
  string chr;
  long pos;
  string ID;
  string ref_seq;
  string alt_seq;
  int QUAL;
  string FILTER;
  int NS;
  long DP;   //Coverage Tumor
  double AF; //Allelic Frequency Tumor
  long DP_N;   //Coverage Normal
  double AF_N; //Allelic Frequency Normal
  double FR;  //Forw-Rev score
  string TG; // Target Name
  long index;
  bool DB;
} vcf_data;

typedef struct {
  string qname;
  long flag;
  string rname;
  long pos;
  long mapq;
  int cigar_pos[MAXCIGAR];
  char cigar_type[MAXCIGAR];
  int n_cigar;
  string mrnm;
  long mpos;
  long isize;
  vector<int> seq;
  string qual;
  long size;
  vector<int> mask_subst;
  int mask_size;
  bool is_indel;
  bool is_MID;
} bam_map;

typedef struct
{
  string chr;
  long pos;
  char wt;
  bool in_partition;
  vector<double> bases_f;
  vector<double> bases_r;
  vector<double> bases_spliced;
  vector<string> insertion;
  map<string,long> insertion_f;
  vector<string> deletion;
  map<string,long> deletion_f;
  vector<string> ins;
  vector<long> ins_f;
  vector<string> del;
  vector<long> del_f;
  vector<long> n_ins;
  vector<long> n_del;
  double coverage;
} pileup;


//prototypes

void cn(int argc, char *argv[]);
void cluster(int argc, char *argv[]);
void bamprocess(int argc, char *argv[]);

int chrom2number(string chr);
double median(vector<double> &data);
void median_smooth(vector<double> &data,long w); 
void segment(vector<double> dat,vector<double> sig,vector<long> &segStart,vector<long> &segEnd,double alpha,long minWd); //in Sclust.cc
void read_vcf_line(string tmp,vcf_data &vcf_record); 
bool read_bam(samfile_t *in,bam1_t *b,bam_map &m);
bool read_bam_iter(samfile_t *in, bam_iter_t iter,bam1_t *b,bam_map &m);
long base2num(char base);
bool is_SNP(long pos,ifstream &in,string &SNP_id,long &pos_SNP,char &allele_A,char &allele_B);
