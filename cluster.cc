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
#include <ctime>

#include "installdir.h"
#include "Sclust.h"
#include "numerics/nr.h"
#include "qp++.h"

using namespace std;

string mcluster_help = "      SYNOPSIS \n \t Sclust cluster <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -i      \t input data filename (prefix)\n \
     \t -o      \t output name (prefix) \n \
     \t -nbins  \t number of histogram bins [100]\n \
     \t -indel  \t include indels for clustering\n \
     \t -lambda \t smoothing parameter [1e-7]\n\n";


typedef struct{
  string mut_id;
  string chr;
  long pos;
  string wt;
  string mut;
  double af_obs;
  long coverage;
  double af_exp;
  long cluster_id;
  double cluster_prob;
  vector<double> p;
} mcluster_data;

//prototypes
void read_mcluster_data(string file,bool use_indels,vector<mcluster_data> &mcluster);
void eval_spline(long n,double x,prec_t *xd,prec_t *g,prec_t *gam,double &y,double &dy);
void find_ccf_clusters(long n,prec_t *x,prec_t *g,prec_t *gam,double sd_noise,vector<double> &cluster_ccf,vector<double> &peak_height);
void assign_muts2clusters(vector<mcluster_data> &mcluster,vector<double> cluster_ccf);

//main program
void cluster(int argc,char *argv[])
{

  int longindex,opt;
  // option definition

  static struct option longopts[]={
    {"help"  , 0, 0,  'h'},
    {"help"  , 0, 0,  '?'},
    {"i"     , 1, 0,    1},
    {"o"     , 1, 0,    2},
    {"nbins" , 1, 0,    3},
    {"lambda", 1, 0,    4},
    {"indel" , 0, 0,    5},
    {0, 0, 0, 0}
  }; 

  string in_name="";
  string out_name="";
  bool use_indels=0;
  long nbins=100;
  double alpha=1e-7;
  double max_x=1.5;
  //parse command line arguments
  while((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
        {
        case 'h':
	  print_header();
	  cerr <<mcluster_help;
	  exit(1);
          break;
        case '?':
	  print_header();
	  cerr <<mcluster_help;
	  cerr << endl;
	  exit(1);  
          break;
        case 1:
	  in_name=(string)optarg;
          break;
        case 2:
	  out_name=(string)optarg;
          break;
	case 3:
	  nbins=atol(optarg);
	  break;
	case 4:
	  alpha=atof(optarg);
	  break;
	case 5:
	  use_indels=1;
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
      cerr <<mcluster_help;
      cerr << endl;
      exit(1);
    }

  if(in_name=="")
    {
      cerr << "Error: please specify input name.\n";
      print_header();
      cerr << mcluster_help;
      exit(1);
    }

  if(out_name=="")
    out_name=in_name;
 
  long i,j;
  string tmp;
  ifstream in;
  ofstream out;
  double dtmp;
  double sd_mean=0,sd_sd=0;
  vector<mcluster_data> mcluster;
  double dx=max_x/(double)nbins,sum=0;
  double sd_noise=0;
  double ccf;
 
  //histogram data
  prec_t *x_l,*x_r,*x,*y,*w;
  x_l=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );
  x_r=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );
  x=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );
  y=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );
  w=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );

  vinit(x_l,nbins);
  vinit(x_r,nbins);
  vinit(x,nbins);
  vinit(y,nbins);
  vector<double> ccf_raw;
  //spline data
  prec_t *g,*gam;
  g=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );
  gam=(prec_t *)mkl_malloc( nbins*sizeof( double ), 64 );

  vinit(g,nbins);
  vinit(gam,nbins);
  vector<double> sd_vect;
  //file pointer for interactive R (plotting)
  FILE *fp;

  //read input data
  read_mcluster_data(in_name,use_indels,mcluster);

  //create histogram
  for(i=0;i<nbins;i++)
    {
      x_l[i] = ((double)i)*dx;
      x_r[i] = ((double)i+1.0)*dx;
      x[i] = ((double)i+0.5)*dx;
    }

  sum=0;
  for(i=0;i<mcluster.size();i++)
    {
      ccf=mcluster[i].af_obs/mcluster[i].af_exp;
      j=(long)(ccf/dx);
      if(j >= 0 && j < nbins)
	{
	  dtmp=sqrt((mcluster[i].af_obs*(1-mcluster[i].af_obs)/(double)mcluster[i].coverage))/mcluster[i].af_exp;
	  sd_vect.push_back(dtmp);
	  ccf_raw.push_back(ccf);
	  sd_mean+=dtmp;
	  sum+=1.0;
	}
    }

  sd_mean*=1.0/sum;
  sum=0;
  for(i=0;i<sd_vect.size();i++)
    {
      sd_sd+=(sd_vect[i]-sd_mean)*(sd_vect[i]-sd_mean);
      sum+=1.0;
    }
  sd_sd*=1.0/sum;

  sd_sd=sqrt(sd_sd);
  if(sd_vect.size()>0)
    sd_noise = median(sd_vect);
  else
    {
      cerr << "Error: no data for sd calculation\n";
      exit(1);
    }

  for(i=0;i<nbins;i++)
    w[i]=sd_noise;

  sum=0;
  for(i=0;i<ccf_raw.size();i++)
    {
      if(sd_vect[i] <= sd_noise+0.5*sd_sd)
	{
	  j=(long)(ccf_raw[i]/dx);
	  y[j]+=1.0;
	  sum+=1.0;
	}
    }

  for(i=0;i<nbins;i++)
    {
      y[i] *= 1.0/(sum*dx);
    }
  
  //solve deconvolution problem
  spline_deconv(x, nbins, y, nbins, w, g, gam, alpha);

  //peak finding
  vector<double> cluster_ccf;
  vector<double> peak_height;
  find_ccf_clusters(nbins,x,g,gam,sd_noise,cluster_ccf,peak_height);
  if(cluster_ccf.size()==0)
    {
      cout << "No cluster in valid CCF range found, please recalibrate copy number analysis.\n";
      exit(1);
    }
  assign_muts2clusters(mcluster,cluster_ccf);
  vector<long> n_muts_cluster(cluster_ccf.size(),0);
  //write cluster assignments
  out.open((out_name+"_cluster_assignments.txt").c_str());
  out << "Mut_ID\tChr\tPosition\tWt\tMut\tCCF\tCoverage\tCluster_Id\tCluster_CCF\tProbability";
  for(i=0;i<cluster_ccf.size();i++)
    out << "\t" << "P" << i;
  out << "\n";

  for(i=0;i<mcluster.size();i++)
    {
      n_muts_cluster[mcluster[i].cluster_id]++;
      out << mcluster[i].mut_id << "\t" << mcluster[i].chr << "\t" << mcluster[i].pos << "\t" << mcluster[i].wt << "\t";
      out << mcluster[i].mut << "\t" << mcluster[i].af_obs/mcluster[i].af_exp << "\t" << mcluster[i].coverage << "\t";
      out << mcluster[i].cluster_id << "\t" << cluster_ccf[mcluster[i].cluster_id] << "\t" << mcluster[i].cluster_prob;
      for(j=0;j<mcluster[i].p.size();j++)
	out << "\t" << mcluster[i].p[j];
      out << "\n"; 
    }
  out.close();
  //write cluster results
  out.open((out_name+"_mclusters.txt").c_str());
  if(!out.is_open())
    {
      cerr << "Error: cannot open " << out_name << "_mclusters.txt.\n";
      exit(1);
    }

  out << "Cluster_ID\tCCF_Cluster\tCluster_Peak_Height\tMutations_In_Cluster\n";
  for(i=0;i<cluster_ccf.size();i++)
    out << i << "\t" << cluster_ccf[i] << "\t" << peak_height[i] << "\t" << n_muts_cluster[i] << endl;
  out.close();

 #ifdef RPLOTTING
  //plot cluster results
  fp=popen("R --slave","w"); // open R

  tmp=(string)INSTALLDIR+"/R/plot_cluster.R";
  fprintf(fp,"source(\"%s\")\n",tmp.c_str());
  fprintf(fp,"pdf(\"%s_mcluster.pdf\",width=6,height=6,useDingbats=FALSE)\n",out_name.c_str());

  //x-values
  fprintf(fp,"x=c(%f",x[0]);
  for(i=1;i<nbins;i++)
    {
      fprintf(fp,",%f",x[i]);
    }
  fprintf(fp,")\n");

  //y-values
  fprintf(fp,"y=c(%f",x[0]);
  for(i=1;i<nbins;i++)
    {
      fprintf(fp,",%f",y[i]);
    }
  fprintf(fp,")\n");

  //cluster
  fprintf(fp,"cluster=numeric(0)\n");
  for(i=0;i<cluster_ccf.size();i++)
    {
      fprintf(fp,"cluster=rbind(cluster,c(%f,%f))\n",cluster_ccf[i],peak_height[i]);
    }

  //curve data
  double xx;
  double y_out,dy_out;
  fprintf(fp,"curve=numeric(0)\n");
  for(i=0;i<=1000;i++)
    {
      xx=x[0]+(double)i*(x[nbins-1]-x[0])/1000.0;
      eval_spline(nbins,xx,x,g,gam,y_out,dy_out);
      fprintf(fp,"curve=rbind(curve,c(%f,%f,%f))\n",xx,y_out,dy_out);
    }

  fprintf(fp,"plot.cluster(x,y,cluster,curve,\"%s\")\n",out_name.c_str());
  fprintf(fp,"m=dev.off()\n");
  pclose(fp);
#endif

  mkl_free(x_l);
  mkl_free(x_r);
  mkl_free(x);
  mkl_free(y);
  mkl_free(g);
  mkl_free(gam);
  mkl_free(w);
}

void find_ccf_clusters(long n,prec_t *x,prec_t *g,prec_t *gam,double sd_noise,vector<double> &cluster_ccf,vector<double> &peak_height)
{
  long i,k;
  vector<double> cl_start_s;
  vector<double> cl_end_s;
  vector<double> tmp_cluster_ccf;
  vector<double> y,dy;
  double tmp_y,tmp_dy;
  double max_y;
  long nc;
  long nsmp=1000;
  double dx;
  double xx;
  cluster_ccf.clear();
  peak_height.clear();
  for(i=0;i<n;i++)
    {
      eval_spline(n,x[i],x,g,gam,tmp_y,tmp_dy);
      y.push_back(tmp_y);
      dy.push_back(tmp_dy);
    }

  //scan for peak regions
  for(i=0;i<n-1;i++)
    {
      if(dy[i] > 0 && dy[i+1] <= 0)
	{
	  cl_start_s.push_back(x[i]);
	  cl_end_s.push_back(x[i+1]);
	}
    }
  nc=cl_start_s.size();
  tmp_cluster_ccf = vector<double> (nc,-1);
  //find exact peak location
  for(k=0;k<nc;k++)
    {
      max_y=-1;
      dx=(cl_end_s[k]-cl_start_s[k])/(double)nsmp;
      for(i=0;i<=nsmp;i++)
	{
	  xx=cl_start_s[k]+dx*(double)i;
	  eval_spline(n,xx,x,g,gam,tmp_y,tmp_dy);
	  if(tmp_y > max_y)
	    {
	      tmp_cluster_ccf[k]=xx;
	      max_y=tmp_y;
	    }
	}
      eval_spline(n,xx,x,g,gam,tmp_y,tmp_dy);
    }

  for(k=nc-1;k>=0;k--)
    {  
      if(tmp_cluster_ccf[k] < 1.0+sd_noise)
	{
	  cluster_ccf.push_back(tmp_cluster_ccf[k]);
	  eval_spline(n,tmp_cluster_ccf[k],x,g,gam,tmp_y,tmp_dy);
	  peak_height.push_back(tmp_y);
	}
    }
}

void eval_spline(long n,double x,prec_t *xd,prec_t *g,prec_t *gam,double &y,double &dy)
{
  long i,k;
  double b,dx;
  double eps = 0.01;

  i=0;
  for(k=0;k<n-1;k++)
    {
      if(x > xd[k] && x <=xd[k+1])
	i=k;
    }
  dx=xd[i+1]-xd[i];
  b=(1+(x-xd[i])/dx)*gam[i+1]+(1+(xd[i+1]-x)/dx)*gam[i];
  y=((x-xd[i])*g[i+1]+(xd[i+1]-x)*g[i])/dx-(1.0/6.0)*(x-xd[i])*(xd[i+1]-x)*b;

  if(y<=eps)
    {
      y=0; dy=0;
      return;
    }
  
  dy=(g[i+1]-g[i])/dx-(1.0/6.0)*(xd[i]+xd[i+1]-2*x)*b-(1.0/6.0)*(x-xd[i])*(xd[i+1]-x)*(gam[i+1]-gam[i])/dx;

}

void assign_muts2clusters(vector<mcluster_data> &mcluster,vector<double> cluster_ccf)
{
  long i,j;
  long n=mcluster.size();
  long nc=cluster_ccf.size();
  double sd,ccf,sum;
  double x;

  double max_p;
  long max_cl;
  
  for(i=0;i<n;i++)
    {
      mcluster[i].cluster_id=0;
      mcluster[i].cluster_prob=1;
    }
  if(nc > 1)
    {
      for(i=0;i<n;i++)
	{
	  vector<double> p(nc,0);
	  if(mcluster[i].af_obs==1)
	    continue;
	  ccf=mcluster[i].af_obs/mcluster[i].af_exp;
	  sd=(1.0/mcluster[i].af_exp)*sqrt(mcluster[i].af_obs*(1-mcluster[i].af_obs)/(double)mcluster[i].coverage);
	  sum=0;
	  for(j=0;j<nc;j++)
	    {
	      x=ccf-cluster_ccf[j];
	      if(j==0)
		{
		  p[j]=0.5+0.5*erf(x/(sqrt(2)*sd));
		}
	      else if(j==nc-1)
		{
		  p[j]=0.5+0.5*erf(-x/(sqrt(2)*sd));
		}
	      else
		{
		  p[j]=1-erf(fabs(x)/(sqrt(2)*sd));
		}
	      sum+=p[j];
	    }
	  max_p=-1; max_cl=0;
	  for(j=0;j<nc;j++)
	    {
	      p[j]=p[j]/sum;
	      if(max_p< p[j])
		{
		  max_p = p[j];
		  max_cl= j;
		}
	    }
	  mcluster[i].cluster_id=max_cl;
	  mcluster[i].cluster_prob=max_p;
	  mcluster[i].p=p;
	}

    }

}

void read_mcluster_data(string file,bool use_indels,vector<mcluster_data> &mcluster)
{
  long i;
  string tmp;
  stringstream line;
  mcluster_data tmp_mcluster;
  ifstream in;

  mcluster.clear();
  
  in.open((file+"_muts_expAF.txt").c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open: " << file << "_muts_expAF.txt\n";
      exit(1);
    }

  getline(in,tmp,'\n');
  while(getline(in,tmp,'\n'))
    {
      line.str(""); line.clear();
      line << tmp;

      getline(line,tmp,'\t');
      tmp_mcluster.mut_id=tmp;
      getline(line,tmp,'\t');
      tmp_mcluster.chr=tmp;
      getline(line,tmp,'\t');
      tmp_mcluster.pos=atol(tmp.c_str());
      getline(line,tmp,'\t');
      tmp_mcluster.wt=tmp;
      getline(line,tmp,'\t');
      tmp_mcluster.mut=tmp;
      getline(line,tmp,'\t');
      tmp_mcluster.af_obs=atof(tmp.c_str());
      getline(line,tmp,'\t');
      tmp_mcluster.coverage=atol(tmp.c_str());
      getline(line,tmp,'\t');
      tmp_mcluster.af_exp=atof(tmp.c_str());
      
      if(tmp_mcluster.af_obs>0)
	{
	  if(use_indels==0)
	    {
	      if(tmp_mcluster.mut_id.find("_SNM") < tmp_mcluster.mut_id.size())
		 mcluster.push_back(tmp_mcluster);
	    }
	  else
	    mcluster.push_back(tmp_mcluster);
	}
    }
  
}
