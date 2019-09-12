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


#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "Sclust.h"

using namespace std;

string commands = "      SYNOPSIS \n \t Sclust [command] <options>\n\n \
     COMMANDS:\n\n \
     \t bamprocess \t processing bam files to generate input data for Sclust\n \
     \t cn         \t compute copy number analysis prior mutation clustering\n \
     \t cluster    \t 1D mutation clustering based on smoothing spline\n\n";


int main(int argc,char *argv[])
{
  string command;

  if(argc==1)
    {
      print_header();
      cerr << commands;
      exit(1);
    }

  command=(string)argv[1];

  if(command=="cn")
    cn(argc,argv); // copy number analysis
  else if(command=="cluster")
    cluster(argc,argv); // mutation clustering
  else if(command=="bamprocess")
    bamprocess(argc,argv);

  else
    {
      print_header();
      cerr << commands;
      exit(1);
    }
}

int chrom2number(string chr)
{
  int i;
  if(chr=="chrX")
    i=23;
  else if(chr=="chrY")
    i=24;
  else
    i=atoi((chr.substr(chr.find("chr")+3)).c_str());
  return(i-1);
}


double meanCount(vector<long> count)
{
  long n=count.size(),i;
  double mean=0;
  
  if(n==0)
    return(0);
  
  for(i=0;i<n;i++)
    mean+=count[i];
  mean*=1.0/(double)n;

  return(mean);
}

long base2num(char base)
  {
    if(toupper(base)=='A')
      return(0);
    else if(toupper(base)=='C')
      return(1);
    else if(toupper(base)=='G')
      return(2);
    else if(toupper(base)=='T')
      return(3);
    else
      return(4);
  }

bool is_SNP(long pos,ifstream &in,string &SNP_id,long &pos_SNP,char &allele_A,char &allele_B)
{
  string tmp;
  if(!in.eof())
    {
      if(pos_SNP < pos)
	in >> SNP_id >> tmp >> pos_SNP >> allele_A >> allele_B;
    }
  
  if(pos==pos_SNP)
    return(1);
  else
    return(0);
    
}

