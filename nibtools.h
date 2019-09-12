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

#include<iostream>
#include<fstream>

using namespace std;

#define MSK 0x6be93d3a

//error messages
const string errormsg[]={"","wrong format"
			 ,"cannot open file"
			 ,"file is not open"
			 ,"position beyond sequence boundary"};

class nib
{
 private:
  long nBases;  //sequence length
  ifstream in; 
  char raw[4]; //binary input
  unsigned long tmp;
  int base_low;
  int base_high;
  bool is_high;

  int bin2ascii(char *out,int inp) //converts binary to ascii
  {
    switch(inp&0xff)
      {
      case 0:  *out='T'; break;
      case 1:  *out='C'; break;
      case 2:  *out='A'; break;
      case 3:  *out='G'; break;
      case 4:  *out='N'; break;
      case 8:  *out='T'; break;
      case 9:  *out='C'; break;
      case 10: *out='A'; break;
      case 11: *out='G'; break;
      default: *out='N'; return 1;
      }
    return 0;
  };
  unsigned char ascii2bin(char inp)
  {
    switch(inp)
      {
      case 'T': return 0; break;
      case 'C': return 1; break;
      case 'A': return 2; break;
      case 'G': return 3; break;
      case 'N': return 4; break;
      case 't': return 8; break;
      case 'c': return 9; break;
      case 'a': return 10; break;
      case 'g': return 11; break;
      default: return 4;
      }
  };
 public:
  int open(string filename);
  void close()
    {
      if(!in.is_open())
	in.close();
    }
  int getBase(char *base,unsigned long pos);
  int nextBase(char *base);
  int write(string seq,string filename);
  unsigned long size()
  {
    if(!in.is_open())
      return 0;
    else
      return nBases;
  }
};
