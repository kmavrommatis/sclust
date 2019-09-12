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

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include "samtools/bam.h"
#include "samtools/bam_endian.h"
#include "samtools/kstring.h"
#include "samtools/sam_header.h"

#include "samtools/sam.h"
#include "Sclust.h"

using namespace std;

char bam_flag2char_table_peiflyne[] = "pPuUrR12sfd\0\0\0\0\0";


inline int bam2nt(int in)
{
  if(in == 1)
    return(0);
  else if(in == 2)
    return(1);
  else if(in == 4)
    return(2);
  else if(in == 8)
    return(3);
  else
    return(4);

}

bool read_bam(samfile_t *in,bam1_t *b,bam_map &m)
{

  if(samread(in,b) >= 0)
    {
      bam_header_t *header;
      int of; char cig;
      stringstream line;

      m.seq.clear();
      m.mask_subst.clear();

      header = in->header;
      of = in->type>>2&3;

      uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
      int i,j,k,l;
      const bam1_core_t *c = &b->core;
      
      m.qname = (string)bam1_qname(b);
      //flag
      if (of == BAM_OFDEC)
	{
	  m.flag = (long)c->flag;
	}
      else if (of == BAM_OFHEX)
	{
	  cerr << "BAM_OFHEX not supported; please contact developer.\n";
	  exit(1);
	}
     else
	{
	  cerr << "BAM_OFSTR not supported; please contact developer.\n";
	  exit(1);
	}

      if (c->tid < 0)
	m.rname = "*";
      else
	m.rname = header->target_name[c->tid];

      m.pos = c->pos +1;
      m.mapq = c->qual;
      
      m.n_cigar = c->n_cigar;
      m.is_MID = 1;
      if (c->n_cigar != 0)
	{
	  m.is_indel=0;
	  for (i = 0; i < c->n_cigar; ++i) 
	    {
	      m.cigar_pos[i]=(bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT);
	      cig = ("MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	      m.cigar_type[i]=cig;
	      if(cig == 'I' || cig == 'D')
		m.is_indel=1;
	      if(cig != 'I' && cig != 'D' && cig != 'M') 
		m.is_MID = 0;
	    }
	}
      else
	m.is_MID = 0;
      
      //mrnm
      if (c->mtid < 0)
	m.mrnm="*";
      else if (c->mtid == c->tid)
	m.mrnm="=";
      else
	{
	  if (header)
	    m.mrnm=header->target_name[c->mtid];
	  else
	    {
	      line.str(""); line.clear();
	      line << c->mtid;
	      m.mrnm=line.str();
	    }
	}
      
      m.mpos=c->mpos + 1;
      m.isize=c->isize;
      
      if (c->l_qseq) 
	{
	  for (i = 0; i < c->l_qseq; ++i)
	    m.seq.push_back(bam2nt(bam1_seqi(s, i)));
	  // //create substitution mask
	  // for (i = 0; i < MAXRL; ++i)
	  //   m.mask_subst[i] = -1;
 
	  j=0; l=0;
	  for(i = 0; i< m.n_cigar;i++)
	    {
	      if(m.cigar_type[i] == 'M')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      m.mask_subst.push_back(l);
		      j++; l++; 
		    }
		}
	      else if(m.cigar_type[i] == 'S' || m.cigar_type[i] == 'H')
		{
		  if(m.cigar_type[i] == 'S')
		    {
		      for(k=0;k<m.cigar_pos[i];k++)
			{
			  if(i!=0)
			    {
			      m.mask_subst.push_back(-1);
			      j++; 
			    }
			  l++; 
			}
		    }
		  else
		    m.mapq=0;
		}
	      else if(m.cigar_type[i] == 'D')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      m.mask_subst.push_back(-1);
		      j++;
		    }
		}
	      else if(m.cigar_type[i] == 'N')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      m.mask_subst.push_back(-1);
		      j++;
		    }
		}
	      else if(m.cigar_type[i] == 'I')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      l++;
		    }
		}
	      else
		{
		  //cerr << "Error: P not yet implemented ---> " << m.qname << "\n";
		  m.mapq = 0;
		}
	    }
	  m.mask_size = m.mask_subst.size();
	  if (t[0] == 0xff)
	    m.qual="*";
	  else
	    {
	      line.str(""); line.clear();
	      for (i = 0; i < c->l_qseq; ++i)
		line << (char)(t[i] + 33);
	      m.qual=line.str();
	    }
	}
      else
	{
	  m.qual="*";
	}
      m.size = c->l_qseq;
      return(1);
    }
  else
    return(0);
}



bool read_bam_iter(samfile_t *in, bam_iter_t iter,bam1_t *b,bam_map &m)
{

  if(bam_iter_read(in->x.bam, iter, b) >= 0)
    {
      bam_header_t *header;
      int of; char cig;
      stringstream line;

      m.seq.clear();
      m.mask_subst.clear();

      header = in->header;
      of = in->type>>2&3;

      uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
      int i,j,k,l;
      const bam1_core_t *c = &b->core;
      
      m.qname = (string)bam1_qname(b);
      //flag
      if (of == BAM_OFDEC)
	{
	  m.flag = (long)c->flag;
	}
      else if (of == BAM_OFHEX)
	{
	  cerr << "BAM_OFHEX not supported; please contact developer.\n";
	  exit(1);
	}
     else
	{
	  cerr << "BAM_OFSTR not supported; please contact developer.\n";
	  exit(1);
	}

      if (c->tid < 0)
	m.rname = "*";
      else
	m.rname = header->target_name[c->tid];

      m.pos = c->pos +1;
      m.mapq = c->qual;
      
      m.n_cigar = c->n_cigar;
      m.is_MID = 1;
      if (c->n_cigar != 0)
	{
	  m.is_indel=0;
	  for (i = 0; i < c->n_cigar; ++i) 
	    {
	      m.cigar_pos[i]=(bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT);
	      cig = ("MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	      m.cigar_type[i]=cig;
	      if(cig == 'I' || cig == 'D')
		m.is_indel=1;
	      if(cig != 'I' && cig != 'D' && cig != 'M') 
		m.is_MID = 0;
	    }
	}
      else
	m.is_MID = 0;
      
      //mrnm
      if (c->mtid < 0)
	m.mrnm="*";
      else if (c->mtid == c->tid)
	m.mrnm="=";
      else
	{
	  if (header)
	    m.mrnm=header->target_name[c->mtid];
	  else
	    {
	      line.str(""); line.clear();
	      line << c->mtid;
	      m.mrnm=line.str();
	    }
	}
      
      m.mpos=c->mpos + 1;
      m.isize=c->isize;
      
      if (c->l_qseq) 
	{
	  for (i = 0; i < c->l_qseq; ++i)
	    m.seq.push_back(bam2nt(bam1_seqi(s, i)));
	  // //create substitution mask
	  // for (i = 0; i < MAXRL; ++i)
	  //   m.mask_subst[i] = -1;
 
	  j=0; l=0;
	  for(i = 0; i< m.n_cigar;i++)
	    {
	      if(m.cigar_type[i] == 'M')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      m.mask_subst.push_back(l);
		      j++; l++; 
		    }
		}
	      else if(m.cigar_type[i] == 'S' || m.cigar_type[i] == 'H')
		{
		  if(m.cigar_type[i] == 'S')
		    {
		      for(k=0;k<m.cigar_pos[i];k++)
			{
			  if(i!=0)
			    {
			      m.mask_subst.push_back(-1);
			      j++; 
			    }
			  l++; 
			}
		    }
		  else
		    m.mapq=0;
		}
	      else if(m.cigar_type[i] == 'D')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      m.mask_subst.push_back(-1);
		      j++;
		    }
		}
	      else if(m.cigar_type[i] == 'N')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      m.mask_subst.push_back(-1);
		      j++;
		    }
		}
	      else if(m.cigar_type[i] == 'I')
		{
		  for(k=0;k<m.cigar_pos[i];k++)
		    {
		      l++;
		    }
		}
	      else
		{
		  //cerr << "Error: P not yet implemented ---> " << m.qname << "\n";
		  m.mapq = 0;
		}
	    }
	  m.mask_size = m.mask_subst.size();
	  if (t[0] == 0xff)
	    m.qual="*";
	  else
	    {
	      line.str(""); line.clear();
	      for (i = 0; i < c->l_qseq; ++i)
		line << (char)(t[i] + 33);
	      m.qual=line.str();
	    }
	}
      else
	{
	  m.qual="*";
	}
      m.size = c->l_qseq;
      return(1);
    }
  else
    return(0);
}
