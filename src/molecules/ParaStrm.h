/** \file parastrm.h
  * \brief Parameter stream
  * \author Martin Bernreuther et al. <bernreuther@hlrs.de> 
  * \version 0.2
  * \date 08.02.2008
  * license: GPL
*/

/***************************************************************************
 *   Copyright (C) 2008 by Martin Bernreuther and colleagues               *
 *   bernreuther@hlrs.de                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   Due to copyleft all future versions of this program must be           *
 *   published as Free Software.                                           *
 ***************************************************************************/
#ifndef PARASTRM_H_
#define PARASTRM_H_

//#include <vector>

#include <cstdlib>
#include <cassert>

/**
Parameter stream

@author Martin Bernreuther
*/
class ParaStrm{
public:
  //enum Eparatype {UNDEF=0,DOUBLE,INT};

  ParaStrm() : m_size(0)/*, m_pos(0)*/ { m_pstrm = NULL; m_readpos = NULL; }
  ~ParaStrm() { if (m_pstrm) free(m_pstrm); }

  //std::size_t numpara() const { return m_ptypes.size(); }
  /// end of stream reached?
  bool eos() const { return m_readpos>=m_pstrm+m_size; }
  /// stream empty?
  bool empty() const { return m_size==0; }
  /// reset reading "pointer" to the beginning of the stream
  void reset_read() { m_readpos=m_pstrm; /*m_pos=0;*/ }

  /// add value to the (end of) stream
  template <class T> ParaStrm& operator <<(T p)
  {
    ptrdiff_t oldsize=m_size;
    // determine relative reading position
    ptrdiff_t readpos=m_readpos-m_pstrm;
    // enlarge vector
    m_size+=sizeof(p);
    if(!m_pstrm)
      m_pstrm=(char*)malloc(m_size);
    else
      m_pstrm=(char*)realloc(m_pstrm,m_size);
    // realloc might change m_pstrm!
    assert(m_pstrm);
    // save value
    *((T*)(m_pstrm+oldsize))=p;
    //m_ptypes.push_back(UNDEF);
    // adjust m_readpos to maybe changed m_pstrm
    m_readpos=m_pstrm+readpos;
    return *this;
  }
  //ParaStrm& operator <<(double p);
  //ParaStrm& operator <<(int p);

  /// read value at actual stream position and move position
  template <class T> ParaStrm& operator >>(T& p)
  {
    assert(m_readpos+sizeof(T)<=m_pstrm+m_size);
    // get value
    p=*((T*)m_readpos);
    // move reading position
    m_readpos=(char*)((T*)m_readpos+1);
    //++m_pos;
    return *this;
  }

private:
  size_t m_size;
  char *m_pstrm, *m_readpos;  // sizeof(char)==1!
  // we also could check for valid types...
  //std::vector<Eparatype> m_ptypes;
  //std::size_t m_pos;
};

#endif /*PARASTRM_H_*/
