// $Id: grid.hpp 30 2008-04-04 12:39:42Z perroe $

#ifndef NRLIB_GRID_HPP
#define NRLIB_GRID_HPP

#include <vector>
#include <sstream>

#include "../exception/exception.hpp"

namespace NRLib2 {

template<class A>
class Grid {
public:
  typedef typename std::vector<A>::iterator iterator;
  typedef typename std::vector<A>::const_iterator const_iterator;
  typedef typename std::vector<A>::reference reference;
  typedef typename std::vector<A>::const_reference const_reference;

  Grid();
    /// \param val Initial cell value. 
  Grid(int ni, int nj, int nk, const A& val = A());
  virtual ~Grid();

    /// All values in the grid are erased when the grid is
    /// resized.
    /// \param val Initial cell value. 
  void Resize(int ni, int nj, int nk, const A& val = A());

    /// \throw IndexOutOfRange if one of the indexes was invalid.
  inline reference operator()(int i, int j, int k);
    /// \throw IndexOutOfRange if the index was invalid.
  inline reference operator()(int index);

    /// \throw IndexOutOfRange if one of the indexes was invalid.
  inline const_reference operator()(int i, int j, int k) const;
    /// \throw IndexOutOfRange if the index was invalid.
  inline const_reference operator()(int index) const;

  iterator begin() {return(data_.begin());}
  iterator end()   {return(data_.end());}

  const_iterator begin() const {return(data_.begin());}
  const_iterator end() const   {return(data_.end());}

  int GetNI() const {return(ni_);}
  int GetNJ() const {return(nj_);}
  int GetNK() const {return(nk_);}
  int GetN()  const {return(static_cast<int>(data_.size()));}

  inline int GetIndex(int i, int j, int k) const;
  void GetIJK(int index, int &i, int &j, int &k) const;

private:
  int ni_;
  int nj_;
  int nk_;

  std::vector<A> data_;
};

template<class A> 
Grid<A>::Grid() 
  : ni_(0),
    nj_(0),
    nk_(0),
    data_()
{}

template<class A> 
Grid<A>::Grid(int ni, int nj, int nk, const A& val) 
  : ni_(ni),
    nj_(nj),
    nk_(nk),
    data_(ni*nj*nk, val)
{}

template<class A>
Grid<A>::~Grid() 
{}

template<class A> 
void Grid<A>::Resize(int ni, int nj, int nk, const A& val)
{
  ni_ = ni;
  nj_ = nj;
  nk_ = nk;

  data_.resize(0); //To avoid copying of elements
  data_.resize(ni_ * nj_ * nk_, val);
}


template<class A> 
typename Grid<A>::reference Grid<A>::operator()(int i, int j, int k)
{
  return(data_[GetIndex(i, j, k)]);
}


template<class A> 
typename Grid<A>::reference Grid<A>::operator()(int index)
{
#ifdef DEBUG
  if(index < 0 || index >= GetN()) {
    std::ostringstream ost;
    ost << "Index error reading grid: ";
    ost << "Index " << index << " is outside valid range ["
        << 0 << ", " << GetN()-1 << "]";
    throw NRLib2::IndexOutOfRange(ost.str());
  }
#endif
  return(data_[index]);
}


template<class A> 
typename Grid<A>::const_reference Grid<A>::operator()(int i, int j, int k) const
{
  return(data_[GetIndex(i, j, k)]);
}


template<class A> 
typename Grid<A>::const_reference Grid<A>::operator()(int index) const
{
#ifdef DEBUG
  if(index < 0 || index >= GetN()) {
    std::ostringstream ost;
    ost << "Index error reading grid: ";
    ost << "Index " << index << " is outside valid range ["
        << 0 << ", " << GetN()-1 << "]";
    throw NRLib2::IndexOutOfRange(ost.str());
  }
#endif
  return(data_[index]);
}


template<class A> 
int Grid<A>::GetIndex(int i, int j, int k) const
{
#ifdef DEBUG
  if(i < 0 || i >= ni_) {
    std::ostringstream ost;
    ost << "Index error reading grid: ";
    ost << "i index " << i << " is outside valid range ["
        << 0 << ", " << ni_-1 << "]";
    throw NRLib2::IndexOutOfRange(ost.str());
  }
  if(j < 0 || j >= nj_) {
    std::ostringstream ost;
    ost << "Index error reading grid: ";
    ost << "j index " << j << " is outside valid range ["
        << 0 << ", " << nj_-1 << "]";
    throw NRLib2::IndexOutOfRange(ost.str());
  }
  if(k < 0 || k >= nk_) {
    std::ostringstream ost;
    ost << "Index error reading grid: ";
    ost << "k index " << k << " is outside valid range ["
        << 0 << ", " << nk_-1 << "]";
    throw NRLib2::IndexOutOfRange(ost.str());
  }
#endif
  return(i+j*ni_+k*ni_*nj_);
}

template<class A> 
void Grid<A>::GetIJK(int index, int &i, int &j, int &k) const
{
#ifdef DEBUG
  if(index > GetN())
    throw NRLib2::IndexOutOfRange();
#endif
  i = (index % ni_);
  j = ((index-i)/ni_ % nj_);
  k = (index-j*ni_-i)/ni_/nj_;
}


}
#endif




  


