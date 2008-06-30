// $Id: grid2d.hpp 33 2008-04-07 12:00:25Z perroe $

#ifndef NRLIB_GRID2D_HPP
#define NRLIB_GRID2D_HPP

#include <vector>
#include <sstream>

#include "../exception/exception.hpp"

namespace NRLib2 {

template<class A>
class Grid2D {
public:
  typedef typename std::vector<A>::iterator iterator;
  typedef typename std::vector<A>::const_iterator const_iterator;
  typedef typename std::vector<A>::reference reference;
  typedef typename std::vector<A>::const_reference const_reference;

  Grid2D();
    /// \param val Initial cell value. 
  Grid2D(int ni, int nj, const A& val = A());
  virtual ~Grid2D();

    /// All values in the grid are erased when the grid is
    /// resized.
    /// \param val Initial cell value. 
  virtual void Resize(int ni, int nj, const A& val = A());

    /// \throw IndexOutOfRange if one of the indexes was invalid.
  inline reference operator()(int i, int j);
    /// \throw IndexOutOfRange if the index was invalid.
  inline reference operator()(int index);

    /// \throw IndexOutOfRange if one of the indexes was invalid.
  inline const_reference operator()(int i, int j) const;
    /// \throw IndexOutOfRange if the index was invalid.
  inline const_reference operator()(int index) const;

  iterator begin() {return(data_.begin());}
  iterator end()   {return(data_.end());}

  const_iterator begin() const {return(data_.begin());}
  const_iterator end() const   {return(data_.end());}

  int GetNI() const {return(ni_);}
  int GetNJ() const {return(nj_);}
  int GetN()  const {return(static_cast<int>(data_.size()));}

  inline int GetIndex(int i, int j) const;
  void GetIJ(int index, int &i, int &j) const;

private:
  int ni_;
  int nj_;

  std::vector<A> data_;
};

template<class A> 
Grid2D<A>::Grid2D() 
  : ni_(0),
    nj_(0),
    data_()
{}

template<class A> 
Grid2D<A>::Grid2D(int ni, int nj, const A& val) 
  : ni_(ni),
    nj_(nj),
    data_(ni*nj, val)
{}

template<class A>
Grid2D<A>::~Grid2D() 
{}

template<class A> 
void Grid2D<A>::Resize(int ni, int nj, const A& val)
{
  ni_ = ni;
  nj_ = nj;

  data_.resize(0); //To avoid copying of elements
  data_.resize(ni_ * nj_, val);
}


template<class A> 
typename Grid2D<A>::reference Grid2D<A>::operator()(int i, int j)
{
  return(data_[GetIndex(i, j)]);
}


template<class A> 
typename Grid2D<A>::reference Grid2D<A>::operator()(int index)
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
typename Grid2D<A>::const_reference Grid2D<A>::operator()(int i, int j) const
{
  return(data_[GetIndex(i, j)]);
}


template<class A> 
typename Grid2D<A>::const_reference Grid2D<A>::operator()(int index) const
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
int Grid2D<A>::GetIndex(int i, int j) const
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
#endif
  return(i+j*ni_);
}

template<class A> 
void Grid2D<A>::GetIJ(int index, int &i, int &j) const
{
#ifdef DEBUG
  if(index > GetN())
    throw NRLib2::IndexOutOfRange();
#endif
  i = (index % ni_);
  j = ((index-i)/ni_ % nj_);
}

} // namespace NRLib2

#endif // NRLIB_GRID2D_HPP 
