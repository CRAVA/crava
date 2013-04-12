/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <exception>

#include <flens/aux_malloc.h>
#include <flens/blas.h>
#include <flens/hacksforgmpxx.h>

namespace flens {

//== ConstFullStorageView ======================================================

template <typename T, StorageOrder Order>
ConstFullStorageView<T, Order>::ConstFullStorageView(const T *data,
                                                     int numRows, int numCols,
                                                     int leadingDimension,
                                                     int firstRow, int firstCol)
    : _data(data),
      _numRows(numRows), _numCols(numCols),
      _leadingDimension(leadingDimension),
      _firstRow(firstRow), _firstCol(firstCol)
{
    //_allocate(data);
}

template <typename T, StorageOrder Order>
ConstFullStorageView<T, Order>::ConstFullStorageView(const ConstView &rhs)
    : _data(rhs._data),
      _numRows(rhs._numRows), _numCols(rhs._numCols),
      _leadingDimension(rhs._leadingDimension),
      _firstRow(rhs._firstRow), _firstCol(rhs._firstCol)
{
    //_allocate(&rhs(_firstRow, _firstCol));
}

template <typename T, StorageOrder Order>
ConstFullStorageView<T, Order>::~ConstFullStorageView()
{
    //free(leadingDimensionStorage());
}

template<StorageOrder>
static inline int GetIndex(int row, int col, int leadingDimension);

template<>
inline int GetIndex<ColMajor>(int row, int col, int leadingDimension)
{
  return row + col*leadingDimension;
}

template<>
inline int GetIndex<RowMajor>(int row, int col, int leadingDimension)
{
  return col + row*leadingDimension;
}

template <typename T, StorageOrder Order>
const T &
ConstFullStorageView<T, Order>::operator()(int row, int col) const
{
    assert(row>=_firstRow);
    assert(row<_firstRow+_numRows);
    assert(col>=_firstCol);
    assert(col<_firstCol+_numCols);

    int index = GetIndex<Order>(row-_firstRow, col-_firstCol, _leadingDimension);
    return _data[index];
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::firstRow() const
{
    return _firstRow;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::firstCol() const
{
    return _firstCol;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::numRows() const
{
    return _numRows;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::numCols() const
{
    return _numCols;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::leadingDimension() const
{
    return _leadingDimension;
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::strideRow() const
{
    return (Order==ColMajor) ? 1
                             : leadingDimension();
}

template <typename T, StorageOrder Order>
int
ConstFullStorageView<T, Order>::strideCol() const
{
    return (Order==ColMajor) ? leadingDimension()
                             : 1;
}

template <typename T, StorageOrder Order>
const T *
ConstFullStorageView<T, Order>::data() const
{
    return _data;
    // return &(this->operator()(_firstRow, _firstCol));
}

//template <typename T, StorageOrder Order>
//const T * const *
//ConstFullStorageView<T, Order>::leadingDimensionStorage() const
//{
//    return (Order==ColMajor) ? &_data[_firstCol]
//                             : &_data[_firstRow];
//}
//
//template <typename T, StorageOrder Order>
//const T **
//ConstFullStorageView<T, Order>::leadingDimensionStorage()
//{
//    return (Order==ColMajor) ? &_data[_firstCol]
//                             : &_data[_firstRow];
//}

template <typename T, StorageOrder Order>
ConstFullStorageView<T, Order>
ConstFullStorageView<T, Order>::view(int fromRow, int fromCol,
                                     int toRow, int toCol,
                                     int firstViewRow, int firstViewCol) const
{
    assert(fromRow>=firstRow());
    assert(fromRow<=toRow);
    assert(toRow<=lastRow());

    assert(fromCol>=firstCol());
    assert(fromCol<=toCol);
    assert(toCol<=lastCol());

    return ConstView(&(this->operator()(fromRow, fromCol)),// data
                     toRow-fromRow+1,                      // # rows
                     toCol-fromCol+1,                      // # cols
                     leadingDimension(),                   // leading dimension
                     firstViewRow,                         // firstRow
                     firstViewCol);                        // firstCol
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
ConstFullStorageView<T, Order>::viewRow(int row, int firstViewIndex) const
{
    assert(row>=firstRow());
    assert(row<=lastRow());

    return ConstArrayView<T>(&(this->operator()(row, _firstCol))-firstViewIndex,
                             numCols(),
                             strideCol(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
ConstFullStorageView<T, Order>::viewCol(int col, int firstViewIndex) const
{
    assert(col>=firstCol());
    assert(col<=lastCol());

    return ConstArrayView<T>(&(this->operator()(_firstRow, col))-firstViewIndex,
                             numRows(),
                             strideRow(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
ConstFullStorageView<T, Order>::viewRow(int row, int fromCol, int toCol,
                                        int firstViewIndex) const
{
    assert(row>=firstRow());
    assert(row<=lastRow());
    assert(fromCol>=firstCol());
    assert(toCol>=fromCol);
    assert(toCol<=lastCol());


    return ConstArrayView<T>(&(this->operator()(row, fromCol))-firstViewIndex,
                             toCol-fromCol+1,
                             strideCol(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
ConstFullStorageView<T, Order>::viewCol(int col, int fromRow, int toRow,
                                        int firstViewIndex) const
{
    assert(col>=firstCol());
    assert(col<=lastCol());
    assert(fromRow>=firstRow());
    assert(toRow>=fromRow);
    assert(toRow<=lastRow());

    return ConstArrayView<T>(&(this->operator()(fromRow, col))-firstViewIndex,
                             toRow-fromRow+1,
                             strideRow(),
                             firstViewIndex);
}

//template <typename T, StorageOrder Order>
//void
//ConstFullStorageView<T, Order>::_allocate(const T *data)
//{
//    assert(!_data);
//    assert(_numRows>0);
//    assert(_numCols>0);
//
//    if (Order==ColMajor) {
//        _data = static_cast<const T **>(calloc(_numCols, sizeof(T *)))
//                - _firstCol;
//        assert(_data+_firstCol);
//
//        _data[_firstCol] = data - _firstRow;
//
//        for (int i=1; i<_numCols; ++i) {
//            _data[_firstCol+i] = _data[_firstCol] + i*_leadingDimension;
//        }
//    }
//    if (Order==RowMajor) {
//        _data = static_cast<const T **>(calloc(_numRows, sizeof(T *)))
//                - _firstRow;
//        assert(_data+_firstRow);
//
//        _data[_firstRow] = data - _firstCol;
//
//        for (int i=1; i<_numRows; ++i) {
//            _data[_firstRow+i] = _data[_firstRow] + i*_leadingDimension;
//        }
//    }
//}

//== FullStorageView ===========================================================

template <typename T, StorageOrder Order>
FullStorageView<T, Order>::FullStorageView(T *data,
                                           int numRows, int numCols,
                                           int leadingDimension,
                                           int firstRow, int firstCol)
    : _data(data),
      _numRows(numRows), _numCols(numCols),
      _leadingDimension(leadingDimension),
      _firstRow(firstRow), _firstCol(firstCol)
{
    //_allocate(data);
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order>::FullStorageView(const FullStorageView<T, Order> &rhs)
    : _data(rhs._data),
      _numRows(rhs._numRows), _numCols(rhs._numCols),
      _leadingDimension(rhs._leadingDimension),
      _firstRow(rhs._firstRow), _firstCol(rhs._firstCol)
{
    //T* data;
    //if (Order==ColMajor) {
    //    data = rhs._data[_firstCol];
    //}
    //else {
    //    data = rhs._data[_firstRow];
    //}

    //_allocate(data);
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order>::~FullStorageView()
{
    //free(leadingDimensionStorage());
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order> &
FullStorageView<T, Order>::operator=(const FullStorage<T, Order> &rhs)
{
    assert(rhs.numRows()==_numRows);
    assert(rhs.numCols()==_numCols);
    copy(Order, _numRows, _numCols,
         rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    return *this;
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order> &
FullStorageView<T, Order>::operator=(const FullStorageView<T, Order> &rhs)
{
    if (this!=&rhs) {
        resize(rhs.numRows(), rhs.numCols(), rhs.firstRow(), rhs.firstCol());
        copy(Order, _numRows, _numCols,
              rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    }
    return *this;
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order> &
FullStorageView<T, Order>::operator=(const ConstFullStorageView<T, Order> &rhs)
{
    resize(rhs.numRows(), rhs.numCols(), rhs.firstRow(), rhs.firstCol());
    copy(Order, _numRows, _numCols,
         rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    return *this;
}

template <typename T, StorageOrder Order>
const T &
FullStorageView<T, Order>::operator()(int row, int col) const
{
    assert(row>=_firstRow);
    assert(row<_firstRow+_numRows);
    assert(col>=_firstCol);
    assert(col<_firstCol+_numCols);

    int index = GetIndex<Order>(row-_firstRow, col-_firstCol, _leadingDimension);
    return _data[index];
}

template <typename T, StorageOrder Order>
T &
FullStorageView<T, Order>::operator()(int row, int col)
{
    assert(row>=_firstRow);
    assert(row<_firstRow+_numRows);
    assert(col>=_firstCol);
    assert(col<_firstCol+_numCols);

    int index = GetIndex<Order>(row-_firstRow, col-_firstCol, _leadingDimension);
    return _data[index];
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order>::operator ConstView() const
{
    return ConstView(_data,                                    // data
                  // &(this->operator()(_firstRow, _firstCol)),// data
                     _numRows,                                 // # rows
                     _numCols,                                 // # cols
                     leadingDimension(),                       // leading dim.
                     _firstRow,                                // first row
                     _firstCol);                               // first col.
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::firstRow() const
{
    return _firstRow;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::firstCol() const
{
    return _firstCol;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::numRows() const
{
    return _numRows;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::numCols() const
{
    return _numCols;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::leadingDimension() const
{
    return _leadingDimension;
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::strideRow() const
{
    return (Order==ColMajor) ? 1
                             : leadingDimension();
}

template <typename T, StorageOrder Order>
int
FullStorageView<T, Order>::strideCol() const
{
    return (Order==ColMajor) ? leadingDimension()
                             : 1;
}

template <typename T, StorageOrder Order>
const T *
FullStorageView<T, Order>::data() const
{
    return _data;
    //return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, StorageOrder Order>
T *
FullStorageView<T, Order>::data()
{
    return _data;
    //return &(this->operator()(_firstRow, _firstCol));
}

//template <typename T, StorageOrder Order>
//const T * const *
//FullStorageView<T, Order>::leadingDimensionStorage() const
//{
//    return (Order==ColMajor) ? &_data[_firstCol]
//                             : &_data[_firstRow];
//}
//
//template <typename T, StorageOrder Order>
//T**
//FullStorageView<T, Order>::leadingDimensionStorage()
//{
//    return (Order==ColMajor) ? &_data[_firstCol]
//                             : &_data[_firstRow];
//}

template <typename T, StorageOrder Order>
ConstFullStorageView<T, Order>
FullStorageView<T, Order>::view(int fromRow, int fromCol,
                                int toRow, int toCol,
                                int firstViewRow, int firstViewCol) const
{
    assert(fromRow>=firstRow());
    assert(fromRow<=toRow);
    assert(toRow<=lastRow());

    assert(fromCol>=firstCol());
    assert(fromCol<=toCol);
    assert(toCol<=lastCol());

    return ConstView(&(this->operator()(fromRow, fromCol)),// data
                     toRow-fromRow+1,                      // # rows
                     toCol-fromCol+1,                      // # cols
                     leadingDimension(),                   // leading dimension
                     firstViewRow,                         // firstRow
                     firstViewCol);                        // firstCol
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order>
FullStorageView<T, Order>::view(int fromRow, int fromCol,
                                int toRow, int toCol,
                                int firstViewRow, int firstViewCol)
{
    assert(fromRow>=firstRow());
    assert(fromRow<=toRow);
    assert(toRow<=lastRow());

    assert(fromCol>=firstCol());
    assert(fromCol<=toCol);
    assert(toCol<=lastCol());

    return View(&(this->operator()(fromRow, fromCol)),// data
                toRow-fromRow+1,                      // # rows
                toCol-fromCol+1,                      // # cols
                leadingDimension(),                   // leading dimension
                firstViewRow,                         // firstRow
                firstViewCol);                        // firstCol
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorageView<T, Order>::viewRow(int row, int firstViewIndex) const
{
    assert(row>=firstRow());
    assert(row<=lastRow());

    return ConstArrayView<T>(&(this->operator()(row, _firstCol))-firstViewIndex,
                             numCols(),
                             strideCol(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorageView<T, Order>::viewRow(int row, int firstViewIndex)
{
    assert(row>=firstRow());
    assert(row<=lastRow());

    return ArrayView<T>(&(this->operator()(row, _firstCol))-firstViewIndex,
                        numCols(),
                        strideCol(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorageView<T, Order>::viewCol(int col, int firstViewIndex) const
{
    assert(col>=firstCol());
    assert(col<=lastCol());

    return ConstArrayView<T>(&(this->operator()(_firstRow, col))-firstViewIndex,
                             numRows(),
                             strideRow(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorageView<T, Order>::viewCol(int col, int firstViewIndex)
{
    assert(col>=firstCol());
    assert(col<=lastCol());

    return ArrayView<T>(&(this->operator()(_firstRow, col))-firstViewIndex,
                        numRows(),
                        strideRow(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorageView<T, Order>::viewRow(int row, int fromCol, int toCol,
                                   int firstViewIndex) const
{
    assert(row>=firstRow());
    assert(row<=lastRow());
    assert(fromCol>=firstCol());
    assert(toCol>=fromCol);
    assert(toCol<=lastCol());

    return ConstArrayView<T>(&(this->operator()(row, fromCol))-firstViewIndex,
                             toCol-fromCol+1,
                             strideCol(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorageView<T, Order>::viewRow(int row, int fromCol, int toCol,
                                   int firstViewIndex)
{
    assert(row>=firstRow());
    assert(row<=lastRow());
    assert(fromCol>=firstCol());
    assert(toCol>=fromCol);
    assert(toCol<=lastCol());

    return ArrayView<T>(&(this->operator()(row, fromCol))-firstViewIndex,
                        toCol-fromCol+1,
                        strideCol(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorageView<T, Order>::viewCol(int col, int fromRow, int toRow,
                                   int firstViewIndex) const
{
    assert(col>=firstCol());
    assert(col<=lastCol());
    assert(fromRow>=firstRow());
    assert(toRow>=fromRow);
    assert(toRow<=lastRow());

    return ConstArrayView<T>(&(this->operator()(fromRow, col))-firstViewIndex,
                             toRow-fromRow+1,
                             strideRow(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorageView<T, Order>::viewCol(int col, int fromRow, int toRow,
                                   int firstViewIndex)
{
    assert(col>=firstCol());
    assert(col<=lastCol());
    assert(fromRow>=firstRow());
    assert(toRow>=fromRow);
    assert(toRow<=lastRow());

    return ArrayView<T>(&(this->operator()(fromRow, col))-firstViewIndex,
                        toRow-fromRow+1,
                        strideRow(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorageView<T, Order>::viewDiag(int d) const
{
    int col = firstCol() + ( (d>0) ? d : 0 );
    int row = firstRow() + ( (d>0) ? 0 : -d );
    return ConstArrayView<T>(&(this->operator()(row,col)) - 1,
                             std::min(numRows(),numCols()) - std::abs(d),
                             leadingDimension()+1);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorageView<T, Order>::viewDiag(int d)
{
    int col = firstCol() + ( (d>0) ? d : 0 );
    int row = firstRow() + ( (d>0) ? 0 : -d );
    return ArrayView<T>(data(),
                        &(this->operator()(row,col)) - 1,
                        std::min(numRows(),numCols()) - std::abs(d),
                        leadingDimension()+1);
}


template <typename T, StorageOrder Order>
void
FullStorageView<T, Order>::resize(int /*numRows*/, int /*numCols*/,
                                  int /*firstRow*/, int /*firstCol*/)
{
//  PERHAPS: resize should not change the size!
//    assert(0); // you can not resize a view
}

//template <typename T, StorageOrder Order>
//void
//FullStorageView<T, Order>::_allocate(T *data)
//{
//    assert(!_data);
//    assert(_numRows>0);
//    assert(_numCols>0);
//
//    if (Order==ColMajor) {
//        _data = static_cast<T **>(calloc(_numCols, sizeof(T *))) - _firstCol;
//        assert(_data+_firstCol);
//
//        _data[_firstCol] = data - _firstRow;
//
//        for (int i=1; i<_numCols; ++i) {
//            _data[_firstCol+i] = _data[_firstCol] + i*_leadingDimension;
//        }
//    }
//    if (Order==RowMajor) {
//        _data = static_cast<T **>(calloc(_numRows, sizeof(T *))) - _firstRow;
//        assert(_data+_firstRow);
//
//        _data[_firstRow] = data - _firstCol;
//
//        for (int i=1; i<_numRows; ++i) {
//            _data[_firstRow+i] = _data[_firstRow] + i*_leadingDimension;
//        }
//    }
//}

//== FullStorage ===============================================================

template <typename T, StorageOrder Order>
FullStorage<T, Order>::FullStorage()
    : _numRows(0), _numCols(0), _firstRow(FLENS_FIRST_INDEX), _firstCol(FLENS_FIRST_INDEX), _data(0)
{
}

template <typename T, StorageOrder Order>
FullStorage<T, Order>::FullStorage(int numRows, int numCols,
                                    int firstRow, int firstCol)
    : _numRows(numRows), _numCols(numCols),
      _firstRow(firstRow), _firstCol(firstCol), _data(0)
{
    _allocate();
}

template <typename T, StorageOrder Order>
FullStorage<T, Order>::FullStorage(const FullStorage<T, Order> &rhs)
    : _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _firstRow(rhs.firstRow()), _firstCol(rhs.firstCol()), _data(0)
{
    if (_numRows != 0 && _numCols != 0) {
        _allocate();
        copy(Order, _numRows, _numCols,
             rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    }
}

template <typename T, StorageOrder Order>
FullStorage<T, Order>::FullStorage(const FullStorageView<T, Order> &rhs)
    : _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _firstRow(rhs.firstRow()), _firstCol(rhs.firstCol()), _data(0)
{
    if (_numRows != 0 && _numCols != 0) {
        _allocate();
        copy(Order, _numRows, _numCols,
             rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    }
}

template <typename T, StorageOrder Order>
FullStorage<T, Order>::FullStorage(const ConstFullStorageView<T, Order> &rhs)
    : _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _firstRow(rhs.firstRow()), _firstCol(rhs.firstCol()), _data(0)
{
    if (_numRows != 0 && _numCols != 0) {
        _allocate();
        copy(Order, _numRows, _numCols,
             rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    }
}


template <typename T, StorageOrder Order>
FullStorage<T, Order>::~FullStorage()
{
    _release();
}

template <typename T, StorageOrder Order>
FullStorage<T, Order> &
FullStorage<T, Order>::operator=(const FullStorage<T, Order> &rhs)
{
    if (this!=&rhs) {
        resize(rhs.numRows(), rhs.numCols(), rhs.firstRow(), rhs.firstCol());
        copy(Order, _numRows, _numCols,
             rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    }
    return *this;
}

template <typename T, StorageOrder Order>
FullStorage<T, Order> &
FullStorage<T, Order>::operator=(const FullStorageView<T, Order> &rhs)
{
    resize(rhs.numRows(), rhs.numCols(), rhs.firstRow(), rhs.firstCol());
    copy(Order, _numRows, _numCols,
         rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    return *this;
}

template <typename T, StorageOrder Order>
FullStorage<T, Order> &
FullStorage<T, Order>::operator=(const ConstFullStorageView<T, Order> &rhs)
{
    resize(rhs.numRows(), rhs.numCols(), rhs.firstRow(), rhs.firstCol());
    copy(Order, _numRows, _numCols,
         rhs.data(), rhs.leadingDimension(), data(), leadingDimension());
    return *this;
}

template<StorageOrder Order>
static inline void GetDataElement(int row, int col, int& first, int& second);

template<>
inline void GetDataElement<ColMajor>(int row, int col,
                                     int& first, int& second)
{
  first = col;
  second = row;
}

template<>
inline void GetDataElement<RowMajor>(int row, int col,
                                     int& first, int& second)
{
  first = row;
  second = col;
}

template <typename T, StorageOrder Order>
inline
const T &
FullStorage<T, Order>::operator()(int row, int col) const
{
    assert(row>=_firstRow);
    assert(row<_firstRow+_numRows);
    assert(col>=_firstCol);
    assert(col<_firstCol+_numCols);

    int first, second;
    GetDataElement<Order>(row, col, first, second);
    return _data[first][second];
}

template <typename T, StorageOrder Order>
inline
T &
FullStorage<T, Order>::operator()(int row, int col)
{
    assert(row>=_firstRow);
    assert(row<_firstRow+_numRows);
    assert(col>=_firstCol);
    assert(col<_firstCol+_numCols);

    int first, second;
    GetDataElement<Order>(row, col, first, second);
    return _data[first][second];
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::firstRow() const
{
    return _firstRow;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::firstCol() const
{
    return _firstCol;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::numRows() const
{
    return _numRows;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::numCols() const
{
    return _numCols;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::leadingDimension() const
{
    return (Order==ColMajor) ? _numRows
                             : _numCols;
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::strideRow() const
{
    return (Order==ColMajor) ? 1
                             : leadingDimension();
}

template <typename T, StorageOrder Order>
int
FullStorage<T, Order>::strideCol() const
{
    return (Order==ColMajor) ? leadingDimension()
                             : 1;
}

template <typename T, StorageOrder Order>
const T *
FullStorage<T, Order>::data() const
{
    return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, StorageOrder Order>
T *
FullStorage<T, Order>::data()
{
    return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, StorageOrder Order>
const T * const *
FullStorage<T, Order>::leadingDimensionStorage() const
{
    return (Order==ColMajor) ? &_data[_firstCol]
                             : &_data[_firstRow];
}

template <typename T, StorageOrder Order>
T**
FullStorage<T, Order>::leadingDimensionStorage()
{
    return (Order==ColMajor) ? &_data[_firstCol]
                             : &_data[_firstRow];
}

template <typename T, StorageOrder Order>
void
FullStorage<T, Order>::resize(int numRows, int numCols,
                               int firstRow, int firstCol)
{
    if ((_numRows!=numRows)
     || (_numCols!=numCols)
     || (_firstRow!=firstRow)
     || (_firstCol!=firstCol)) {
        _release();
        _numRows = numRows;
        _numCols = numCols;
        _firstRow = firstRow;
        _firstCol = firstCol;
        _allocate();
    }
}

template <typename T, StorageOrder Order>
void
FullStorage<T, Order>::resizeOrClear(int numRows, int numCols,
                                     int firstRow, int firstCol)
{
    if ((_numRows!=numRows)
     || (_numCols!=numCols)) {
        _release();
        _numRows = numRows;
        _numCols = numCols;
        _firstRow = firstRow;
        _firstCol = firstCol;
        _allocate();
    } else {
        _firstRow = firstRow;
        _firstCol = firstCol;
        std::fill_n(this->data(), _numRows*_numCols, T());
    }
}

template <typename T, StorageOrder Order>
void
FullStorage<T, Order>::initialize(T val)
{
  std::fill_n(this->data(), _numRows*_numCols, val);
}

template <typename T, StorageOrder Order>
void
FullStorage<T, Order>::shiftIndexTo(int firstRow, int firstCol)
{
    if (Order==RowMajor) {
        _data[firstRow] -= firstCol - _firstCol;
        for (int i=1; i<_numRows; ++i) {
            _data[firstRow+i] = _data[firstRow] + i*_numCols;
        }
        _data -= firstRow - _firstRow;
    }
    if (Order==ColMajor) {
        _data[_firstCol] -= firstRow - _firstRow;
        for (int i=1; i<_numCols; ++i) {
            _data[_firstCol+i] = _data[_firstCol] + i*_numRows;
        }
        _data -= firstCol - _firstCol;
    }
    _firstRow = firstRow;
    _firstCol = firstCol;
}

template <typename T, StorageOrder Order>
ConstFullStorageView<T, Order>
FullStorage<T, Order>::view(int fromRow, int fromCol,
                            int toRow, int toCol,
                            int firstViewRow, int firstViewCol) const
{
    assert(fromRow>=firstRow());
    assert(fromRow<=toRow);
    assert(toRow<=lastRow());

    assert(fromCol>=firstCol());
    assert(fromCol<=toCol);
    assert(toCol<=lastCol());

    return ConstView(&(this->operator()(fromRow, fromCol)),// data
                     toRow-fromRow+1,                      // # rows
                     toCol-fromCol+1,                      // # cols
                     leadingDimension(),                   // leading dimension
                     firstViewRow,                         // firstRow
                     firstViewCol);                        // firstCol
}

template <typename T, StorageOrder Order>
FullStorageView<T, Order>
FullStorage<T, Order>::view(int fromRow, int fromCol,
                            int toRow, int toCol,
                            int firstViewRow, int firstViewCol)
{
    assert(fromRow>=firstRow());
    assert(fromRow<=toRow);
    assert(toRow<=lastRow());

    assert(fromCol>=firstCol());
    assert(fromCol<=toCol);
    assert(toCol<=lastCol());

    return View(&(this->operator()(fromRow, fromCol)),// data
                toRow-fromRow+1,                      // # rows
                toCol-fromCol+1,                      // # cols
                leadingDimension(),                   // leading dimension
                firstViewRow,                         // firstRow
                firstViewCol);                        // firstCol
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorage<T, Order>::viewRow(int row, int firstViewIndex) const
{
    assert(row>=firstRow());
    assert(row<=lastRow());

    return ConstArrayView<T>(&(this->operator()(row, _firstCol))-firstViewIndex,
                             numCols(),
                             strideCol(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorage<T, Order>::viewRow(int row, int firstViewIndex)
{
    assert(row>=firstRow());
    assert(row<=lastRow());

    return ArrayView<T>(&(this->operator()(row, _firstCol))-firstViewIndex,
                        numCols(),
                        strideCol(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorage<T, Order>::viewCol(int col, int firstViewIndex) const
{
    assert(col>=firstCol());
    assert(col<=lastCol());

    return ConstArrayView<T>(&(this->operator()(_firstRow, col))-firstViewIndex,
                             numRows(),
                             strideRow(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorage<T, Order>::viewCol(int col, int firstViewIndex)
{
    assert(col>=firstCol());
    assert(col<=lastCol());

    return ArrayView<T>(&(this->operator()(_firstRow, col))-firstViewIndex,
                        numRows(),
                        strideRow(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorage<T, Order>::viewRow(int row, int fromCol, int toCol,
                                int firstViewIndex) const
{
    assert(row>=firstRow());
    assert(row<=lastRow());
    assert(fromCol>=firstCol());
    assert(toCol>=fromCol);
    assert(toCol<=lastCol());

    return ConstArrayView<T>(&(this->operator()(row, fromCol))-firstViewIndex,
                             toCol-fromCol+1,
                             strideCol(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorage<T, Order>::viewRow(int row, int fromCol, int toCol,
                                int firstViewIndex)
{
    assert(row>=firstRow());
    assert(row<=lastRow());
    assert(fromCol>=firstCol());
    assert(toCol>=fromCol);
    assert(toCol<=lastCol());

    return ArrayView<T>(&(this->operator()(row, fromCol))-firstViewIndex,
                        toCol-fromCol+1,
                        strideCol(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorage<T, Order>::viewCol(int col, int fromRow, int toRow,
                                int firstViewIndex) const
{
    assert(col>=firstCol());
    assert(col<=lastCol());
    assert(fromRow>=firstRow());
    assert(toRow>=fromRow);
    assert(toRow<=lastRow());

    return ConstArrayView<T>(&(this->operator()(fromRow, col))-firstViewIndex,
                             toRow-fromRow+1,
                             strideRow(),
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorage<T, Order>::viewCol(int col, int fromRow, int toRow,
                    int firstViewIndex)
{
    assert(col>=firstCol());
    assert(col<=lastCol());
    assert(fromRow>=firstRow());
    assert(toRow>=fromRow);
    assert(toRow<=lastRow());

    return ArrayView<T>(&(this->operator()(fromRow, col))-firstViewIndex,
                        toRow-fromRow+1,
                        strideRow(),
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
ConstArrayView<T>
FullStorage<T, Order>::viewDiag(int d, int firstViewIndex) const
{
    int col = firstCol() + ( (d>0) ? d : 0 );
    int row = firstRow() + ( (d>0) ? 0 : -d );
    return ConstArrayView<T>(&(this->operator()(row,col)) - 1,
                             std::min(numRows(),numCols()) - std::abs(d),
                             leadingDimension()+1,
                             firstViewIndex);
}

template <typename T, StorageOrder Order>
ArrayView<T>
FullStorage<T, Order>::viewDiag(int d, int firstViewIndex)
{
    int col = firstCol() + ( (d>0) ? d : 0 );
    int row = firstRow() + ( (d>0) ? 0 : -d );
    return ArrayView<T>(&(this->operator()(row,col)) - 1,
                        std::min(numRows(),numCols()) - std::abs(d),
                        leadingDimension()+1,
                        firstViewIndex);
}

template <typename T, StorageOrder Order>
void
FullStorage<T, Order>::_allocate()
{
    assert(!_data);
    assert(_numRows>0);
    assert(_numCols>0);

    if (Order==ColMajor) {
        _data = static_cast<T **>(flens_malloc(_numCols * sizeof(T *))) - _firstCol;
        if (!(_data+_firstCol))
            throw std::bad_alloc();

        _data[_firstCol] = static_cast<T *>(flens_malloc(_numCols * _numRows * sizeof(T)))
                         - _firstRow;
        if (!(_data[_firstCol]+_firstRow))
            throw std::bad_alloc();

        for (int i=1; i<_numCols; ++i) {
            _data[_firstCol+i] = _data[_firstCol] + i*_numRows;
        }
    }
    if (Order==RowMajor) {
        _data = static_cast<T **>(flens_malloc(_numRows * sizeof(T *))) - _firstRow;
        if (!(_data+_firstRow))
            throw std::bad_alloc();

        _data[_firstRow] = static_cast<T *>(flens_malloc(_numCols * _numRows * sizeof(T)))
                         - _firstCol;
        if (!(_data[_firstRow]+_firstCol))
            throw std::bad_alloc();

        for (int i=1; i<_numRows; ++i) {
            _data[_firstRow+i] = _data[_firstRow] + i*_numCols;
        }
    }
    Initializer<FullStorage<T,Order> >::initialize(*this);
}

template <typename T, StorageOrder Order>
void
FullStorage<T, Order>::_release()
{
    if (_data) {
        flens_free(data());
        flens_free(leadingDimensionStorage());
        _data=0;
    }
}

//==============================================================================

} // namespace flens
