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

namespace flens {

// == SyMatrix =================================================================

template <typename FS>
SyMatrix<FS>::SyMatrix()
{
}

template <typename FS>
SyMatrix<FS>::SyMatrix(int dim, StorageUpLo upLo, int firstIndex)
    : _fs(dim, dim, firstIndex, firstIndex), _upLo(upLo)
{
}

template <typename FS>
SyMatrix<FS>::SyMatrix(const FS &fs, StorageUpLo upLo)
    : _fs(fs), _upLo(upLo)
{
}

template <typename FS>
SyMatrix<FS>::SyMatrix(const SyMatrix<FS> &rhs)
    : SymmetricMatrix<SyMatrix<FS> >(),
      _fs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SyMatrix<FS>::SyMatrix(const SyMatrix<RHS> &rhs)
    : _fs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
SyMatrix<FS>::SyMatrix(const TrMatrix<FS> &rhs)
    : _fs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

template <typename FS>
template <typename RHS>
SyMatrix<FS>::SyMatrix(const TrMatrix<RHS> &rhs)
    : _fs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

template <typename FS>
template <typename RHS>
SyMatrix<FS>::SyMatrix(const Matrix<RHS> &rhs)
{
    copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
SyMatrix<FS> &
SyMatrix<FS>::operator=(const SyMatrix<FS> &rhs)
{
    copy(rhs, *this);
    return *this;
}

template <typename FS>
SyMatrix<FS> &
SyMatrix<FS>::operator=(const TrMatrix<FS> &rhs)
{
    copy(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
SyMatrix<FS> &
SyMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    copy(rhs.impl(), *this);
    return *this;
}

template <typename FS>
SyMatrix<FS> &
SyMatrix<FS>::operator+=(const SyMatrix<FS>& rhs)
{
    axpy(T(1), rhs, *this);
    return *this;
}

template <typename FS>
SyMatrix<FS> &
SyMatrix<FS>::operator*=(T alpha)
{
    scal(alpha, *this);
    return *this;
}

template <typename FS>
SyMatrix<FS> &
SyMatrix<FS>::operator/=(T alpha)
{
    scal(T(1)/alpha, *this);
    return *this;
}

template <typename FS>
const typename SyMatrix<FS>::T &
SyMatrix<FS>::operator()(int row, int col) const
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
    } else {
        assert(col<=row);
    }
#endif
    return _fs(row, col);
}

template <typename FS>
typename SyMatrix<FS>::T &
SyMatrix<FS>::operator()(int row, int col)
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
    } else {
        assert(col<=row);
    }
#endif
    return _fs(row, col);
}

// -- views --------------------------------------------------------------------

// general views
template <typename FS>
typename SyMatrix<FS>::ConstGeneralView
SyMatrix<FS>::general() const
{
    return ConstGeneralView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                            _fs.lastCol(), _fs.lastRow(),
                            firstRow(), firstCol()));
}

template <typename FS>
typename SyMatrix<FS>::GeneralView
SyMatrix<FS>::general()
{
    return GeneralView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                       _fs.lastCol(), _fs.lastRow(),
                       firstRow(), firstCol()));
}

// triangular views
template <typename FS>
typename SyMatrix<FS>::ConstTriangularView
SyMatrix<FS>::triangular() const
{
    return ConstTriangularView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                               _fs.lastCol(), _fs.lastRow(),
                               firstRow(), firstCol()), upLo(), NonUnit);
}

template <typename FS>
typename SyMatrix<FS>::TriangularView
SyMatrix<FS>::triangular()
{
    return TriangularView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                          _fs.lastCol(), _fs.lastRow(),
                          firstRow(), firstCol()), upLo(), NonUnit);
}

// rectangular view

template <typename FS>
typename SyMatrix<FS>::ConstGeneralView
SyMatrix<FS>::operator()(const Range &rows, const Range &cols,
                         int firstViewRow, int firstViewCol) const
{
    assert(rows.stride()==1);
    assert(cols.stride()==1);
    if (_upLo == Upper) {
      assert(cols.lastIndex() <= rows.firstIndex());
    }
    else {
      assert(rows.lastIndex() >= cols.firstIndex());
    }

    return _fs.view(rows.firstIndex(), cols.firstIndex(),
                    rows.lastIndex(), cols.lastIndex(),
                    firstViewRow, firstViewCol);
}

template <typename FS>
typename SyMatrix<FS>::GeneralView
SyMatrix<FS>::operator()(const Range &rows, const Range &cols,
                         int firstViewRow, int firstViewCol)
{
    assert(rows.stride()==1);
    assert(cols.stride()==1);
    if (_upLo == Upper) {
      assert(cols.firstIndex() >= rows.lastIndex());
    }
    else {
      assert(rows.firstIndex() >= cols.lastIndex());
    }

    return _fs.view(rows.firstIndex(), cols.firstIndex(),
                    rows.lastIndex(), cols.lastIndex(),
                    firstViewRow, firstViewCol);
}

template <typename FS>
typename SyMatrix<FS>::ConstSymmetricView
SyMatrix<FS>::diagSub(const Range & rows,
                      int firstViewIndex) const
{
    assert(rows.stride() == 1);

    ConstFSView view = _fs.view(rows.firstIndex(), rows.firstIndex(),
                                rows.lastIndex(),  rows.lastIndex(),
                                firstViewIndex, firstViewIndex);

    return ConstSymmetricView(view, _upLo);
}

template <typename FS>
typename SyMatrix<FS>::SymmetricView
SyMatrix<FS>::diagSub(const Range & rows,
                      int firstViewIndex)
{
    assert(rows.stride() == 1);

    FSView view = _fs.view(rows.firstIndex(), rows.firstIndex(),
                           rows.lastIndex(),  rows.lastIndex(),
                           firstViewIndex, firstViewIndex);

    return SymmetricView(view, _upLo);
}


// -- methods ------------------------------------------------------------------

// for BLAS/LAPACK
template <typename FS>
StorageUpLo
SyMatrix<FS>::upLo() const
{
    return _upLo;
}

template <typename FS>
int
SyMatrix<FS>::dim() const
{
    assert(_fs.numRows()==_fs.numCols());

    return _fs.numRows();
}

template <typename FS>
int
SyMatrix<FS>::leadingDimension() const
{
    return _fs.leadingDimension();
}

template <typename FS>
const typename SyMatrix<FS>::T *
SyMatrix<FS>::data() const
{
    return _fs.data();
}

template <typename FS>
typename SyMatrix<FS>::T *
SyMatrix<FS>::data()
{
    return _fs.data();
}

// for element access
template <typename FS>
int
SyMatrix<FS>::firstRow() const
{
    return _fs.firstRow();
}

template <typename FS>
int
SyMatrix<FS>::lastRow() const
{
    return _fs.lastRow();
}

template <typename FS>
int
SyMatrix<FS>::firstCol() const
{
    return _fs.firstCol();
}

template <typename FS>
int
SyMatrix<FS>::lastCol() const
{
    return _fs.lastCol();
}

template <typename FS>
Range
SyMatrix<FS>::rows()  const
{
    return _(firstRow(), lastRow());
}

template <typename FS>
Range
SyMatrix<FS>::cols() const
{
    return _(firstCol(), lastCol());
}

template <typename FS>
void
SyMatrix<FS>::resize(int dim, StorageUpLo upLo,
                     int firstIndex)
{
  _upLo = upLo;
  _fs.resize(dim, dim, firstIndex, firstIndex);
}

template <typename FS>
void
SyMatrix<FS>::initialize(T val)
{
  _fs.initialize(val);
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const FS &
SyMatrix<FS>::engine() const
{
    return _fs;
}

template <typename FS>
FS &
SyMatrix<FS>::engine()
{
    return _fs;
}

// == SbMatrix =================================================================

template <typename BS>
SbMatrix<BS>::SbMatrix()
{
}

template <typename BS>
SbMatrix<BS>::SbMatrix(int dim, StorageUpLo upLo, int numOffDiags,
                       int firstIndex)
    : _bs(dim, dim,
          (upLo==Lower) ? numOffDiags : 0,
          (upLo==Upper) ? numOffDiags : 0,
          firstIndex),
      _upLo(upLo)
{
}

template <typename BS>
SbMatrix<BS>::SbMatrix(const BS &bs, StorageUpLo upLo)
    : _bs(bs), _upLo(upLo)
{
}

template <typename BS>
SbMatrix<BS>::SbMatrix(const SbMatrix<BS> &rhs)
    : SymmetricMatrix<SbMatrix<BS> >(),
      _bs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename BS>
template <typename RHS>
SbMatrix<BS>::SbMatrix(const SbMatrix<RHS> &rhs)
    : _bs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename BS>
SbMatrix<BS>::SbMatrix(const TbMatrix<BS> &rhs)
    : _bs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

template <typename BS>
template <typename RHS>
SbMatrix<BS>::SbMatrix(const TbMatrix<RHS> &rhs)
    : _bs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

// -- operators ----------------------------------------------------------------

template <typename BS>
SbMatrix<BS> &
SbMatrix<BS>::operator*=(T alpha)
{
    scal(alpha, *this);
    return *this;
}

template <typename BS>
SbMatrix<BS> &
SbMatrix<BS>::operator/=(T alpha)
{
    scal(T(1)/alpha, *this);
    return *this;
}

template <typename BS>
const typename SbMatrix<BS>::T &
SbMatrix<BS>::operator()(int row, int col) const
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
        assert(col-row<=numOffDiags());
    } else {
        assert(col<=row);
        assert(row-col<=numOffDiags());
    }
#endif
    return _bs(row, col);
}

template <typename BS>
typename SbMatrix<BS>::T &
SbMatrix<BS>::operator()(int row, int col)
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
        assert(col-row<=numOffDiags());
    } else {
        assert(col<=row);
        assert(row-col<=numOffDiags());
    }
#endif
    return _bs(row, col);
}

// -- views --------------------------------------------------------------------

template <typename BS>
typename SbMatrix<BS>::ConstVectorView
SbMatrix<BS>::diag(int d) const
{
    return _bs.viewDiag(d);
}

template <typename BS>
typename SbMatrix<BS>::VectorView
SbMatrix<BS>::diag(int d)
{
    return _bs.viewDiag(d);
}

// general views
template <typename BS>
typename SbMatrix<BS>::ConstGeneralView
SbMatrix<BS>::general() const
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return ConstGeneralView(view);
}

template <typename BS>
typename SbMatrix<BS>::GeneralView
SbMatrix<BS>::general()
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return GeneralView(view);
}

// triangular views
template <typename BS>
typename SbMatrix<BS>::ConstTriangularView
SbMatrix<BS>::triangular() const
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return ConstTriangularView(view, upLo());
}

template <typename BS>
typename SbMatrix<BS>::TriangularView
SbMatrix<BS>::triangular()
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return TriangularView(view, upLo());
}

// -- methods ------------------------------------------------------------------

// for BLAS/LAPACK
template <typename BS>
StorageUpLo
SbMatrix<BS>::upLo() const
{
    return _upLo;
}

template <typename BS>
int
SbMatrix<BS>::dim() const
{
    assert (_bs.numRows()==_bs.numCols());

    return _bs.numRows();
}

template <typename BS>
int
SbMatrix<BS>::numOffDiags() const
{
    return (_upLo==Upper) ? _bs.numSuperDiags()
                          : _bs.numSubDiags();
}

template <typename BS>
int
SbMatrix<BS>::leadingDimension() const
{
    return _bs.leadingDimension();
}

template <typename BS>
const typename SbMatrix<BS>::T *
SbMatrix<BS>::data() const
{
    return _bs.data();
}

template <typename BS>
typename SbMatrix<BS>::T *
SbMatrix<BS>::data()
{
    return _bs.data();
}

// for element access
template <typename BS>
int
SbMatrix<BS>::firstIndex() const
{
    return _bs.firstRow();
}

template <typename BS>
int
SbMatrix<BS>::lastIndex() const
{
    assert(_bs.lastRow()==_bs.lastCol());
    return _bs.lastRow();
}

template <typename BS>
Range
SbMatrix<BS>::indices()  const
{
    return _(firstIndex(), lastIndex());
}

template <typename BS>
Range
SbMatrix<BS>::diags() const
{
    return (_upLo==Upper) ? _(0, numOffDiags())
                          : _(-numOffDiags(),0);
}

// -- implementation -----------------------------------------------------------

template <typename BS>
const BS &
SbMatrix<BS>::engine() const
{
    return _bs;
}

template <typename BS>
BS &
SbMatrix<BS>::engine()
{
    return _bs;
}

// == SpMatrix =================================================================

template <typename PS>
SpMatrix<PS>::SpMatrix()
{
}

template <typename PS>
SpMatrix<PS>::SpMatrix(int dim, int firstIndex)
    : _ps(dim, firstIndex)
{
}

template <typename PS>
SpMatrix<PS>::SpMatrix(const PS &ps)
    : _ps(ps)
{
}

// -- operators ----------------------------------------------------------------

template <typename PS>
const typename SpMatrix<PS>::T &
SpMatrix<PS>::operator()(int row, int col) const
{
    return _ps(row, col);
}

template <typename PS>
typename SpMatrix<PS>::T &
SpMatrix<PS>::operator()(int row, int col)
{
    return _ps(row, col);
}

// -- views --------------------------------------------------------------------

// triangular views
template <typename PS>
typename SpMatrix<PS>::ConstTriangularView
SpMatrix<PS>::triangular(UnitDiag unitDiag) const
{
    return ConstTriangularView(_ps.view(), unitDiag);
}

template <typename PS>
typename SpMatrix<PS>::TriangularView
SpMatrix<PS>::triangular(UnitDiag unitDiag)
{
    return TriangularView(_ps.view(), unitDiag);
}

// -- methods ------------------------------------------------------------------

// for BLAS/LAPACK
template <typename PS>
StorageUpLo
SpMatrix<PS>::upLo() const
{
    return StorageInfo<PS>::upLo;
}

template <typename PS>
int
SpMatrix<PS>::dim() const
{
    return _ps.dim();
}

template <typename PS>
const typename SpMatrix<PS>::T *
SpMatrix<PS>::data() const
{
    return _ps.data();
}

template <typename PS>
typename SpMatrix<PS>::T *
SpMatrix<PS>::data()
{
    return _ps.data();
}

// for element access
template <typename PS>
int
SpMatrix<PS>::firstIndex() const
{
    return _ps.firstIndex();
}

template <typename PS>
int
SpMatrix<PS>::lastIndex() const
{
    return _ps.lastIndex();
}

template <typename PS>
Range
SpMatrix<PS>::indices() const
{
    return _(firstIndex(), lastIndex());
}

// -- implementation -----------------------------------------------------------

template <typename PS>
const PS &
SpMatrix<PS>::engine() const
{
    return _ps;
}

template <typename PS>
PS &
SpMatrix<PS>::engine()
{
    return _ps;
}

} // namespace flens
