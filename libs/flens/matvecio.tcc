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

#define ENTER    static bool enter = false; assert(!enter); enter=true;
#define LEAVE    enter=false;

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

#include <flens/range.h>

namespace flens {

/*
using std::complex;

template <typename T>
std::ostream &
operator<<(std::ostream &out, const std::complex<T> &c)
{
    out.precision(6);
    out.setf(std::ios::fixed);

    out.width(14); out << "(" << c.real();
    if (c.imag()>=0) {
        out << "+";
    }
    out.width(12); out << c.imag() << "i) ";
    return out;
}
*/

template <typename I>
std::ostream &
operator<<(std::ostream &out, const Matrix<I> &A)
{
    ENTER;
    out << A.impl();
    LEAVE;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const Vector<I> &x)
{
    ENTER;
    out << x.impl();
    LEAVE;
    return out;
}

template <typename A>
std::ostream &
operator<<(std::ostream &out, const TinyVector<A> &x)
{
    const int N = x.length();
    out << "[0.." << N-1 << "]" << std::endl;
    out.setf(std::ios::fixed);
    for (int i=0; i<N; ++i) {
        out << x(i) << "  ";
    }
    out << std::endl;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const DenseVector<I> &x)
{

/*
    out << std::endl << "[";
    (x.firstIndex()==1) ? out << x.length()
                        : out << x.firstIndex() << ".." << x.lastIndex();
    out << "] ";
*/
    out << "[ ";
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
//        out.width(11);
        out << x(i) << " ";
        if (i<x.lastIndex()) {
            out << " ";
        }
    }
    out << "]" << std::endl;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const TinyGeMatrix<I> &A)
{
    out << "[" << A.numRows() << ", " << A.numCols() << "]" << std::endl;
    out.setf(std::ios::fixed);
    for (int i=0; i<A.numRows(); ++i) {
        for (int j=0; j<A.numCols(); ++j) {
            out.width(11);
            out << A(i,j) << " ";
        }
        out << ";" << std::endl;
    }
    return out;
}


template <typename I>
std::ostream &
operator<<(std::ostream &out, const GeMatrix<I> &A)
{
    out << "[";
    (A.firstRow()==1) ? out << A.numRows()
                      : out << A.firstRow() << ".." << A.lastRow();
    out << ", ";
    (A.firstCol()==1) ? out << A.numCols()
                      : out << A.firstCol() << ".." << A.lastCol();
    out << "]" << std::endl;
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out.width(11);
            out << A(i,j) << " ";
        }
        out << ";" << std::endl;
    }
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const GbMatrix<I> &A)
{
    if ((A.firstRow()==1) && (A.firstCol()==1)) {
        out << "[" <<  A.numRows() << ", " << A.numCols() << "]" << std::endl;
    } else {
        out << "[" << A.firstRow() << ".." << A.lastRow()
           << ", " << A.firstCol() << ".." << A.lastCol()
           << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out.width(11);
            if ((j-i<=A.numSuperDiags()) && (j-i>=-A.numSubDiags())) {
                out << A(i,j);
            } else {
                out << 0;
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const SyMatrix<I> &A)
{
    if ((A.firstRow()==1) && (A.firstCol()==1)) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out.width(11);
            out << ((A.upLo()==Upper) ? A(std::min(i,j), std::max(i,j))
                                      : A(std::max(i,j), std::min(i,j)));
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename E>
std::ostream &
operator<<(std::ostream &out, const SbMatrix<E> &A)
{
    if (A.firstIndex()==1) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (int j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            out.width(11);
            int I = (A.upLo()==Upper) ? std::min(i,j) : std::max(i,j);
            int J = (A.upLo()==Upper) ? std::max(i,j) : std::min(i,j);
            if (std::abs(I-J)>A.numOffDiags()) {
                out << 0;
            } else {
                out << A(I,J);
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename E>
std::ostream &
operator<<(std::ostream &out, const SpMatrix<E> &A)
{
    if (A.firstIndex()==1) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (int j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            out.width(11);
            int I = (A.upLo()==Upper) ? std::min(i,j) : std::max(i,j);
            int J = (A.upLo()==Upper) ? std::max(i,j) : std::min(i,j);
            out << A(I,J);
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const HeMatrix<I> &A)
{
    if ((A.firstRow()==1) && (A.firstCol()==1)) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out.width(22);
            if (i==j) {
                out << real(A(i,j));
            }
            if (i<j) {
                out << ((A.upLo()==Upper) ? A(i,j) : conj(A(j,i)));
            }
            if (i>j) {
                out << ((A.upLo()==Upper) ? conj(A(j,i)) : A(i,j));
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename E>
std::ostream &
operator<<(std::ostream &out, const HbMatrix<E> &A)
{
    if (A.firstIndex()==1) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (int j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            out.width(22);
            int I = (A.upLo()==Upper) ? std::min(i,j) : std::max(i,j);
            int J = (A.upLo()==Upper) ? std::max(i,j) : std::min(i,j);
            if (std::abs(I-J)>A.numOffDiags()) {
                out << 0;
                continue;
            }
            if (i==j) {
                out << real(A(i,j));
            }
            if (i<j) {
                out << ((A.upLo()==Upper) ? A(i,j) : conj(A(j,i)));
            }
            if (i>j) {
                out << ((A.upLo()==Upper) ? conj(A(j,i)) : A(i,j));
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename E>
std::ostream &
operator<<(std::ostream &out, const HpMatrix<E> &A)
{
    if (A.firstIndex()==1) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (int j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            out.width(22);
            if (i==j) {
                out << real(A(i,j));
            }
            if (i<j) {
                out << ((A.upLo()==Upper) ? A(i,j) : conj(A(j,i)));
            }
            if (i>j) {
                out << ((A.upLo()==Upper) ? conj(A(j,i)) : A(i,j));
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const TrMatrix<I> &A)
{
    if ((A.firstRow()==1) && (A.firstCol()==1)) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstRow() << ".." << A.lastRow()
            << ", " << A.firstCol() << ".." << A.lastCol()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out.width(11);
            if (i==j) {
                (A.unitDiag()==Unit) ? out << 1
                                     : out << A(i,j);
                continue;
            }
            if (((i>j) && (A.upLo()==Lower))
             || ((i<j) && (A.upLo()==Upper))) {
                out << A(i,j);
            } else {
                out << 0;
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename E>
std::ostream &
operator<<(std::ostream &out, const TbMatrix<E> &A)
{
    if (A.firstIndex()==1) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (int j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            out.width(11);
            if ((std::abs(i-j)<=A.numOffDiags())
              && (((A.upLo()==Upper) && (i<=j))
               || ((A.upLo()==Lower) && (i>=j)))) {
                out << A(i,j);
            } else {
                out << 0;
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename E>
std::ostream &
operator<<(std::ostream &out, const TpMatrix<E> &A)
{
    if (A.firstIndex()==1) {
        out << "[" << A.dim() << ", " << A.dim()
            << "]" << std::endl;
    } else {
        out << "["  << A.firstIndex() << ".." << A.lastIndex()
            << ", " << A.firstIndex() << ".." << A.lastIndex()
            << "]" << std::endl;
    }
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (int i=A.firstIndex(); i<=A.lastIndex(); ++i) {
        for (int j=A.firstIndex(); j<=A.lastIndex(); ++j) {
            out.width(11);
            if (((A.upLo()==Upper) && (i<=j))
             || ((A.upLo()==Lower) && (i>=j))) {
                out << A(i,j);
            } else {
                out << 0;
            }
        }
        out << ";" << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const SparseGeMatrix<I> &A)
{
    typedef typename SparseGeMatrix<I>::const_iterator It;

    out << "[" << A.numRows() << ", " << A.numCols() << "]" << std::endl;
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (It it=A.begin(); it!=A.end(); ++it) {
        out << " (" << it->first.first << ", " << it->first.second << "): ";
        out.width(11);
        out << it->second;

        out << std::endl;
    }
    return out;
}

template <typename I>
std::ostream &
operator<<(std::ostream &out, const SparseSymmetricMatrix<I> &A)
{
    typedef typename SparseSymmetricMatrix<I>::const_iterator It;

    out << "[" << A.numRows() << ", " << A.numCols() << "]" << std::endl;
//    out.precision(5);
    out.setf(std::ios::fixed);
    for (It it=A.begin(); it!=A.end(); ++it) {
        out << " (" << it->first.first << ", " << it->first.second << "): ";
        out.width(11);
        out << it->second
            << std::endl;
    }
    return out;
}

//-- file stream output --------------------------------------------------------

template <typename I>
std::ofstream &
operator<<(std::ofstream &out, const DenseVector<I> &x)
{
    out.setf(std::ios::fixed);
    out << x.range() << std::endl;
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        out << x(i) << std::endl;
    }
    return out;
}

template <typename I>
std::ofstream &
operator<<(std::ofstream &out, const GeMatrix<I> &A)
{
    out.setf(std::ios::fixed);
    out << A.rows() << " " << A.cols() << std::endl;
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out << A(i,j) << " ";
        }
        out << std::endl;
    }
    out << std::endl;
    return out;
}

template <typename T>
std::ofstream &
operator<<(std::ofstream &out, const SparseGeMatrix<CRS<T,CRS_General> > &A)
{
    out.setf(std::ios::fixed);
    const CRS<T,CRS_General> &crs=A.engine();
    out << crs.numRows() << " " << crs.numCols() << std::endl;
    out << crs.values << std::endl;
    out << crs.columns << std::endl;
    out << crs.rows << std::endl;
    return out;
}

//-- file stream input --------------------------------------------------------

template <typename I>
std::ifstream &
operator>>(std::ifstream &in, DenseVector<I> &x)
{
    Range range;
    in >> range;
    x.resize(range);
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        in >> x(i);
    }
    return in;
}

template <typename T>
std::ifstream &
operator>>(std::ifstream &in, GeMatrix<FullStorage<T,ColMajor> > &A)
{
    Range rows, cols;
    in >> rows >> cols;
    A.resize(rows,cols);
    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            in >> A(i,j);
        }
    }
    return in;
}

template <typename T>
std::ifstream &
operator>>(std::ifstream &in, SparseGeMatrix<CRS<T,CRS_General> > &A)
{
    int numRows, numCols;
    in >> numRows >> numCols;
    A.resize(numRows,numCols);
    CRS<T,CRS_General> &crs=A.engine();
    in >> crs.values;
    in >> crs.columns;
    in >> crs.rows;
    return in;
}

//-- binary output -------------------------------------------------------------

template <typename T>
void
binary(std::ostream &out, const T &t)
{
    out.write(reinterpret_cast<const char *>(&t), sizeof(T));
}

template <typename I>
void
binary(std::ostream &out, const GeMatrix<I> &A)
{
    typedef typename I::ElementType T;
    StorageOrder order = StorageInfo<I>::order;

    if (order==RowMajor) {
        for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
            out.write(reinterpret_cast<const char *>(A(i,_).data()), sizeof(T)*A.numCols());
        }
    } else {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            out.write(reinterpret_cast<const char *>(A(_,j).data()), sizeof(T)*A.numRows());
        }
    }
}

template <typename I>
void
binary(std::ostream &out, const DenseVector<I> &x)
{
    typedef typename I::ElementType T;

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        out.write(reinterpret_cast<const char *>(&x(i)), sizeof(T));
    }
}

template <typename T>
void
binary(std::ostream &out, const DenseVector<Array<T> > &x)
{
    out.write(reinterpret_cast<const char *>(x.data()), sizeof(T)*x.length());
}

//-- binary input --------------------------------------------------------------
template <typename T>
void
binary(std::istream &in, T &t)
{
    in.read(reinterpret_cast<char *>(&t), sizeof(T));
}

template <typename I>
void
binary(std::istream &in, GeMatrix<I> &A)
{
    typedef typename I::ElementType T;
    StorageOrder order = StorageInfo<I>::order;

    if (order==RowMajor) {
        for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
            in.read(reinterpret_cast<char *>(A(i,_).data()), sizeof(T)*A.numCols());
        }
    } else {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            in.read(reinterpret_cast<char *>(A(_,j).data()), sizeof(T)*A.numRows());
        }
    }
}

template <typename I>
void
binary(std::istream &in, DenseVector<I> x)
{
    typedef typename I::ElementType T;

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        in.read(reinterpret_cast<char *>(&x(i)), sizeof(T));
    }
}

template <typename T>
void
binary(std::istream &in, DenseVector<Array<T> > &x)
{
    in.read(reinterpret_cast<char *>(x.data()), sizeof(T)*x.length());
}

} // namespace flens
