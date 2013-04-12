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

#ifndef FLENS_MATVEC_IO
#define FLENS_MATVEC_IO 1

#include <fstream>
#include <iostream>

#include <flens/crs.h>
#include <flens/densevector.h>
#include <flens/generalmatrix.h>
#include <flens/hermitianmatrix.h>
#include <flens/matvec.h>
#include <flens/sparsematrix.h>
#include <flens/symmetricmatrix.h>
#include <flens/tinymatrix.h>
#include <flens/tinyvector.h>
#include <flens/triangularmatrix.h>

namespace flens {

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const Matrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const Vector<I> &x);

template <typename A>
    std::ostream &
    operator<<(std::ostream &out, const TinyVector<A> &x);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const TinyGeMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const GeMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const GbMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const SyMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const SbMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const SpMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const HeMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const HbMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const HpMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const TrMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const TbMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const TpMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const SparseGeMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const SparseSymmetricMatrix<I> &A);

template <typename I>
    std::ostream &
    operator<<(std::ostream &out, const DenseVector<I> &x);

template <typename T>
    void
    binary(std::ostream &out, const T &t);

template <typename I>
    void
    binary(std::ostream &out, const GeMatrix<I> &A);

template <typename I>
    void
    binary(std::ostream &out, const DenseVector<I> &x);

template <typename T>
    void
    binary(std::ostream &out, const DenseVector<Array<T> > &x);

//-- file stream ouput ---------------------------------------------------------

template <typename I>
    std::ofstream &
    operator<<(std::ofstream &out, const DenseVector<I> &x);

template <typename I>
    std::ofstream &
    operator<<(std::ofstream &out, const GeMatrix<I> &A);

template <typename T>
    std::ofstream &
    operator<<(std::ofstream &out,
               const SparseGeMatrix<CRS<T,CRS_General> > &A);

//-- file stream input ---------------------------------------------------------

template <typename I>
    std::ifstream &
    operator>>(std::ifstream &in, DenseVector<I> &x);

template <typename T>
    std::ifstream &
    operator>>(std::ifstream &in, GeMatrix<FullStorage<T,ColMajor> > &A);

template <typename T>
    std::ifstream &
    operator>>(std::ifstream &in, SparseGeMatrix<CRS<T,CRS_General> > &A);

//-- binary input --------------------------------------------------------------
template <typename T>
    void
    binary(std::istream &in, T &t);

template <typename I>
    void
    binary(std::istream &in, GeMatrix<I> &A);

template <typename I>
    void
    binary(std::istream &in, DenseVector<I> x);

template <typename T>
    void
    binary(std::istream &in, DenseVector<Array<T> > &x);

} // namespace flens

#include <flens/matvecio.tcc>

#endif // FLENS_MATVEC_IO
