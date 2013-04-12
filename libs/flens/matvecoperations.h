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

#ifndef FLENS_MATVECOPERATIONS_H
#define FLENS_MATVECOPERATIONS_H 1

#include <flens/matvecclosures.h>

namespace flens {

//-- VectorClosure -------------------------------------------------------------

// -x
template <typename V>
const VectorClosure<OpMult,
                    typename V::Impl,
                    Scalar<typename V::ElementType> >
operator-(const Vector<V> &x);

// alpha * x
template <typename V>
const VectorClosure<OpMult,
                    typename V::Impl,
                    Scalar<typename V::ElementType> >
operator*(const typename V::ElementTypeRef alpha, const Vector<V> &x);

// x * alpha
template <typename V>
const VectorClosure<OpMult,
                    typename V::Impl,
                    Scalar<typename V::ElementType> >
operator*(const Vector<V> &x, const typename V::ElementTypeRef alpha);

// x / alpha
template <typename V>
const VectorClosure<OpMult,
                    typename V::Impl,
                    Scalar<typename V::ElementType> >
operator/(const Vector<V> &x, const typename V::ElementTypeRef alpha);

// x + y
template <typename V1, typename V2>
const VectorClosure<OpAdd, typename V1::Impl, typename V2::Impl>
operator+(const Vector<V1> &x, const Vector<V2> &y);

// x - y
template <typename V1, typename V2>
const VectorClosure<OpSub, typename V1::Impl, typename V2::Impl>
operator-(const Vector<V1> &x, const Vector<V2> &y);

// x * y' TO BE GENERALIZED
template <typename V1, typename V2>
const typename Promotion<typename Vector<V1>::ElementType,
                         typename Vector<V2>::ElementType>::Type
operator*(const Vector<V1> &x, const Vector<V2> &y);

// A * x
template <typename M, typename V>
const VectorClosure<OpMult, typename M::Impl, typename V::Impl>
operator*(const Matrix<M> &A, const Vector<V> &x);

// x * A
template <typename M, typename V>
const VectorClosure<OpMult, typename V::Impl, typename M::Impl>
operator*(const Vector<V> &x, const Matrix<M> &A);

//-- MatrixClosure -------------------------------------------------------------

// -A
template <typename M>
const MatrixClosure<OpMult,
                    typename M::Impl,
                    Scalar<typename M::ElementType> >
operator-(const Matrix<M> &x);

// alpha * A
template <typename M>
const MatrixClosure<OpMult,
                    typename M::Impl,
                    Scalar<typename M::ElementType> >
operator*(const typename M::ElementTypeRef alpha, const Matrix<M> &A);

// A * alpha
template <typename M>
const MatrixClosure<OpMult,
                    typename M::Impl,
                    Scalar<typename M::ElementType> >
operator*(const Matrix<M> &A, const typename M::ElementTypeRef alpha);

// A / alpha
template <typename M>
const MatrixClosure<OpMult,
                    typename M::Impl,
                    Scalar<typename M::ElementType> >
operator/(const Matrix<M> &A, const typename M::ElementTypeRef alpha);

// A + B
template <typename M1, typename M2>
const MatrixClosure<OpAdd, typename M1::Impl, typename M2::Impl>
operator+(const Matrix<M1> &A, const Matrix<M2> &B);

// A - B
template <typename M1, typename M2>
const MatrixClosure<OpSub, typename M1::Impl, typename M2::Impl>
operator-(const Matrix<M1> &A, const Matrix<M2> &B);

// A^T
template <typename M>
const MatrixClosure<OpTrans,
                    typename M::Impl,
                    typename M::Impl>
transpose(const Matrix<M> &A);

// A^H
template <typename M>
const MatrixClosure<OpConjTrans,
                    typename M::Impl,
                    typename M::Impl>
conjugateTranspose(const Matrix<M> &A);

// A * B
template <typename M1, typename M2>
const MatrixClosure<OpMult, typename M1::Impl, typename M2::Impl>
operator*(const Matrix<M1> &A, const Matrix<M2> &B);

} // namespace flens

#include <flens/matvecoperations.tcc>

#endif // FLENS_MATVECOPERATIONS_H
