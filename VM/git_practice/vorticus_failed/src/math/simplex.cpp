/*
 * Copyright (C) 2015 Matthias Kirchhart
 *
 * This file is part of vorticus.
 * vorticus is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * vorticus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#include "math/simplex.h"

#include <iostream>
#include <cblas.h>
#include <lapacke.h>

namespace math
{

namespace
{


using arma::uvec;
using arma::uword;

enum class simplex_phase { ONE, TWO };
int revised_simplex( simplex_phase phase,
                     const matrix &A, const vector &b, const vector &c,
                     uvec &B, uvec &NB, vector &solution, real &z );

void set_basis_matrix( matrix &Bmat, const matrix &A, const uvec &B );
void set_basis_col( vector &v, const matrix &A, const uword j );

void compute_simplex_multiplier( simplex_phase phase, vector &pi,
                                 const matrix &LU, const std::vector<int> &P,
                                 const vector &c,  const uvec &B );

double compute_target_value( simplex_phase phase, const vector &c,
                             const uvec &B, const vector &beta );

void compute_reduced_cost( simplex_phase phase, vector &d,
                           const vector &c, const vector &pi,
                           const matrix &A, const uvec &NB );

void decompose_lu( matrix &A, std::vector<int> &P );
void     solve_lu( const matrix &LU, const std::vector<int> &P, vector &b );

void fletcher_matthews( matrix &LU, std::vector<int> &P,
                        uword col_idx, const double *const column );
void ipiv_to_permutation( std::vector<int> &P );
void permutation_to_ipiv( std::vector<int> &P );

}


int nnsolve( const matrix &A, const vector &b, vector &result )
{
    using arma::uvec;
    using arma::span;
    using arma::uword;

    const uword m { A.n_rows };
    const uword n { A.n_cols };

    vector c( nullptr, n, false, true ); // Dummy cost-vector.

    uvec free ( n );
    uvec basis( m );
    for ( uword i = 0; i < n; ++i ) free (i) = i;
    for ( uword i = 0; i < m; ++i ) basis(i) = i + n;

    real z;
    return revised_simplex( simplex_phase::ONE, A, b, c, basis, free, result, z );
}

namespace
{

/*!
 * \brief The revised simplex method.
 * \param  A The \f$\mathbb{R}^{m\times n}\f$ matrix of coefficients.
 * \param  b The \f$\mathbb{R}^{m}_{+}\f$ right-hand side vector.
 * \param  c The \f$\mathbb{R}^{n}\f$ cost coefficient vector.
 * \param  B The list of \f$m\f$ initial     basic variable indeces.
 * \param NB The list of \f$n\f$ initial non-basic variable indeces.
 * \param solution The vector where the result will be stored.
 * \param  z The value of the corresponding cost \f$c^\topx\f$.
 * \returns 0 if the optimal solution has been found, 1 if the problem is unbounded.
 *
 * This is an implementation of the revised simplex method. Given a matrix
 * \f$ A\in\mathbb{R}^{m\times n}\f$ with \f$n >= m\f$, a right-hand side
 * vector \f$b\in\mathbb{R}^{m}_{+}\f$ and a cost vector \f$ c\in\mathbb{R}^{n}\f$,
 * this method solves the linear programming problem:
 * \f[
 * \min c^\top \vv{x}, \text{s.t.} \\
 * A\cdot x = b, \\
 * x \geq 0.
 * \f]
 * It is assumed that \f$A\f$ has full rank.
 *
 * In order to solve the problem, an initial basic feasible solution needs to
 * be given. This is done by passing basis and non-basis vectors. The non-basis
 * vector is a list of indeces of dimension \f$n\f$, where the basis-vector
 * contains of \f$m\f$ indeces. For the initial basic feasible solution, all
 * entries in the solution corresponding to non-basic variables are assumed to
 * be zero, the remaining variables are then result of a linear system of
 * equations. In order to obtain such an initial basic feasible solution, a
 * two-phase method can be employed.
 */
int revised_simplex( simplex_phase phase,
                     const matrix &A, const vector &b, const vector &c,
                     uvec &B, uvec &NB, vector &solution, real &z )
{
    using arma::dot;
    using arma::uword;

    const uword m {  B.size() };
    const uword n { NB.size() };
    
    vector alpha( m );
    vector beta ( m );
    vector pi   ( m );
    vector d    ( n );
    matrix LU   ( m, m );
    std::vector<int> P(m);

    // LU contains the LU-decomposition of Bmat := A.cols(B);
    if ( phase == simplex_phase::ONE )
    {
        // In phase I we start with an identity matrix, we thus do not need
        // to perform an LU-decomposition.
        LU.fill( arma::fill::eye );
        for ( uword i = 0; i < m; ++i )
            P[ i ] = i + 1; // Plus one due to Fortran convention.
    } 
    else
    {
        // Compute the initial LU-decomposition of the basis matrix.
        set_basis_matrix( LU, A, B );
        decompose_lu( LU, P );
    }

    while ( true )
    {
        // Current basic solution and value of target functional.
        beta = b; solve_lu( LU, P, beta );                    // β = inv(Bmat)*b.
        z = compute_target_value( phase, c, B, beta );        // z = dot( c(B), β );

        // Pricing.
        compute_simplex_multiplier( phase, pi, LU, P, c, B ); // π = inv(Bmat.t())*c(B);
        compute_reduced_cost( phase, d, c, pi, A, NB );       // d = c(NB) - A.cols(NB).t()*π;

        // Dantzig's rule.
        uword q; d.max(q); 

        if ( d(q) <= 1e-9 )
        {
            solution.resize( n ); solution.zeros();
            for ( uword i = 0; i < m; ++i )
            {
                if ( B(i) < n )
                    solution( B(i) ) = beta(i);
            }
            return 0; // Optimal solution has been found.
        }

        // Ratio test.
        set_basis_col( alpha, A, NB(q) );
        solve_lu( LU, P, alpha );         // α = inv(Bmat)*A(NB(q));

        real theta =  std::numeric_limits<real> ::max();
        uword p    =  std::numeric_limits<uword>::max();
        for ( uword i = 0; i < m; ++i )
        {
            if ( alpha(i) >= 1e-9 && theta > beta(i)/alpha(i) )
            {
                p = i;
                theta = beta(i)/alpha(i);
            }
        }

        if ( p == std::numeric_limits<uword>::max() )
        {
            // All alpha zero or negative.
            return 1; // Unbounded.
        }

        // Update.
        uword real_p =  B(p);
        uword real_q = NB(q);

        if ( m <= 10 )
        {
            // For small problem sizes refactorisation is faster.
            NB(q) = real_p;
             B(p) = real_q;
            set_basis_matrix( LU, A, B );
            decompose_lu( LU, P );
        }
        else
        {
            // Use the Fletcher-Matthews update for larger m.
            NB(q) = real_p;
            for ( uword i = p; i < m - 1; ++i )
                B(i) = B(i+1);
            B( m - 1 ) = real_q;
            set_basis_col( alpha, A, real_q );
            fletcher_matthews( LU, P, p, alpha.memptr() );
        }
    }
}

/*!
 * \brief Essentially performs Bmat = A.cols(B).
 *
 * This function sets Bmat equal to the specified columns of A. In a phase I
 * problem, however, the identity matrix of A is not explicitly stored. We thus
 * need to check for each column if it is part of the non-stored identity
 * matrix.
 */
void set_basis_matrix( matrix &Bmat, const matrix &A, const uvec &B )
{
    const uword m = A.n_rows;
    const uword n = A.n_cols;

    for ( uword j = 0; j < m; ++j )
    {
        if ( B(j) < n )
        {
            // B(j) is a normal variable.
            std::memcpy( Bmat.memptr() + m*j, A.memptr() + m*B(j), m*sizeof(double) ); 
        }
        else
        {
            // B(j) is an artificial/slack variable.
            Bmat.col(j).zeros();
            Bmat( B(j) - n, j ) = 1;
        }
    }
}

/*!
 * \brief Essentially performs v = A.col(j).
 *
 * This function sets v equal to the specified column of A. In a phase I
 * problem, however, the identity matrix of A is not explicitly stored. We thus
 * need to check for each column if it is part of the non-stored identity
 * matrix.
 */
void set_basis_col( vector &v, const matrix &A, const uword j )
{
    const uword m = A.n_rows;
    const uword n = A.n_cols;

    if ( j < n )
    {
        // j is a normal variable.
        std::memcpy( v.memptr(), A.memptr() + m*j, m*sizeof(double) );
    }
    else
    {
        v.zeros();
        v( j - n ) = 1;
    }
}

/*!
 * \brief Essentially returns dot( c(B), β );
 *
 * This functions returns the current value of the target functional. In a
 * phase I problem the cost-vector is implicitly given, however, so the general
 * call dot( c(B), beta ) does not work in this case.
 */
double compute_target_value( simplex_phase phase, const vector &c,
                             const uvec &B, const vector &beta )
{
    const uword m = B.size();
    const uword n = c.size();

    double result { 0 };
    if ( phase == simplex_phase::ONE )
    {
        for ( uword i = 0; i < m; ++i )
        {
            if ( B(i) >= n )
                result -= beta(i);
        }
    }
    else result = dot( c(B), beta );

    return result;
}

/*!
 * \brief Essentially computes π = inv(Bmat^T)*c(B).
 *
 * Computes the simplex multiplier π, takes care of implicitly given cost-vector
 * in a phase I problem.
 */
void compute_simplex_multiplier( simplex_phase phase, vector &pi,
                                 const matrix &LU, const std::vector<int> &P,
                                 const vector &c,  const uvec &B )
{
    const uword m = B.size();
    const uword n = c.size();

    if ( phase == simplex_phase::ONE )
    {
        for ( uword i = 0; i < m; ++i )
        {
            if ( B(i) < n ) pi(i) =  0;
            else            pi(i) = -1;
        }
    }
    else pi = c(B);

    LAPACKE_dgetrs( LAPACK_COL_MAJOR, 'T', m, 1, LU.memptr(), m, P.data(),
                    pi.memptr(), m );
}

/*!
 * \brief Essentially computes d = c(NB) - A.cols(NB).t()*π;
 *
 * Computes the values of reduced cost d, takes care of implicitly given cost
 * vector in phase I problem.
 */
 void compute_reduced_cost( simplex_phase phase, vector &d,
                            const vector &c, const vector &pi,
                            const matrix &A, const uvec &NB )
{
    const uword n = c.size();

    if ( phase == simplex_phase::ONE )
    {
        for ( uword i = 0; i < n; ++i )
        {
            if ( NB(i) < n )
            {
                d(i) = - dot( pi, A.col(NB(i)) );
            }
            else
            {
                d(i) = -1 - pi( NB(i) - n );
            }
        }
    }
    else d = c(NB) - A.cols(NB).t()*pi;
}

inline
void solve_lu( const matrix &LU, const std::vector<int> &P, vector &b )
{
    LAPACKE_dgetrs( LAPACK_COL_MAJOR, 'N', LU.n_cols, 1, LU.memptr(),
                    LU.n_cols, P.data(), b.memptr(), LU.n_rows );
}

inline
void decompose_lu( matrix &A, std::vector<int> &P )
{
    LAPACKE_dgetrf( LAPACK_COL_MAJOR, A.n_rows, A.n_cols, A.memptr(), A.n_rows,
                    P.data() );
}

void fletcher_matthews( matrix &LU, std::vector<int> &P,
                        uword col_idx, const double *const column )
{
    using std::abs;
    using std::max;
    using arma::span;

    const auto m = LU.n_rows;
    const auto n = LU.n_cols;
 
    vector subdiag( n - 1 );
    for ( uword c = col_idx; c < n - 1; ++c )
    {
        subdiag( c ) = LU( c+1, c+1 );
        cblas_dcopy( c+1, LU.memptr() + (c+1)*m, 1,
                          LU.memptr() + (c  )*m, 1 );
    }
    cblas_dcopy( m, column, 1,
                    LU.memptr() + (n-1)*m, 1 );

    ipiv_to_permutation( P );
    for ( uword c = col_idx; c < n - 1; ++c )
    {
        const double r   = subdiag(c) / LU(c,c);
        const double D   = subdiag(c) + LU(c+1,c)*LU(c,c);
        const double b11 = LU(c,c)/D;
        const double b21 = subdiag(c)/D;
        const double l21 = LU(c+1,c);

        bool no_permute = max(1.,abs(r)) <= max(1.,abs(LU(c+1,c))) * 
                                            max(abs(b11),abs(b21));

        if ( no_permute )
        {
            // Update L.
            LU(c+1,c) += r;
            cblas_daxpy( m - c - 2,  r, LU.memptr() + (c+1)*m + c + 2, 1,
                                        LU.memptr() + (c  )*m + c + 2, 1 );

            // Update U.
            cblas_daxpy( n - c - 2, -r, LU.memptr() + (c+1)*m + c,     m,
                                        LU.memptr() + (c+1)*m + c + 1, m );
        }
        else
        {
            // Update L.
            LU(c+1,c) = LU(c,c)/D;

            for ( uword i = c + 2; i < m; ++i )
            {
                // In-place! Unlike DGEMM.
                double tmp = LU(i,c+1);
                LU(i,c+1) = LU(i,c) - l21*tmp;
                LU(i,c  ) = tmp + b11*LU(i,c+1);
            }

            cblas_dswap( c, LU.memptr() + c,     m,
                            LU.memptr() + c + 1, m );
            std::swap( P[c], P[c+1] );

            // Update U.
            LU(c,c) = D;
            for ( uword j = c + 1; j < n - 1; ++j )
            {
                // In-place! Unlike DGEMM.
                double tmp = LU(c,j);
                LU(c  ,j) = LU(c+1,j)+l21*tmp;
                LU(c+1,j) = tmp - b11*LU(c,j);
            }
        }
    }
    permutation_to_ipiv( P );
    
    for ( uword i = 0; i < P.size(); ++i )
    {
        // P[i] *minus one* due to Fortran convention.
        std::swap( LU(i,n-1), LU( P[i] - 1, n-1 ) );
    }
    cblas_dtrsv( CblasColMajor, CblasLower, CblasNoTrans, CblasUnit,
                 m, LU.memptr(), m, LU.memptr() + (n-1)*m, 1 );
}

void ipiv_to_permutation( std::vector<int> &P )
{
    std::vector<int> tmp( P.size() );
    for ( uword i = 0; i < P.size(); ++i ) tmp[ i ] = i;
    for ( uword i = 0; i < P.size(); ++i ) std::swap( tmp[ i ], tmp[ P[i] - 1 ] );
    P.swap(tmp);
}

void permutation_to_ipiv( std::vector<int> &P )
{
    // Convert PP to P.
    std::vector<int>     PP ( P.size() );
    std::vector<int>  invPP ( P.size() );
    std::vector<int>  result( P.size() );
    for ( uword i = 0; i < P.size(); ++i )
    {
        invPP[ i ] = PP[ i ] = i;
    }

    for ( uword i = 0; i < P.size(); ++i )
    {
        result[ i ] = invPP[ P[ i ] ] + 1;
        int   i_row =    PP[ i ];
        std::swap( invPP[ i_row ], invPP[ P[ i ] ] );
        std::swap( PP[ i ],  PP[ result[ i ] - 1 ] );
    }
    P.swap(result);
}

} // End anonymous namespace.

} // End namespace math.

