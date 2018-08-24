/*
 * Copyright (C) 2016 Matthias Kirchhart
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
 * details.  *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#include "math/rev_simplex.h"

#include <memory>
#include <limits>
#include <cstring>

#include "math/lu.h"

namespace
{

double dot( int n, const double *a, const double *b );

int simplex_one( int m, int n, const double *A, const double *b,
                 int *B, int *NB, double *solution, double &z,
                 double *LU, int *P,
                 double *d, double *alpha, double *beta, double *pi );

double compute_target_value_one( int m, int n, int *B, const double *beta );

void compute_simplex_multiplier_one( int m, int n, double *pi, int *B,
                                     const double *LU, const int *P );

void compute_reduced_cost_one( int m, int n, double *d, const double *pi,
                               const double *A, const int *NB );

void set_basis_col_one( int m, int n, double *v, const double *A, const int j );
void set_basis_matrix_one( int m, int n, double *Bmat, const double *A, const int *B );
}

namespace math
{

int rev_nnsolve( int m, int n, const double *A, const double *b, double *result )
{
    // Statically allocated memory for small-sized problems.
    // Currently supported size: m <= 16, n <= 256.
    double    LU_static[ 256 ]; 
    double     d_static[ 256 ];
    int       NB_static[ 256 ];
    double alpha_static[  16 ];
    double  beta_static[  16 ];
    double    pi_static[  16 ];
    int        B_static[  16 ];
    int        P_static[  16 ];

    // Dynamically allocated memory for bigger problems.
    std::unique_ptr<double[]>    LU_dynamic; if ( m > 16  )        LU_dynamic.reset(new double[ m*m ]);
    std::unique_ptr<double[]>     d_dynamic; if ( n > 256 )         d_dynamic.reset(new double[ n ]);
    std::unique_ptr<int[]>       NB_dynamic; if ( n > 256 )        NB_dynamic.reset(new int   [ n ]);
    std::unique_ptr<double[]> alpha_dynamic; if ( m > 16 )      alpha_dynamic.reset(new double[ m ]);
    std::unique_ptr<double[]>  beta_dynamic; if ( m > 16 )       beta_dynamic.reset(new double[ m ]);
    std::unique_ptr<double[]>    pi_dynamic; if ( m > 16 )         pi_dynamic.reset(new double[ m ]);
    std::unique_ptr<int[]>        B_dynamic; if ( m > 16 )          B_dynamic.reset(new int   [ m ]);
    std::unique_ptr<int[]>        P_dynamic; if ( m > 16 )          P_dynamic.reset(new int   [ m ]);

    // Classic pointer to pass on to worker routine.
    double    *LU    = ( m <= 16  )     ?    LU_static :    LU_dynamic.get();
    double    *d     = ( n <= 256 )     ?     d_static :     d_dynamic.get();
    int       *NB    = ( n <= 256 )     ?    NB_static :    NB_dynamic.get();
    double    *alpha = ( m <= 16 )      ? alpha_static : alpha_dynamic.get();
    double    *beta  = ( m <= 16 )      ?  beta_static :  beta_dynamic.get();
    double    *pi    = ( m <= 16 )      ?    pi_static :    pi_dynamic.get();
    int       *B     = ( m <= 16 )      ?     B_static :     B_dynamic.get();
    int       *P     = ( m <= 16 )      ?     P_static :     P_dynamic.get();
    
    for ( int i = 0; i < n; ++i ) NB[i] = i;
    for ( int i = 0; i < m; ++i )  B[i] = i + n;

    double z;
    return simplex_one( m, n, A, b, B, NB, result, z, LU, P, d, alpha, beta, pi );
}

}


namespace
{

int simplex_one( int m, int n, const double *A, const double *b,
                 int *B, int *NB, double *solution, double &z,
                 double *LU, int *P,
                 double *d, double *alpha, double *beta, double *pi )
{
    std::fill( LU, LU + m*m, 0. );
    for ( int i = 0; i < m; ++i )
    {
        P [ i ] = i + 1;
        LU[ i + m*i ] = 1;
    }

    for ( int count = 0; count < 200; ++count )
    {
        // Current basic solution and value of target functional.
        std::memcpy( beta, b, sizeof(double)*m );
        math::lu_solve( m, LU, P, beta );                // β = inv(Bmat)*b.
        z = compute_target_value_one( m, n, B, beta );   // z = dot( c(B), β );

        // Pricing.
        compute_simplex_multiplier_one( m, n, pi, B, LU, P ); // π = inv(Bmat.t())*c(B);
        compute_reduced_cost_one( m, n, d, pi, A, NB );       // d = c(NB) - A.cols(NB).t()*π;

        // Dantzig's rule.
        int q = 0; double dmax = d[0]; 
        for ( int i = 1; i < n; ++i )
        {
            if ( dmax < d[ i ] )
            {
                   q = i;
                dmax = d[ i ];
            }
        }

        if ( dmax <= 1e-9 )
        {
            std::fill( solution, solution + n, 0. );
            for ( int i = 0; i < m; ++i )
            {
                if ( B[i] < n )
                    solution[ B[i] ] = beta[i];
            } 
            return 0; // Optimal solution has been found.
        }

        // Ratio test.
        set_basis_col_one( m, n, alpha, A, NB[q] );
        math::lu_solve( m, LU, P, alpha ); // α = inv(Bmat)*A(NB(q));

        double theta =  std::numeric_limits<double> ::max();
        int     p    =  std::numeric_limits<int>    ::max();
        for ( int i = 0; i < m; ++i )
        {
            if ( alpha[i] >= 1e-9 && theta > beta[i]/alpha[i] )
            {
                p = i;
                theta = beta[i]/alpha[i];
            }
        }

        if ( p == std::numeric_limits<int>::max() )
        {
            // All alpha zero or negative.
            return 1; // Unbounded.
        }

        // Update.
        int real_p =  B[p];
        int real_q = NB[q];
        NB[q] = real_p;
         B[p] = real_q;

        set_basis_matrix_one( m, n, LU, A, B );
        math::lu( m, LU, P );
    }
    return 2;
}

/*!
 * \brief Essentially returns dot( c(B), β );
 *
 * This functions returns the current value of the target functional. In a
 * phase I problem the cost-vector is implicitly given, however, so the general
 * call dot( c(B), beta ) does not work in this case.
 */
double compute_target_value_one( int m, int n, int *B, const double *beta )
{
    double result { 0 };
    for ( int i = 0; i < m; ++i )
    {
        if ( B[i] >= n )
            result -= beta[i];
    }

    return result;
}

/*!
 * \brief Essentially computes π = inv(Bmat^T)*c(B).
 *
 * Computes the simplex multiplier π, takes care of implicitly given cost-vector
 * in a phase I problem.
 */
void compute_simplex_multiplier_one( int m, int n, double *pi, int *B,
                                     const double *LU, const int *P )
{
    for ( int i = 0; i < m; ++i )
    {
        if ( B[i] < n ) pi[i] =  0;
        else            pi[i] = -1;
    }
    math::lu_solve_t( m, LU, P, pi );
}

/*!
 * \brief Essentially computes d = c(NB) - A.cols(NB).t()*π;
 *
 * Computes the values of reduced cost d, takes care of implicitly given cost
 * vector in phase I problem.
 */
void compute_reduced_cost_one( int m, int n, double *d, const double *pi,
                               const double *A, const int *NB )
{
    for ( int i = 0; i < n; ++i )
    {
        if ( NB[i] < n )
        {
            d[i] = - dot( m, pi, A + m*(NB[i]) );
        }
        else
        {
            d[i] = -1 - pi[ NB[i] - n ];
        }
    }
}

/*!
 * \brief Computes the dot-product between two vectors.
 */
inline
double dot( int n, const double *a, const double *b )
{
    double result = 0;
    for ( int i = 0; i < n; ++i )
        result += a[i]*b[i];
    return result;
}

/*!
 * \brief Essentially performs v = A.col(j).
 *
 * This function sets v equal to the specified column of A. In a phase I
 * problem, however, the identity matrix of A is not explicitly stored. We thus
 * need to check for each column if it is part of the non-stored identity
 * matrix.
 */
void set_basis_col_one( int m, int n, double *v, const double *A, const int j )
{
    if ( j < n )
    {
        // j is a normal variable.
        std::memcpy( v, A + m*j, m*sizeof(double) );
    }
    else
    {
        std::fill( v, v + m, 0. );
        v[ j - n ] = 1;
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
void set_basis_matrix_one( int m, int n, double *Bmat, const double *A, const int *B )
{
    for ( int j = 0; j < m; ++j )
    {
        if ( B[j] < n )
        {
            // B[j] is a normal variable.
            std::memcpy( Bmat + m*j, A + m*(B[j]), m*sizeof(double) ); 
        }
        else
        {
            // B[j] is an artificial/slack variable.
            std::fill( Bmat + m*j, Bmat + m*(j+1), 0. );
            Bmat[ (B[j] - n) + m*j ] = 1;
        }
    }
}

}

