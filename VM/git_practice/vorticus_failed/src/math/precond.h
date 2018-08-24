/*
 * Copyright (C) 2016 Matthias Kirchhart
 *
 * This file is part of vorticus.
 *
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
#ifndef MATH_PRECOND_H
#define MATH_PRECOND_H

#include "math/spmat.h"

namespace math
{

struct diag_precond
{
    const vector diag;

    vector operator*( const vector &rhs ) const
    {
        return rhs / diag;
    }
};

template <typename Mat>
struct jacobi_precond
{
    const Mat    &A;
    const vector &diag;
    const real   omega    { 1 };
    const size_t max_iter { 1 };

    vector operator*( const vector &rhs )
    {
        vector x ( diag.size(), arma::fill::zeros );
        for ( size_t iter = 0; iter < max_iter; ++iter )
        {
            x = (1-omega)*x + omega*(rhs - (A*x - diag % x)) / diag;
        }
        return x;
    }
};

/*
struct jacobi_precond
{
    const crs_matrix &A;
    const real omega { 1 };
    const size_t max_iter { 2 };

    vector operator*( const vector &rhs )
    {
        vector x     ( A.rows(), arma::fill::zeros ),
               x_prev( A.rows(), arma::fill::zeros );

        for ( size_t iter = 0; iter < max_iter; ++iter)
        {
            #pragma omp parallel for schedule(dynamic)
            for ( size_t i = 0; i < rhs.size(); ++i )
            {
                real sigma = 0; 
                real diag  = 0;
                for ( size_t k = A.row_ptr_[i]; k < A.row_ptr_[i+1]; ++k )
                {
                    size_t j = A.col_idx_[k];
                    real   v = A.values_[k];
                    if ( i == j ) diag = v;
                    else sigma = sigma + v*x_prev(j);
                }
                sigma = (rhs(i)-sigma)/diag; 
                x(i) = omega*sigma + (1-omega)*x_prev(i);
            }
            x_prev = x;
        }
        return std::move(x);
    }
};

struct ssor_precond
{
    const crs_matrix &A;
    const real omega { 1.6 };
    const size_t max_iter { 2 };

    vector operator*( const vector &rhs )
    {
        vector xhalf( A.rows(), arma::fill::zeros ),
               x    ( A.rows(), arma::fill::zeros );


        for ( size_t iter = 0; iter < max_iter; ++iter )
        {

        for ( size_t i = 0; i < A.rows(); ++i )
        {
            real sigma = 0;
            real diag = 0;
            for ( size_t k = A.row_ptr_[i]; k < A.row_ptr_[i+1]; ++k )
            {
                size_t j = A.col_idx_[k];
                real   v = A.values_[k];
                if ( j == i ) diag = v;
                else sigma = (j<i) ? (sigma + v*xhalf(j)) : (sigma + v*x(j));
            }
            sigma = ( rhs(i) - sigma ) / diag;
            xhalf(i) = omega*sigma + (1-omega)*x(i);
        }

        for ( size_t i = A.rows(); i-- > 0; )
        {
            real sigma = 0;
            real diag  = 0;
            for ( size_t k = A.row_ptr_[i]; k < A.row_ptr_[i+1]; ++k )
            {
                size_t j = A.col_idx_[k];
                real   v = A.values_[k];
                if ( j == i ) diag = v;
                else sigma = (j<i) ? (sigma + v*xhalf(j)) : (sigma + v*x(j));
            }
            sigma = ( rhs(i) - sigma ) / diag;
            x(i) = omega*sigma + (1-omega)*xhalf(i);
        }

        }
        return std::move(x);
    }
};
*/

}

#endif

