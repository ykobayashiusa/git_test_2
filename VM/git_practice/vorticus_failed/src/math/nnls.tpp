/*
 * Copyright (C) 2015 Matthias Kirchhart
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

namespace math
{

template <typename flt>
arma::Col<flt> nnls( const arma::Mat<flt> &A, const arma::Col<flt> &b, flt tol )
{
    using idx_t = arma::uword;
    const idx_t m { A.n_rows };
    const idx_t n { A.n_cols };

    // Initialisation.
    arma::uvec P( n ); P.zeros();
    arma::uvec R( n ); R.ones();

    arma::Col<flt> x( n ); x.zeros();
    arma::Col<flt> s( n ); s.zeros();
    arma::Col<flt> w = A.t()*( b - A*x );
    idx_t j; w.max(j);

    while ( max(R) == 1 && w(j) > tol )
    {
        P(j) = 1;
        R(j) = 0;

        arma::uvec Pvec = find( P == 1 );

        s.zeros();
        s(Pvec) = solve( A.cols(Pvec), b );

        while ( min(s(Pvec)) <= 0 )
        {
            arma::uvec tmpvec = find( P == 1 && s <= 0 );
            flt alpha = min( x(tmpvec)/(x(tmpvec)-s(tmpvec)) );
            x += alpha*(s-x);

            tmpvec = find( abs(x) < tol );
            P( tmpvec ).zeros();
            R( tmpvec ).ones();

            Pvec = find( P == 1 );
            s.zeros();
            s(Pvec) = solve( A.cols(Pvec), b );
        }

        x = s;
        w = A.t()*(b-A*x); 
        w.max(j);
    }

    return x;
}

}

