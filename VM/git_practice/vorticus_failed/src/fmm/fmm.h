/*
 * Copyright (C) 2017 Matthias Kirchhart
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
#ifndef FMM_FMM_H
#define FMM_FMM_H

#include "fmm/tree.h"

////////////////////////////
// ‘Public’ interface     //
////////////////////////////


namespace fmm
{

template <typename strategy, typename source_iter, typename target_iter>
void fmm( strategy& S, source_iter sbegin, source_iter send,
                       target_iter tbegin, target_iter tend );
}








/////////////////////////////////
// Implementation interface.   //
/////////////////////////////////


namespace fmm
{

template <typename strategy, typename source_box>
void upward_pass( source_box *n, strategy *S );

template <typename strategy, typename source_box, typename target_box>
void dual_tree_traversal( target_box *A, source_box *B, strategy *S );

template <typename strategy, typename source_box, typename target_box>
void interact( target_box *A, source_box *B, strategy *S );

template <typename strategy, typename target_box>
void downward_pass( target_box *n, strategy *S );

}




/////////////////////////
//  Implementation.    //
/////////////////////////


namespace fmm
{

template <typename strategy, typename source_iter, typename target_iter>
void fmm( strategy& S, source_iter sbegin, source_iter send,
                       target_iter tbegin, target_iter tend )
{
    using Mcoeff_t = typename strategy::Mcoeff_t;
    using Lcoeff_t = typename strategy::Mcoeff_t;
    using source_tree = tree<Mcoeff_t,source_iter>;
    using target_tree = tree<Lcoeff_t,target_iter>;

    if ( sbegin == send || tbegin == tend ) return;

    source_tree sources { sbegin, send, S, S.source_leaf_max() };
    target_tree targets { tbegin, tend, S, S.target_leaf_max() };

    #pragma omp parallel
    #pragma omp single
    {
        upward_pass        ( sources.root.get(), &S );
        dual_tree_traversal( targets.root.get(), sources.root.get(), &S );
        downward_pass( targets.root.get(), &S );
    }
}

template <typename strategy, typename source_box>
void upward_pass( source_box *n, strategy *S )
{
    if ( n->is_leaf() )
    {
        S->p2m( n->coeffs, n->centre, n->begin, n->end );
    }
    else
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                #pragma omp task untied
                upward_pass( n->children[ i ].get(), S );
            }
        }
        #pragma omp taskwait

        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                S->m2m( n->coeffs,  n->children[ i ]->coeffs,
                        n->centre - n->children[ i ]->centre );
            }
        }
    }
}

template <typename strategy, typename source_box, typename target_box>
void dual_tree_traversal( target_box *A, source_box *B, strategy *S )
{
    if ( ( A->length > B->length && ! A->is_leaf() ) ||
         ( ! A->is_leaf() && B->is_leaf() ) )
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( A->children[ i ] )
            {
                #pragma omp task untied
                interact( A->children[ i ].get(), B, S );
            }
        }
        #pragma omp taskwait
    }
    else if ( ! B->is_leaf() )
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( B->children[ i ] )
            {
                interact( A, B->children[ i ].get(), S );
            }
        }
    }
    else interact( A, B, S );
}

template <typename strategy, typename source_box, typename target_box>
void interact( target_box *A, source_box *B, strategy *S )
{
    using std::min;

    if ( S->mac( A->centre, A->length, B->centre, B->length ) )
    {
        const size_t n = A->end - A->begin;
        const size_t m = B->end - B->begin;

        const real t_m2l = S->time_m2l();
        const real t_m2p = S->time_m2p( n );
        const real t_p2l = S->time_p2l( m );
        const real t_p2p = S->time_p2p( n, m );
        const real t_min = min( min(t_m2l,t_m2p), min(t_p2l,t_p2p) );

        if      ( t_m2l == t_min )
        {
            S->m2l( A->coeffs, B->coeffs, A->centre - B->centre );
        }
        else if ( t_m2p == t_min )
        {
            S->m2p( A->begin, A->end, B->coeffs, B->centre );
        }
        else if ( t_p2l == t_min )
        {
            S->p2l( A->coeffs, A->centre, B->begin, B->end );
        }
        else
        {
            S->p2p( A->begin, A->end, B->begin, B->end );
        }
    }
    else if ( A->is_leaf() && B->is_leaf() )
    {
        S->p2p( A->begin, A->end, B->begin, B->end );
    }
    else
    {
        dual_tree_traversal( A, B, S );
    }
}

template <typename strategy, typename target_box>
void downward_pass( target_box *n, strategy *S )
{
    if ( n->is_leaf() )
    {
        S->l2p( n->begin, n->end, n->coeffs, n->centre );
    }
    else
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                S->l2l( n->children[ i ]->coeffs, n->coeffs, n->children[ i ]->centre - n->centre );
            }
        }

        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                #pragma omp task untied
                downward_pass( n->children[ i ].get(), S );
            }
        }
        #pragma omp taskwait
    }
}

}

#endif

