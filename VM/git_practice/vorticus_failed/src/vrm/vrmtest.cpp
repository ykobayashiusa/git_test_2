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


#include "fmm/fmm.h"
#include "misc/stopwatch.h"
#include "math/constants.h"

#include "vrm/particle.h"
#include "vrm/accessor.h"
#include "vrm/laplacian.h"
#include "vrm/particle_strategy.h"

#include "misc/stable_vector.h"
#include "geometry/mesh_index.h"

#include <sstream>

namespace vrm
{

template <uint stages>
using field_t  = misc::stable_vector<particle<stages>>;

template <uint stages>
void simulate( field_t<stages>& field, const real nu,
               const real h );

template <uint stages>
void biot_savart2d( field_t<stages>& field, accessor<stages> get, real epsilon );

template <uint stages>
void print_statistics( field_t<stages>& field, real nu, real t  );

template <uint stages>
void error_estimate( field_t<stages>& field, real nu, real t, real epsilon  );

template <uint stages>
void print( field_t<stages>& field, real nu, real t, real epsilon, int count );

void heat_test( real h )
{
    const real nu = 1./50.;

    field_t<4> field;

    particle<4> p;
    p.G.z = 2*math::pi;
    field.push_back( p );

    simulate( field, nu, h );
}


template <uint stages>
void simulate( field_t<stages>& field, const real nu, const real h )
{
    using geometry::point;
    using iterator = typename field_t<stages>::iterator;

    real t = 0; 
    const real r      = 0.5;
    const real max_dt = std::min( 0.125*((r*h)*(r*h)/(4*nu)), 0.01 );
    const real   T    = 1;

    accessor<stages> access;
    //const real epsilon = std::pow( h, 1./3.5 ) * ( 0.1 / std::pow( 0.1, 1./3.5 ) );
    const real epsilon = 3*h;

    int n = 0; int count = 0;
    while ( t < T )
    {
        real dt = max_dt;
        if ( (t + dt) > T ) dt = T - t;
        access.dt(dt);

        bool print_data = ( n++ % 50 ) == 0;
        if ( print_data )
        {
            std::cout << "t = " << t << ", N = " << field.size() << ".\n";
            print_statistics(field, nu, t );
            print( field, nu, t, epsilon, count++ );
        }

        real time_biot_savart  = 0;
        real time_lapl_matrix  = 0;
        for ( size_t s = 0; s < stages; ++s )
        {
            access.stage(s);
            geometry::mesh_index<iterator,accessor<stages>> find( field.begin(), field.end(), 3*h, access );
            auto make = [&field,&find]( geometry::point pos ) -> iterator
            {
                particle<stages> p; p.x = pos;
                field.push_back(p);
                find.insert(field.end() - 1);
                return field.end() - 1;
            };


            misc::stopwatch clock;
            stiffness_matrix<iterator> A = laplacian_2d( field.begin(), field.end(), access, find, make, h, 1 );
            time_lapl_matrix += clock.elapsed() / stages;

            for ( auto &entry: A )
            {
                iterator i = std::get<0>(entry);
                iterator j = std::get<1>(entry);
                real     v = std::get<2>(entry);
                if ( j != iterator {} )
                    access.dG(*j) += (nu*v)*access.G(*i);
            }

            clock.reset();
            biot_savart2d( field, access, epsilon );
            time_biot_savart += clock.elapsed() / stages;
            
            if ( s == 0 )
            {
                real max_u = 0;
                for ( particle<stages> &p: field )
                {
                    max_u = std::max( max_u, access.u(p).r() );
                }
                dt = std::min( dt, (0.125*h)/max_u );
                access.dt(dt);
            }
        }
        if ( print_data ) std::cout << "Time for Laplacian:   " << time_lapl_matrix << std::endl;
        if ( print_data ) std::cout << "Time for Biot-Savart: " << time_biot_savart << std::endl;
        if ( print_data ) std::cout << "Ratio: " << time_biot_savart / time_lapl_matrix << std::endl;

        for ( particle<stages> &p: field )
            access.advance(p);

        t += dt;
        if ( print_data ){ std::cout << "\n\n"; std::cout.flush(); }
    }
    
    std::cout << "t = " << t << ", N = " << field.size() << ".\n";
    print_statistics(field, nu, t );
    error_estimate( field, nu, t, epsilon );
    print( field, nu, t, epsilon, count );
}

template <uint stages>
void biot_savart2d( field_t<stages>& field, accessor<stages> get, real epsilon )
{
    particle_strategy<16,accessor<stages>> S(get); S.epsilon = epsilon; S.theta_max = 0.8;
    fmm::fmm( S, field.begin(), field.end(),
                 field.begin(), field.end() );
}

template <uint stages>
void print_statistics( field_t<stages>& field, real nu, real t  )
{
    using geometry::point;
    constexpr real MAX =  std::numeric_limits<real>::max();
    constexpr real MIN = -std::numeric_limits<real>::max();

    real  circulation      {};
    point  linear_momentum {};
    real  angular_momentum {};
    point min { MAX, MAX, MAX };
    point max { MIN, MIN, MIN };
    for ( particle<stages> p: field )
    {
         circulation     += p.G.z;
         linear_momentum += p.x*p.G.z;
        angular_momentum += scal_prod( p.x, p.x )*p.G.z;
        min = min_components( min, p.x );
        max = max_components( max, p.x );
    }
    std::cout << "Circulation:      " <<      circulation << std::endl;
    std::cout << "Linear Momentum:  " <<  linear_momentum << std::endl;
    std::cout << "Angular Momentum: " << angular_momentum << std::endl;
    std::cout << "Angular Momentum Error: " << (4*2*math::pi*t*nu) - angular_momentum << std::endl;
    std::cout << "Bounding box: " << min << ' ' << max << std::endl;
}

template <uint stages>
void error_estimate( field_t<stages>& field, real nu, real t, real epsilon  )
{
    size_t  N = 1000;
    real h = 3./N;
    std::cout << "h-Error estimate: " << h << ", epsilon: " << epsilon << ".\n";
    misc::stable_vector<particle<stages>> probes( N * N );

    #pragma omp parallel for
    for ( size_t i = 0; i < N; ++i )
    for ( size_t j = 0; j < N; ++j )
    {
        particle<stages> &p = probes[ i*N + j ];
        p.x.x = -1.5 + h/2 + i*h;
        p.x.y = -1.5 + h/2 + j*h;
        p.k_x[ 0 ] = geometry::point {};
    }

    accessor<stages> get; get.stage(0);
    particle_strategy<16,accessor<stages>> S(get); S.epsilon = epsilon; S.theta_max = 0;
    fmm::fmm( S, field.begin(), field.end(),
                 probes.begin(), probes.end() );

    real err = 0, norm = 0;
    #pragma omp parallel for reduction(+:err,norm)
    for ( size_t i = 0; i < N*N; ++i )
    {
        geometry::point u_true, x = probes[ i ].x;
        if ( x.r() != 0 )
        {
            u_true = ( -std::expm1( -x.r2() / (4*nu*t) ) / x.r() ) * x.e_phi();
        }
        geometry::point u_probe = probes[ i ].k_x[0];
        err  += h*h*(u_true - u_probe).r2();
        norm += h*h*(u_true).r2();
    }
    err = std::sqrt( err ); norm = std::sqrt( norm );
    std::cout << "Approximate relative L_2-Error: " << err/norm << std::endl;
}

template <uint stages>
void print( field_t<stages>& field, real nu, real t, real epsilon, int count )
{
    accessor<stages> get; get.stage(0);
    particle_strategy<16,accessor<stages>> S(get); S.epsilon = epsilon; S.theta_max = 0;
    fmm::fmm( S, field.begin(), field.end(),
                 field.begin(), field.end() );

    std::stringstream filename; filename << "output" << count << ".vtp";
    std::ofstream str( filename.str() );
    str << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    str << "<PolyData>\n";
    str << "<Piece NumberOfPoints=\"" << field.size() << "\">\n";
    str << "<Points>\n";
    str << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\">\n";
    for ( size_t i = 0; i < field.size(); ++i ) str << field[ i ].x << "\n";
    str << "</DataArray>\n";
    str << "</Points>\n";
    str << "<PointData>\n";
    str << "<DataArray type=\"Float32\" Name=\"G\" NumberOfComponents=\"3\">\n";
    for ( size_t i = 0; i < field.size(); ++i ) str << field[ i ].G << "\n";
    str << "</DataArray>\n";
    str << "<DataArray type=\"Float32\" Name=\"u\" NumberOfComponents=\"3\">\n";
    for ( size_t i = 0; i < field.size(); ++i ) str << field[ i ].k_x[ 0 ] << "\n";
    str << "</DataArray>\n";
    str << "<DataArray type=\"Float32\" Name=\"err\" NumberOfComponents=\"3\">\n";
    for ( size_t i = 0; i < field.size(); ++i )
    {
        geometry::point x = field[ i ].x; 
        geometry::point u_true = ( -std::expm1( -x.r2() / (4*nu*t) ) / x.r() ) * x.e_phi();
        str << field[ i ].k_x[0] - u_true << std::endl;
    }
    str << "</DataArray>\n";
    str << "<DataArray type=\"Float32\" Name=\"rel_err\" NumberOfComponents=\"1\">\n";
    for ( size_t i = 0; i < field.size(); ++i )
    {
        geometry::point x = field[ i ].x; 
        geometry::point u_true = ( -std::expm1( -x.r2() / (4*nu*t) ) / x.r() ) * x.e_phi();
        str << (field[ i ].k_x[0] - u_true).r() / u_true.r() << std::endl;
    }
    str << "</DataArray>\n";
    str << "<DataArray type=\"Float32\" Name=\"u_true\" NumberOfComponents=\"3\">\n";
    for ( size_t i = 0; i < field.size(); ++i )
    {
        geometry::point x = field[ i ].x; 
        geometry::point u_true = ( -std::expm1( -x.r2() / (4*nu*t) ) / x.r() ) * x.e_phi();
        str << u_true << std::endl;
    }
    str << "</DataArray>\n";
    str << "</PointData>\n";
    str << "</Piece>\n";
    str << "</PolyData>\n";
    str << "</VTKFile>\n" ;
}

}

int main( int argc, char *argv[] )
{
    std::stringstream str( argv[1] );
    double h; str >> h;
    vrm::heat_test( h );
} 

