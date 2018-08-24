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
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * vorticus; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */

namespace fem
{

template <size_t order>
vtkwriter<order>::vtkwriter( const_grid_iterator<order> begin, 
                             const_grid_iterator<order> end ):
numb_ { begin, end }
{}

template <size_t order> inline
void vtkwriter<order>::register_scalar( const scalar &f, std::string name )
{
    scalars_[ name ] = &f;
}

template <size_t order> inline
void vtkwriter<order>::register_vector( const vector &v, std::string name )
{
    vectors_[ name ] = &v;
}

template <size_t order>
void vtkwriter<order>::write( std::ostream& out )
{
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << numb_.size() << "\" "
        << "NumberOfCells=\"" << numb_.cell_nodes().size() << "\">\n";

    out << "<Points>\n";
    out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

    for ( point p: numb_.nodes() ) out << p << '\n';

    out << "</DataArray>\n";
    out << "</Points>\n";

    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
    for ( auto nums: numb_.cell_nodes() )
    {
        out << nums.second[0] << ' '
            << nums.second[1] << ' '
            << nums.second[2] << ' '
            << nums.second[3] << '\n';
    }
    out << "</DataArray>\n";

    out << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
    for ( size_t i { 1 }; i <= numb_.cell_nodes().size(); ++i ) out << 4*i << '\n';
    out << "</DataArray>\n";


    out << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
    for ( size_t i { 0 }; i < numb_.cell_nodes().size(); ++i ) out << 10 << '\n';
    out << "</DataArray>\n";

    out << "</Cells>\n";

    if ( scalars_.size() != 0 || vectors_.size() != 0)
    {

    out << "<PointData";
    if ( scalars_.size() != 0 ) out << " Scalars=\"" << scalars_.begin()->first << "\"";
    if ( vectors_.size() != 0 ) out << " Vectors=\"" << vectors_.begin()->first << "\"";
    out << ">\n";

    for ( const auto &s: scalars_ )
    {
        out << "<DataArray type=\"Float32\" Name=\""
            << s.first << "\" Format=\"ascii\">\n";

        for ( point p: numb_.nodes() )
            out << s.second->eval( p ) << "\n";

        out << "</DataArray>\n";
    }
    for ( const auto &v: vectors_ )
    {
        out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\""
            << v.first << "\" Format=\"ascii\">\n";

        for ( point p: numb_.nodes() )
            out << v.second->eval( p ) << "\n";

        out << "</DataArray>\n";
    }
    out << "</PointData>\n";
    }

    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

}

