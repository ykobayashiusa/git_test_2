/*
 * Copyright (C) 2014 Matthias Kirchhart
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

namespace bem
{

template <uint order>
vtkwriter<order>::vtkwriter( const typename multigrid<order>::grid& triang ):
triang_ { triang }
{
    for ( const auto& t: triang_ )
    {
        node_numbers_[ t.get_node(0) ] = 0; 
        node_numbers_[ t.get_node(1) ] = 0; 
        node_numbers_[ t.get_node(2) ] = 0; 
    }

    size_t count { 0 };
    for ( auto& no: node_numbers_ )
    {
        no.second = count++;
        node_numbers_inv_[ no.second ] = no.first;
    }
}

template <uint order>
void vtkwriter<order>::register_scalar( string name, const scalar_function& f )
{
    scalars_[ name ] = &f;
}

template <uint order>
void vtkwriter<order>::register_vector( string name, const vector_function& f )
{
    vectors_[ name ] = &f;
}

template <uint order>
void vtkwriter<order>::write( std::ostream& out )
{
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << node_numbers_.size() << "\" "
        << "NumberOfCells=\"" << triang_.size() << "\">\n";

    if ( scalars_.size() || vectors_.size() )
    {
        out << "<PointData ";
        if ( scalars_.size() ) out << "Scalars=\"" << scalars_.begin()->first <<"\" ";
        if ( vectors_.size() ) out << "Vectors=\"" << vectors_.begin()->first <<"\" ";
        out  << ">\n";

        for ( const auto& pair: scalars_ )
        {
            string name { pair.first };
            const scalar_function *f { pair.second };

            out << "<DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
            for ( size_t i { 0 }; i < node_numbers_inv_.size(); ++i )
            {
                out << f->eval( node_numbers_inv_[ i ] ) << '\n';
            } 
            out << "</DataArray>\n";
        }

        for ( const auto& pair: vectors_ )
        {
            string name { pair.first };
            const vector_function *f { pair.second };

            out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"" << name << "\" format=\"ascii\">\n";
            for ( size_t i { 0 }; i < node_numbers_inv_.size(); ++i )
            {
                out << f->eval( node_numbers_inv_[ i ] ) << '\n';
            }
            out << "</DataArray>\n";
        }
        out << "</PointData>\n"; 
    }

    out << "<Points>\n";
    out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

    for ( auto p: node_numbers_ ) out << p.first << '\n';

    out << "</DataArray>\n";
    out << "</Points>\n";

    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
    for ( auto t: triang_ )
    {
        out << node_numbers_[ t.get_node(0) ] << ' ' 
            << node_numbers_[ t.get_node(1) ] << ' ' 
            << node_numbers_[ t.get_node(2) ] << '\n';
    }
    out << "</DataArray>\n";

    out << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
    for ( size_t i { 1 }; i <= triang_.size(); ++i ) out << 3*i << '\n';
    out << "</DataArray>\n";


    out << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
    for ( size_t i { 0 }; i < triang_.size(); ++i ) out << 5 << '\n';
    out << "</DataArray>\n";

    out << "</Cells>\n";
    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

}

