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
#include "fem/gmsh_reader.h" 

#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <boost/regex.hpp>

#include <iostream>

using namespace std;
using namespace fem;
using namespace geometry;

namespace
{

////////////////////////////////////////////////////////////
// Data structures for the contents of a gmsh MSH-file.   //
////////////////////////////////////////////////////////////

/*!
 * \brief Stores information about the format of the file.
 *
 * Instants of this structure essentially store the contents of the $MeshFormat
 * section of a MSH file. This information contains the version of the file
 * format, if the file is written in text or as binary, and if so, the byte
 * ordering.
 */
struct mesh_format
{
    double version; //!< File-format version.
    bool   binary;  //!< If the file is written in binary mode.
    bool   swap;    //!< If we need to swap the byte order in binary mode.
};

/*!
 * \brief Stores information about a single element in a MSH file.
 *
 * Instants of this data structure store all the information given about an
 * element in the MSH file. This includes the type, a finite number of integer
 * tags, and the nodes it consists of. The first tag, if existant, denotes the
 * physical group the element belongs to.
 */
struct element
{
    int type;
    vector<int> tags;
    vector<int> nodes;
}; 

/*!
 * \brief Type for storing information about the available element types.
 *
 * Instants of this structure encode the number of nodes and the name of
 * element types supported by the MSH format.
 */
struct element_type
{
    string name;
    size_t nodes;
};

/*!
 * \brief The set of sections that are parsed by the reader.
 *
 * Only the sections named in this set are parsed when reading a MSH file.
 * Every other section is ignored.
 */
const set<string> known_sections
{
    string { "MeshFormat" },
    string { "Nodes" },
    string { "Elements" },
//    string { "PhysicalNames" }
};


/*!
 * \brief Map of all element types supported by the MSH format.
 */
const map<int,element_type> type_info
{
    make_pair(  1, element_type { "  2-node line",                      2 } ),
    make_pair(  2, element_type { "  3-node triangle",                  3 } ),
    make_pair(  3, element_type { "  4-node quadrangle",                4 } ),
    make_pair(  4, element_type { "  4-node tetrahedron",               4 } ), 
    make_pair(  5, element_type { "  8-node hexahedron",                8 } ),
    make_pair(  6, element_type { "  6-node prism",                     6 } ), 
    make_pair(  7, element_type { "  5-node pyramid",                   5 } ), 
    make_pair(  8, element_type { "  3-node 2nd order line",            3 } ),
    make_pair(  9, element_type { "  6-node 2nd order triangle",        6 } ),
    make_pair( 10, element_type { "  9-node 2nd order quadrangle",      9 } ),
    make_pair( 11, element_type { " 10-node 2nd order tetrahedron",    10 } ),
    make_pair( 12, element_type { " 27-node 2nd order hexahedron",     27 } ),
    make_pair( 13, element_type { " 18-node 2nd order prism",          18 } ),
    make_pair( 14, element_type { " 14-node 2nd order pyramid",        14 } ), 
    make_pair( 15, element_type { "  1-node point",                     1 } ), 
    make_pair( 16, element_type { "  8-node 2nd order quadrangle",      8 } ), 
    make_pair( 17, element_type { " 20-node 2nd order hexahedron",     20 } ),
    make_pair( 18, element_type { " 15-node 2nd order prism",          15 } ),
    make_pair( 19, element_type { " 13-node 2nd order pyramid",        13 } ),
    make_pair( 20, element_type { "  9-node 3rd order inc. triangle",   9 } ),
    make_pair( 21, element_type { " 10-node 3rd order triangle",       10 } ),
    make_pair( 22, element_type { " 12-node 4th order inc. triangle",  12 } ),
    make_pair( 23, element_type { " 15-node 4th order triangle",       15 } ),
    make_pair( 24, element_type { " 15-node 5th order inc. triangle",  15 } ),
    make_pair( 25, element_type { " 21-node 5th order triangle",       21 } ), 
    make_pair( 26, element_type { "  4-node 3rd order edge",            4 } ),
    make_pair( 27, element_type { "  5-node 4th order edge",            5 } ),
    make_pair( 28, element_type { "  6-node 5th order edge",            6 } ),
    make_pair( 29, element_type { " 20-node 3rd order tetrahedron",    20 } ),
    make_pair( 30, element_type { " 35-node 4th order tetrahedron",    35 } ),
    make_pair( 31, element_type { " 56-node 5th order tetrahedron",    56 } ),
    make_pair( 92, element_type { " 64-node 3rd order hexahedron",     64 } ),
    make_pair( 93, element_type { "125-node 4th order hexahedron",    125 } )
};


/////////////////////////////////////
// Methods for reading MSH files.  //
/////////////////////////////////////

/*!
 * \brief  Reads the contents of an entire file into a single string.
 * \param  filename The name of the file to be read.
 * \return A string object containing the contents of the file.
 * \throw  runtime_error If an error occured while opening or reading the file.
 */
string read_file( const string filename );


/*!
 * \brief Extracts the known sections from a MSH-file.
 * \param input A string containing the contents of a MSH-file.
 *
 * This function takes the contents of a MSH-file produced by Gmsh. These
 * files are composed of several named sections. This method extracts all
 * known sections from the contents of a file and stores them as
 * (NameOfSection,ContentsOfSection)-pairs in a map. All other sections
 * or file-contents are ignored.
 */
map<string,string> get_sections( string input );


/*!
 * \brief Reads a binary int from the stream.
 * \return The int read from the stream.
 * \param  str  The stream the int should be read from.
 * \param  swap If the byte order needs to be swapped.
 */
int read_int( std::istream &str, bool swap = false );


/*!
 * \brief Reads a binary int from the stream.
 * \return The int read from the stream.
 * \param  str  The stream the int should be read from.
 * \param  swap If the byte order needs to be swapped.
 */
double read_double( std::istream &str, bool swap = false );



/*!
 * \brief Determines the mesh_format information from the file contents.
 */
mesh_format get_format( const map<string,string> &sections );

/*!
 * \brief Reads all nodes from the file contents.
 */
map<int,point> get_nodes ( const map<string,string> &sections, const mesh_format format );

/*!
 * \brief Reads all elements from the file contents.
 */
map<int,element> get_elements( const map<string,string> &sections, const mesh_format format );

}




namespace
{

map<string,string> get_sections( string input )
{
    using boost::regex;
    using boost::smatch;
    using boost::sregex_iterator;

    regex section_pattern
    {
        // Header of a section.
        R"(^\$(\w*)\s*$)"
                            
        // Contents of the section. 
        // Ignore the leading and trailing line-breaks.
        R"(\n(.*?)\n)"

        // End of a section, matching the section’s header.
        R"(^\$End\1\s*$)"
    };

    map<string,string> sections;
    for ( sregex_iterator p (input.begin(),input.end(),section_pattern);
                          p != sregex_iterator(); ++p )
    {
        smatch section = *p;
        string name    = section[1];
        string content = section[2];

        if ( known_sections.count(name) )
            sections[name] = content; 
    }
    return move(sections);
}

string read_file( string filename )
{
    string contents;
    std::ifstream str( filename, ios::in | ios::binary );

    if ( str )
    {
        str.seekg(0,ios::end); contents.resize(str.tellg());
        str.seekg(0,ios::beg); str.read(&contents[0],contents.size());
    }
    else throw runtime_error ( "Could not open file." );
    
    if ( str ) return move(contents); 
    else throw runtime_error( "Error while reading file." );
}

int read_int( std::istream &str, bool swap )
{
    union int_conv
    {
        int  as_int;
        char as_char[ sizeof(int) ];
    };

    int_conv conv;
    str.read( conv.as_char, sizeof(int) );

    if ( !str ) throw std::runtime_error( "Error while reading int in binary." );
    if ( swap ) std::reverse( conv.as_char, conv.as_char + sizeof(int) );

    return conv.as_int;
}

double read_double( std::istream &str, bool swap )
{
    union double_conv
    {
        double as_double;
        char   as_char[ sizeof(double) ];
    };

    double_conv conv;
    str.read( conv.as_char, sizeof(double) );

    if ( !str ) throw std::runtime_error( "Error while reading double in binary." );
    if ( swap ) std::reverse( conv.as_char, conv.as_char + sizeof(double) );

    return conv.as_double;
}

mesh_format get_format( const map<string,string> &sections )
{
    auto iter = sections.find( "MeshFormat" );
    if ( iter == sections.end() )
    {
        throw runtime_error( "Could not find MeshFormat section." );
    }
    stringstream str( iter->second );
  
    double version; int file_type, data_type; 
    string line; getline( str, line );
    if ( ! str ) throw runtime_error( "Could not read first line of MeshFormat section." ); 
    stringstream linestr( line );
    linestr >> version >> file_type >> data_type;

    if ( ! linestr )
        throw runtime_error ( "Error while parsing MeshFormat section." );
    if ( version != 2.2 )
        throw runtime_error( "Only MSH file version 2.2 is supported." );
    if ( file_type != 0 && file_type != 1 )
        throw runtime_error( R"(File-type must be either "0" (ascii) or "1" (binary).)" );

 
    bool swap { false }; 
    if ( file_type == 1 )
    { 
        if ( data_type != sizeof(double) )
            throw runtime_error( "data-type field does not equal sizeof(double)." );

        union int_conv
        {
            int  as_int;
            char as_char[ sizeof(int) ];
        };
        int_conv one;
        one.as_int = read_int( str, false );

        if ( one.as_int != 1 )
        {
            reverse( one.as_char, one.as_char + sizeof(int) );
            if ( one.as_int == 1 ) swap = true;
            else throw runtime_error( "Could not determine endianness of file." );
        }
    }

    return mesh_format { version, (file_type == 1), swap };
}

map<int,point> get_nodes( const map<string,string> &sections, const mesh_format format )
{
    auto iter = sections.find( "Nodes" );
    if ( iter == sections.end() )
    {
        throw runtime_error( "Could not find Nodes section." );
    }
    stringstream str( iter->second );

    string line; getline(str,line);
    if ( ! str ) throw runtime_error( "Could not read first line of Nodes section." ); 
    stringstream linestr(line);
    

    size_t N; linestr >> N;
    if ( ! linestr ) throw runtime_error( "Could not read number of nodes." );
    if ( N > iter->second.size() )
        throw runtime_error( "Number of nodes incorrect." );
   
    map<int,point> result;
    int num;
    double x, y, z;
    if ( format.binary )
    {
        for ( size_t i = 0; i < N; ++i )
        {
            num = read_int   ( str, format.swap );
              x = read_double( str, format.swap );
              y = read_double( str, format.swap );
              z = read_double( str, format.swap );
            result[ num ] = point { x, y, z };
        }
    }
    else
    {
        for ( size_t i = 0; i < N; ++i )
        {
            str >> num >> x >> y >> z;
            if ( ! str ) throw runtime_error( "Error while reading node." );
            result[ num ] = point { x, y, z };
        }
    }

    return move(result);
}

map<int,element> get_elements( const map<string,string> &sections, const mesh_format format )
{
    auto iter = sections.find( "Elements" );
    if ( iter == sections.end() )
    {
        throw runtime_error( "Could not find Elements section." );
    }
    stringstream str( iter->second );

    string line; getline(str,line);
    if ( ! str ) throw runtime_error( "Could not read first line of Elements section." ); 
    stringstream linestr(line);

    size_t N; linestr >> N;
    if ( ! linestr ) throw runtime_error( "Could not read number of elements." );
    if ( N > iter->second.size() ) throw runtime_error( "Number of elements incorrect." );

    map<int,element> result;
    if ( format.binary )
    {
        size_t i { 0 };
        while ( i < N )
        {
            size_t type  = read_int( str, format.swap );
            size_t Nsect = read_int( str, format.swap );
            size_t ntags = read_int( str, format.swap );

            if ( ntags > str.str().size() ) throw runtime_error( "Number of tags too large" );
            auto iter = type_info.find( type );
            if ( iter == type_info.end() ) throw runtime_error( "Encountered element of unknown type." );
            size_t nodes = iter->second.nodes;

            for ( size_t j = 0; j < Nsect; ++i, ++j )
            {
                element e; e.type = type;
                int number = read_int( str, format.swap );

                e.tags.reserve(ntags);
                for ( size_t k = 0; k < ntags; ++k )
                    e.tags.push_back( read_int( str, format.swap ) );

                e.nodes.reserve(nodes);
                for ( size_t k = 0; k < nodes; ++k )
                    e.nodes.push_back( read_int( str, format.swap ) );

                result[ number ] = move(e);
            }
        }
    }
    else
    {
        for ( size_t i = 0; i < N; ++i )
        {   
            element e;
            int number;
            size_t ntags;
            getline(str,line);
            if ( ! str ) throw runtime_error( "Could not read line in Elements section." ); 
            linestr.str(line); linestr.seekg(0,ios::beg);

            linestr >> number >> e.type >> ntags;
            if ( ! linestr ) throw runtime_error( "Error while reading Elements section." );
            if ( ntags > line.size() ) throw runtime_error( "Number of tags too large" );

            auto iter = type_info.find( e.type );
            if ( iter == type_info.end() ) throw runtime_error( "Encountered element of unknown type." );
            size_t nodes = iter->second.nodes;

            e.tags.reserve(ntags);
            for ( size_t j = 0; j < ntags; ++j )
            {
                int tag;
                linestr >> tag;
                if ( ! linestr ) throw runtime_error( "Could not read all element tags." );
                e.tags.push_back(tag);
            }

            e.nodes.reserve(nodes);
            for ( size_t j = 0; j < nodes; ++j )
            {
                int node;
                linestr >> node;
                if ( ! linestr ) throw runtime_error( "Could not read all element nodes." );
                e.nodes.push_back(node);
            }
            result[ number ] = move(e);
        }
    }
    return move(result);
}

}


namespace fem
{

multigrid<1> gmsh_reader::make_grid()
{
    using std::map;
    using std::pair;
    using std::string;
    using geometry::point;

    map<string,string> contents { get_sections( read_file(filename) ) };
    mesh_format        format   { get_format  (contents) };        contents.erase( "MeshFormat" );
    map<int,point>     nodes    { get_nodes   (contents,format) }; contents.erase( "Nodes" );
    map<int,element>   elements { get_elements(contents,format) }; contents.erase( "Elements" );
    contents.clear();

    multigrid<1> mg;
    for ( pair<const int,element> &p: elements )
    {
        element &e { p.second };
        switch ( e.type )
        {
        case  1: break; // 2-node line.
        case  2: break; // 3-node triangle.
        case  4:        // 4-node tetrahedron.
        {
            shapefcts3d::coeffs<1,point> t;
            std::sort( e.nodes.begin(), e.nodes.end() );
            for ( size_t i = 0; i < 4; ++i )
            {
                auto iter = nodes.find( e.nodes[i] );
                if ( iter == nodes.end() )
                {
                    throw runtime_error( "Element referring to non-existant node." );
                }
                t[i] = iter->second;
            }
            add( mg, t ); 
        } break;
        case 15: break; // 1-node point.
        default: throw runtime_error( "Element type not supported yet." );
        }
    }

    return std::move(mg);
}

}

