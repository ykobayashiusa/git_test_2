#include <fstream>
#include <sstream>
#include <iostream>
#include "libvorticus.h"

using namespace fem;

void temporal_BC ( node_numbering<1,1> &numb, math::hash_matrix &kk, math::vector &ff )
{
	// temporal input data for boundary conditions valid only for cube.msh
	std::vector<int> bcdof;
	for ( int i=0; i<numb.size(); ++i )
		if (    numb.nodes()[i].x == 0.5 || numb.nodes()[i].x == -0.5 ||
			numb.nodes()[i].y == 0.5 || numb.nodes()[i].y == -0.5 ||	
			numb.nodes()[i].z == 0.5 || numb.nodes()[i].z == -0.5 )
			bcdof.push_back(i);
	for ( int i=0; i<bcdof.size(); ++i ) std::cout<<"bcdof : "<<bcdof[i]<<std::endl;
	

	std::unordered_map<int,std::vector<int>> bcdof_;
	for ( int i=0; i<6; ++i ) bcdof_.insert({i,{}});
	for ( int i=0; i<bcdof_.size(); ++i )
		std::cout<<"bcdof__first : "<<std::next( bcdof_.begin(), i )->first<<std::endl;
	
	// sorting nodes at the boundary into 6 faces
	for ( int i=0; i<numb.size(); ++i )
		if ( numb.nodes()[i].x ==  0.5 ) bcdof_.at(0).push_back(i); // front 
	for ( int i=0; i<numb.size(); ++i )
		if ( numb.nodes()[i].x == -0.5 ) bcdof_.at(1).push_back(i); // back
	for ( int i=0; i<numb.size(); ++i )
		if ( numb.nodes()[i].y ==  0.5 ) bcdof_.at(2).push_back(i); // right
	for ( int i=0; i<numb.size(); ++i )
		if ( numb.nodes()[i].y == -0.5 ) bcdof_.at(3).push_back(i); // left
	for ( int i=0; i<numb.size(); ++i )
		if ( numb.nodes()[i].z ==  0.5 ) bcdof_.at(4).push_back(i); // top
	for ( int i=0; i<numb.size(); ++i )
		if ( numb.nodes()[i].z == -0.5 ) bcdof_.at(5).push_back(i); // bottom
	
	for ( int j=0; j<bcdof_.size(); ++j )	
	for ( int i=0; i<bcdof_.at(0).size(); ++i )
		std::cout<<"dcdof.at("<<j<<") : "<<bcdof_.at(j)[i]<<std::endl;

	double B_C_=0; // temporal coefficient
	std::cout<<"bcdof__size : "<<bcdof_.size()<<std::endl;
	std::cout<<"bcdof__second_size : "<<bcdof_.at(0).size()<<std::endl;
	std::vector<double> bcval( bcdof.size() );
	std::cout<<"bcval.size() : "<<bcval.size()<<std::endl;
	
	
	for (int j=0; j<bcdof_.size(); ++j ) // loop for 6 faces
	for (int i=0; i<bcdof_.at(j).size(); ++i ) // loop for every node in j-face
	{
		int ii=0;
		while ( bcdof_.at(j)[i] != bcdof[ii] ) ii=ii+1; // selecting nodes of j-face from bcdof
			std::cout<<"ii : "<<ii<<", bcdof(ii) : "<<bcdof[ii]<<", bcdof_.at("<<j<<")[i]: "<<bcdof_.at(j)[i]<<std::endl;
		if ( j==0 )  bcval[ii] = numb.nodes()[ bcdof[ii] ].x*B_C_; // B.C. value at front face
		if ( j==1 )  bcval[ii] = numb.nodes()[ bcdof[ii] ].y*B_C_; // B.C. value at back face
		if ( j==2 )  bcval[ii] = numb.nodes()[ bcdof[ii] ].z*B_C_; // B.C. value at right face
		if ( j==3 )  bcval[ii] = numb.nodes()[ bcdof[ii] ].x*B_C_; // B.C. value at left face
		if ( j==4 )  bcval[ii] = numb.nodes()[ bcdof[ii] ].y*B_C_; // B.C. value at top face
		if ( j==5 )  bcval[ii] = numb.nodes()[ bcdof[ii] ].z*B_C_; // B.C. value at bottom face
	}
	for (int i=0; i<bcval.size(); ++i )
	std::cout<<"bcval : "<<bcval[i]<<std::endl;	
	

	// apply boundary condition without destroying the symmetry of the system matrix
	for ( int ic = 0; ic < bcdof.size(); ++ic ) // loop for all constraints
	{
		for ( int i = 1; i < ff.size(); ++i ) // loop for number of equations in system
		{
			ff(i) = ff(i) - bcval[ic]*kk( i, bcdof[ic] ); // modify column using constrained value
			kk(bcdof[ic],i) = 0; // set all the bcdof[ic]-th row to zero
			kk(i,bcdof[ic]) = 0; // set all the bcdof[ic]-th column to zero
		}
		kk(bcdof[ic],bcdof[ic]) = 1; // set the bcdof[ic]-th diagonal to unity
		ff(bcdof[ic])		= bcval[ic]; // put the constrained value in the column
	}
}



void mesh_output ( node_numbering<1,1> &numb )	// node_numbering output for checking in Matlab
{
	for (auto nums = numb.cell_nodes().begin(); nums != numb.cell_nodes().end(); ++nums ) {
		for ( int i=0; i<nums->second.size(); ++i ) std::cout <<
		       "nodes(" << std::distance( numb.cell_nodes().begin(), nums )+1 << ","<<i+1<<")=" << nums->second[i]+1<< "; ";
	std::cout<<std::endl;}

	for (auto nums = 0; nums<numb.size(); ++nums) std::cout << 
		       "gcoord("<< nums+1 <<",1)=" << numb.nodes()[ nums ].x << "; "<<
		       "gcoord("<< nums+1 <<",2)=" << numb.nodes()[ nums ].y << "; "<<
		       "gcoord("<< nums+1 <<",3)=" << numb.nodes()[ nums ].z << "; "<<std::endl;

	for (auto nums = numb.cell_nodes().begin(); nums != numb.cell_nodes().end(); ++nums )
		for ( int i=0; i<nums->second.size(); ++i ) std::cout << 
		       "gcoord("<< nums->second[i]+1 <<",1)=" << numb.nodes()[ nums->second[i] ].x << "; "<<
		       "gcoord("<< nums->second[i]+1 <<",2)=" << numb.nodes()[ nums->second[i] ].y << "; "<<
		       "gcoord("<< nums->second[i]+1 <<",3)=" << numb.nodes()[ nums->second[i] ].z << "; "<<std::endl;

	std::cout<<"norm of numb.nodes() " << numb.nodes()[9].r()  <<std::endl;
	std::cout<<"numb.cell_nodes() : " <<numb.cell_nodes().begin()->second[0]<<std::endl;
	std::cout << "numb_test : " << numb.nodes()[
		std::next( numb.cell_nodes().begin(), std::distance( numb.cell_nodes().begin(), numb.cell_nodes().end() ) - 1 )->second[0]
		].x << std::endl;
	std::cout<<"numb_size : "<<numb.cell_nodes().begin()->second.size()<<std::endl;
}



void notes ()
{
	math::matrix test_matrix { 
		{1,2,3},
		{4,5,6}     	}; std::cout<<test_matrix<<std::endl;
}


int test_133 ( math::matrix kk, int &nnel )
{
	int test_133_=133;
	return test_133_;
}


