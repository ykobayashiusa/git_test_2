#include "fem/multigrid.h"
#include "fem/vtkwriter.h"
#include "fem/gmsh_reader.h"
#include "libvorticus.h" 

#include <fstream>
#include <sstream>
#include <iostream>

using namespace fem;


//math::matrix make_elem_matrix () 

void test_function ( multigrid<1>::grid_iterator &it )
{
	it->set_ref();
}

void test_function_2 ( multigrid<1>::grid_iterator &it ) 
{
	std::cout<<"meow"<<std::endl;
}



int main(int argc, char *argv[] )
{

	multigrid<1> mg = read_gmsh( argv[1] );
	std::cout<<"done reading mesh." <<std::endl;

	for (size_t lvl = 0; lvl < 1; ++lvl)
	{
		//size_t lvl=0;
	
		// marking mesh, reqired for mesh refinement 
		for ( auto it = mg.grid_begin(lvl); it !=mg.grid_end(lvl); ++it )
		test_function(it);
		//	it->set_ref();
	
		// mesh refinement based on multilevel refinement algorithm
		mg.adapt();


		node_numbering<1,1>	numb( mg.grid_begin(lvl), mg.grid_end(lvl) );
		grid_function<real,1,1>	numb_fct( numb );

		for ( size_t i=0; i<numb.size(); ++i )
			numb_fct( numb(i) ) = i;


	// writing out .vtu file for mesh visualization using paraview
	vtkwriter<1> writer( mg.grid_begin(lvl), mg.grid_end(lvl) );
	writer.register_scalar( numb_fct, "number" );

	std::stringstream filename; filename << "lol_" << lvl <<".vtu";
	std::ofstream file( filename.str() );
	writer.write( file );

	std::cout<<filename.str() <<std::endl;

	

	std::cout<< "number of points ( total number of nodes ) : " << numb.size() << std::endl;
	std::cout<< "number of cells ( elements(tetrahedra) ) :  "<< numb.cell_nodes().size() <<std::endl;
	

	// node_numbering output for checking in Matlab
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







	int nnel = numb.cell_nodes().begin()->second.size();  // number of nodes per element
	int ndof=1; // degree of freedom per node ( always have one degree for conventional FEM )

	math::matrix kk ( ndof*numb.size(), ndof*numb.size(), arma::fill::zeros );
		std::cout << "matrix : " << kk.n_rows  << std::endl;
//	math::hash_matrix kk_ ( ndof*numb.size(), ndof*numb.size() ); 
	//auto kk_ = math::hash_matrix ( ndof*numb.size(), ndof*numb.size() ); 
	//std::cout<<"hash_matix : "<<kk_<<std::endl;


	// creating rhs vector and grid_function 
	grid_function<real,1,1>	b( numb );
	math::vector ff ( ndof*numb.size(), arma::fill::zeros );
	for (auto nums = 0; nums<numb.size(); ++nums) 
	{
		ff(nums) = 2* numb.nodes()[ nums ].y; // vorticity for Channel flow 
		b( numb(nums) ) = ff(nums);
	}
	for (int i=0; i<numb.size(); ++i) std::cout <<
		"rhs_vector ff : "  	<< ff(i) << ", "
		"grid_function b : " 	<< b(numb(i)) << ", " 
		"grid_function eval : "	<< b.eval(numb(i)) << std::endl; 


	// temporal input data for boundary conditions valid only for cube.msh
	std::vector<int> bcdof;
	for ( int i=0; i<numb.size(); ++i )
		if (    numb.nodes()[i].x == 0.5 || numb.nodes()[i].x == -0.5 ||
			numb.nodes()[i].y == 0.5 || numb.nodes()[i].y == -0.5 ||	
			numb.nodes()[i].z == 0.5 || numb.nodes()[i].z == -0.5 )
			bcdof.push_back(i);
	for ( int i=0; i<bcdof.size(); ++i ) std::cout<<"bcdof : "<<bcdof[i]<<std::endl;
	
	/*
	std::unordered_map<int,std::vector<double>> testing_;
	testing_.insert({1,{}});
	testing_.at(1).push_back(4);
	//testing_.begin()->second.push_back(3); 
	std::cout<<"testing_ : "<<testing_.size()<<std::endl;
	std::cout<<"testing_ : "<<testing_.begin()->second.size()<<std::endl;
	std::cout<<"testing_ : "<<testing_.begin()->second[0]<<std::endl;
*/

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
	
	
	
	for (int j=0; j<bcdof_.size(); ++j ) // for 6 faces
	for (int i=0; i<bcdof_.at(j).size(); ++i ) // for every node in j-face
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


	for ( auto it = mg.grid_begin(lvl); it !=mg.grid_end(lvl); ++it )
	{
		std::cout << "get_node : " << it->get_node(1).x << std::endl; // replace with numb.nodes()[ nums->second[1] ].x
		std::cout << "get_node_num : " << std::next( numb.cell_nodes().begin(), std::distance( mg.grid_begin(lvl), it ) )->second[1]<<std::endl; // replace with nums->second[1]
	}
	// loop for the total number of elements
	for ( auto nums = numb.cell_nodes().begin(); nums != numb.cell_nodes().end(); ++nums )
	//for (auto nums = numb.cell_nodes().begin(); nums != std::next( numb.cell_nodes().begin(), 1 ); ++nums )
	{
		
		//test_function_2(it);
		math::matrix xbar ( nnel, nnel, arma::fill::zeros );
		for ( int i=0; i < nnel; ++i )
		{
			xbar(i,0) = 1;
			xbar(i,1) = numb.nodes()[ nums->second[i] ].x;   // nums->second[i] : returns global node number of the i-th node consisting an element (tetrahedron)
			xbar(i,2) = numb.nodes()[ nums->second[i] ].y;
			xbar(i,3) = numb.nodes()[ nums->second[i] ].z;
		}	
		
		math::matrix xinv { xbar.i() };   //		xinv.print("xinv : ");
		//double vol = 0.166666666666666667* det(xbar); // 1/6*det(xbar)
		double vol = det(xbar)/6;  
		//double vol = 0.1667* det(xbar); // 1/6*det(xbar)
		std::cout<<"vol : "<<vol<<std::endl;
		//xbar.print("xbar :");

		// compute element matrix
		math::matrix k ( nnel, nnel, arma::fill::zeros );
			k(0,0) = xinv(1,0)*xinv(1,0) + xinv(2,0)*xinv(2,0) + xinv(3,0)*xinv(3,0); 
			k(0,1) = xinv(1,0)*xinv(1,1) + xinv(2,0)*xinv(2,1) + xinv(3,0)*xinv(3,1); 
			k(0,2) = xinv(1,0)*xinv(1,2) + xinv(2,0)*xinv(2,2) + xinv(3,0)*xinv(3,2); 
			k(0,3) = xinv(1,0)*xinv(1,3) + xinv(2,0)*xinv(2,3) + xinv(3,0)*xinv(3,3); 
			k(1,0) = k(0,1); 
			k(1,1) = xinv(1,1)*xinv(1,1) + xinv(2,1)*xinv(2,1) + xinv(3,1)*xinv(3,1); 
			k(1,2) = xinv(1,1)*xinv(1,2) + xinv(2,1)*xinv(2,2) + xinv(3,1)*xinv(3,2); 
			k(1,3) = xinv(1,1)*xinv(1,3) + xinv(2,1)*xinv(2,3) + xinv(3,1)*xinv(3,3); 
			k(2,0) = k(0,2); 
			k(2,1) = k(1,2); 
			k(2,2) = xinv(1,2)*xinv(1,2) + xinv(2,2)*xinv(2,2) + xinv(3,2)*xinv(3,2); 
			k(2,3) = xinv(1,2)*xinv(1,3) + xinv(2,2)*xinv(2,3) + xinv(3,2)*xinv(3,3); 
			k(3,0) = k(0,3); 
			k(3,1) = k(1,3); 
			k(3,2) = k(2,3); 
			k(3,3) = xinv(1,3)*xinv(1,3) + xinv(2,3)*xinv(2,3) + xinv(3,3)*xinv(3,3); 
			k = vol*k;	//	k.print("k : ");


		// compute element vector 
		math::vector f ( nnel, arma::fill::zeros );
			f(0) = 0;
			f(1) = 0;
			f(2) = 0;
			f(3) = 0;

		// compute system dofs associated with each element
		math::vector index ( ndof*nnel );  int k_ = 0;   
		for ( int i=0; i < nnel; ++i )
		{
			int start = ( nums->second[i] )*ndof;
			for ( int j=0; j < ndof; ++j )
			{
				index(k_) = start + j;	k_ = k_ + 1; 
			}
		}
		//index.print("index : ");
			

		

		// assemble into system matrix and vector
		for ( int i = 0; i < ndof*nnel; ++i ) // loop for element rows
		{
			int ii = index[i]; // address for the system rows
			//ff(ii) = ff(ii) + f(i); // assembly into the system vector 
			for (int j = 0; j < index.size(); ++j ) // loop for element columns
			{
				int jj = index[j]; // address for the system column
				kk(ii,jj) = kk(ii,jj) + k(i,j); // assembly into the system matrix
			}
		}
		std::cout << "iteration : " << std::distance( numb.cell_nodes().begin(), nums )+1 << std::endl;
	}
	kk.print("system matrix kk :");


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
	//kk.print("system matrix kk with B.C. :");
	
	// creating grid_function for u
	grid_function<real,1,1>	u( numb );
	std::cout<<"grid_function before : "<<u(numb(13))<<std::endl;
	
	//solving system matrix equation
	math::vector	x ( ndof*numb.size(), arma::fill::zeros );
	//math::diag_precond        Q { M.diagonal() };
	math::cg( kk, ff, x, 1e-12 );
	std::cout<<"x : "<<x<<std::endl;

	//math::vector fsol { solve(kk,ff) };  std::cout<<fsol<<std::endl;

	std::cout<<"grid_function after : "<<u(numb(13))<<std::endl;
	std::cout<<u.eval(numb(13))<<std::endl;

	
	// playing around with grid_function
	for (int i = 0; i < x.size(); ++i )
	{
		std::cout<<
			"x : "<<x(i)<<", "<<
			"grid_function u : "<<u(numb(i))<<", "<<
			"u.eval : "<<u.eval(numb(i))<<", ";
		u(numb(i)) = x(i);
		std::cout<<
			"u_after : "<< u(numb(i))<<", "<<
			"u_after.eval : "<<u.eval(numb(i))<<std::endl;
	}
	grid_function<real,1,1> u_( numb );
	for (int i = 0; i < x.size(); ++i )
		std::cout<<
			"u : "<< u(numb(i))<<", "<<
			"u_ : "<< u(numb(i))<<", "<<std::endl;





	math::matrix test { 
		{1,2,3},
		{4,5,6}     	}; std::cout<<test<<std::endl;

	// playing around with multigrid
	for ( auto it = mg.grid_begin(lvl); it != std::next( mg.grid_begin(lvl), 3 )  ; ++it )
	for (int num=0; num<4; ++num)
	std::cout<< "get_node(" << num << ") of element " << std::distance( mg.grid_begin(lvl), it ) << " : " << it->get_node(num) <<std::endl;
	
	std::cout<< "get_node(0) : " << std::next( mg.grid_begin(lvl), std::distance( mg.grid_begin(lvl), mg.grid_end(lvl) )-1 )->get_node(0) <<std::endl; 

	//std::cout<<"mg_id() : "<<mg.grid_begin(lvl)->id() <<std::endl;
	std::cout<<"mg_centre() : "<<mg.grid_begin(lvl)->centre()  <<std::endl;

	
	std::cout<<"mg.grid_begin : "<<mg.grid_begin(lvl)->get_node(4)<<std::endl;
	geoid id = mg.grid_begin(lvl)->id();
	std::cout<<"id.dim : "<<id.dim<<", id.lvl : "<<id.lvl<<", id.pos : "<<id.pos<<std::endl;
	//astd::cout<<"mg.grid_begin : "<<mg.grid_begin(lvl)<<std::endl;


	vtkwriter<1> writer_( mg.grid_begin(lvl), mg.grid_end(lvl) );
	writer_.register_scalar( u_, "u value" );

	std::stringstream filename_; filename_ << "sol_" << lvl <<".vtu";
	std::ofstream file_( filename_.str() );
	writer_.write( file_ );

	std::cout<<filename_.str() <<std::endl;
	}
}

