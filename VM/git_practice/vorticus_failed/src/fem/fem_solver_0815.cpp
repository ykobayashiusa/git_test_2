#include "fem/multigrid.h"
#include "fem/vtkwriter.h"
#include "fem/gmsh_reader.h"
#include "libvorticus.h" 

#include <fstream>
#include <sstream>
#include <iostream>

using namespace fem;

int main(int argc, char *argv[] )
{
	
	multigrid<1> mg = read_gmsh( argv[1] );
	std::cout<<"done reading mesh." <<std::endl;

	for (size_t lvl = 0; lvl < 1; ++lvl)
	{
		//size_t lvl=0;
	
		// marking mesh, reqired for mesh refinement 
		for ( auto it = mg.grid_begin(lvl); it !=mg.grid_end(lvl); ++it )
		it->set_ref();
	
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

	std::cout<<"numb.nodes() " << numb.nodes()[9].x  <<std::endl;
	std::cout<<"numb.cell_nodes() : " <<numb.cell_nodes().begin()->second[0]<<std::endl;
	std::cout << "numb_test : " << numb.nodes()[
		std::next( numb.cell_nodes().begin(), std::distance( numb.cell_nodes().begin(), numb.cell_nodes().end() ) - 1 )->second[0]
		].x << std::endl;
	std::cout<<"numb_size : "<<numb.cell_nodes().begin()->second.size()<<std::endl;





	int ndof=1;
	//std::cout<<"ndof ? : "<<unity::shapefcts::dofs_per_node<1>()<<std::endl;

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


	// loop for the total number of elements
	for (auto nums = numb.cell_nodes().begin(); nums != numb.cell_nodes().end(); ++nums )
	//for (auto nums = numb.cell_nodes().begin(); nums != std::next( numb.cell_nodes().begin(), 1 ); ++nums )
	{
		auto xbar = math::matrix( nums->second.size(), nums->second.size(), arma::fill::zeros );
		for ( int i=0; i<nums->second.size(); ++i )
		{
			xbar(i,0) = 1;
			xbar(i,1) = numb.nodes()[ nums->second[i] ].x;   // nums->second[i] : returns global node number of the i-th node consisting an element (tetrahedron)
			xbar(i,2) = numb.nodes()[ nums->second[i] ].y;
			xbar(i,3) = numb.nodes()[ nums->second[i] ].z;
		}	
		
		math::matrix xinv { xbar.i() };   //		xinv.print("xinv : ");
		double vol = 0.166666666666666667* det(xbar); // 1/6*det(xbar)
		//double vol = 0.1667* det(xbar); // 1/6*det(xbar)
		std::cout<<"vol : "<<vol<<std::endl;
		//xbar.print("xbar :");

		// compute element matrix
		auto k = math::matrix( nums->second.size(), nums->second.size(), arma::fill::zeros );
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


		// compute system dofs associated with each element
		int k_ = 0;   math::vector index( ndof*nums->second.size() ); 
		for ( int i=0; i<nums->second.size(); ++i )
		{
			int start = ( nums->second[i] )*ndof;
			for ( int j=0; j<ndof; ++j )
			{
				index(k_) = start + j;	k_ = k_ + 1; 
			}
		}
		//index.print("index : ");
			
		// assemble into system matrix
		for ( int i=0; i<4; ++i )
		{
			int ii = index[i];
			for (int j=0; j<4; ++j )
			{
				int jj = index[j];
				kk(ii,jj) = kk(ii,jj) + k(i,j);
			}
		}
		std::cout << "iteration : " << std::distance( numb.cell_nodes().begin(), nums )+1 << std::endl;
	}
	kk.print("system matrix kk :");


	// creating grid_function for u
	grid_function<real,1,1>	u( numb );
	std::cout<<"grid_function before : "<<u(numb(13))<<std::endl;
	
	//solving system matrix equation
	math::vector x ( ndof*numb.size(), arma::fill::zeros );
	math::cg( kk, ff, x, 1e-12 );
	std::cout<<"x : "<<x<<std::endl;

	//math::vector fsol { solve(kk,ff) };  std::cout<<fsol<<std::endl;

	std::cout<<"grid_function after : "<<u(numb(13))<<std::endl;
	std::cout<<u.eval(numb(13))<<std::endl;

	


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



	}
}

