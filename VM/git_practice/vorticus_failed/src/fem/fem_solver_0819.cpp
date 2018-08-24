#include "fem/multigrid.h"
#include "fem/vtkwriter.h"
#include "fem/gmsh_reader.h"
#include "libvorticus.h" 
#include "geometry/tetrahedral_quadrules.cpp"
#include "fem/Y_temporal_BC.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/iterator/filter_iterator.hpp>


using namespace fem;

math::matrix make_elem_mat ( multigrid<1>::grid_iterator &it, int &nnel ) // compute element matrix
{
	math::matrix xbar ( nnel, nnel, arma::fill::zeros );
	for ( int i=0; i < nnel; ++i )
	{
		xbar(i,0) = 1;
		xbar(i,1) = it->get_node(i).x;   // it->get_node(i) : returns global node number of the i-th node consisting an element (tetrahedron)
		xbar(i,2) = it->get_node(i).y;
		xbar(i,3) = it->get_node(i).z;
	}	//xbar.print("xbar :");	
	
	math::matrix xinv { xbar.i() };   // xinv.print("xinv : ");
	double vol = det(xbar)/6; std::cout<<"vol : "<<vol<<std::endl;

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
	return k;
}

math::matrix make_elem_mat_iso ( multigrid<1>::grid_iterator &it, int &nnel ) // compute element matrix for isoparametric element
{
	
}

math::matrix make_elem_rhs ( multigrid<1>::grid_iterator &it, int &nnel ) // compute element vector
{
	math::vector f ( nnel, arma::fill::zeros );	
		f(0) = 0;
		f(1) = 0;
		f(2) = 0;
		f(3) = 0;
	return f;
}

math::vector make_index ( multigrid<1>::grid_iterator &it, node_numbering<1,1> &numb, int &ndof, int &nnel )
{
	math::vector index ( ndof*nnel ); int k_ = 0;   
	for ( int i=0; i < nnel; ++i )
	{
		//int start = ( std::next( numb.cell_nodes().begin(), std::distance( mg.grid_begin(lvl), it ) )->second[i] )*ndof;
		int start =  numb.cell_nodes().at( it->id() )[i] * ndof;
		for ( int j=0; j < ndof; ++j )
		{
			index(k_) = start + j;	k_ = k_ + 1; 
		}
	}	//index.print("index : ");
	return index;
}

math::matrix assemble_kk ( math::matrix &kk, math::matrix &k, math::vector &ff, math::vector &f, math::vector &index, int &ndof, int &nnel )
{
	// assemble into system matrix and vector
	for ( int i = 0; i < ndof*nnel; ++i ) 		// loop for element rows
	{
		int ii = index[i]; 			// address for the system rows
		//ff(ii) = ff(ii) + f(i); 		// assembly into the system vector 
		for (int j = 0; j < index.size(); ++j )	// loop for element columns
		{
			int jj = index[j]; 		// address for the system column
			kk(ii,jj) = kk(ii,jj) + k(i,j);	// assembly into the system matrix
		}
	}
	return kk;
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
		it->set_ref();
	
		mg.adapt(); // mesh refinement based on multilevel refinement algorithm

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
	mesh_output ( numb ); // node_numbering output for checking in Matlab

	std::cout<< "number of points ( total number of nodes ) : " << numb.size() << std::endl;
	std::cout<< "number of cells ( elements(tetrahedra) ) :  "<< numb.cell_nodes().size() <<std::endl;
	


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


	// creating system matrix and vector
	for ( auto it = mg.grid_begin(lvl); it !=mg.grid_end(lvl); ++it ) // loop for the total number of elements
	{
		math::matrix k { make_elem_mat(it, nnel) };	// compute element matrix
		//math::matrix k { make_elem_mat_iso (it, nnel) };// compute element matrix for isoparametric element
		math::vector f { make_elem_rhs(it, nnel) };	// compute element vector 
		math::vector index { make_index( it, numb, ndof, nnel ) }; 	// compute system dofs associated with each element
		kk = assemble_kk ( kk, k, ff, f, index, ndof, nnel );	// assemble into system matrix and vector
		
		std::cout << "iteration : " << std::distance( mg.grid_begin(lvl), it )+1 << std::endl;
	}
	kk.print("system matrix kk :");

	//temporal_BC ( numb, kk, ff ); // apply boundary condition without destroying the symmetry of the system matrix
	kk.print("system matrix kk with B.C. :");

	// creating grid_function for u
	grid_function<real,1,1>	u( numb );
	std::cout<<"grid_function before : "<<u(numb(13))<<std::endl;
	
	//solving system matrix equation
	math::vector	x ( ndof*numb.size(), arma::fill::zeros );
	std::cout<<"grid_function before : "<<u(numb(13))<<std::endl;

	//math::diag_precond        Q { M.diagonal() };
//	math::cg( kk, ff, x, 1e-12 );
//	math::cg( kk, ff, x, 1e-6 );
	std::cout<<"x : "<<x<<std::endl;
	//math::vector fsol { solve(kk,ff) };  std::cout<<fsol<<std::endl;

	std::cout<<"grid_function after : "<<u(numb(13))<<std::endl;
	std::cout<<u.eval(numb(13))<<std::endl;

	
	// studying grid_function
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
			"u_ : "<< u_(numb(i))<<", "<<std::endl;


	vtkwriter<1> writer_( mg.grid_begin(lvl), mg.grid_end(lvl) );
	writer_.register_scalar( u, "u value" );

	std::stringstream filename_; filename_ << "sol_" << lvl <<".vtu";
	std::ofstream file_( filename_.str() );
	writer_.write( file_ );

	std::cout<<filename_.str() <<std::endl;


// end of FEM code.


//////////////////////////////////
// notes & tests for code study //
//////////////////////////////////

	// studying multigrid
	for ( auto it = mg.grid_begin(lvl); it != std::next( mg.grid_begin(lvl), 3 ) ; ++it )
	for (int num=0; num<4; ++num)
	std::cout<< "get_node(" << num << ") of element " << std::distance( mg.grid_begin(lvl), it ) << " : " << it->get_node(num) <<std::endl;
	
	std::cout<< "get_node(0) : " << std::next( mg.grid_begin(lvl), std::distance( mg.grid_begin(lvl), mg.grid_end(lvl) )-1 )->get_node(0) <<std::endl; 

	std::cout<<"mg_id().pos : "<<mg.grid_begin(lvl)->id().pos << std::endl;
	std::cout<<"mg_centre() : "<<mg.grid_begin(lvl)->centre() << std::endl;
	

	std::array<bary3d,4> bary_ {bary3d(0,0,0),bary3d(1,0,0),bary3d(0,1,0),bary3d(0,0,1)};
	// checking 
	//for ( auto it = mg.grid_begin(lvl); it !=mg.grid_end(lvl); ++it )
	for ( auto it = mg.grid_begin(lvl); it != std::next( mg.grid_begin(lvl), 4 ) ; ++it )
	{

		for ( int i=0; i<nnel; ++i )
		{
			//using shapefcts3d::dNdxi;
			//std::array<real,4> dNdxi(bary_[i]);
			
			//std::array<real,4> shapefcts3d::dNdxi(bary_[i]);

			bary3d bary_11(1,0,0);
			using shapefcts3d::vals;
			using shapefcts3d::dNdxi;
			vals<1> dxi={1,2,3,4};
			//vals<1> dxi=shapefcts3d::dNdxi(bary_[i]);
			//vals<1> dNdxi;
			//std::array<real,4>=dNdxi(bary3d p);
			//std::array<real,4>=dNdxi(bary3d p(1,0,0));
			//dNdxi(bary_11) dxi_;
			dNdxi<1>(bary_11);
			dNdxi<1>(bary3d(1,0,0));
			//size_t shapefcts3d::num<1>();

			//using  vals = std::array<real,num<1>()>;
			//using  vals = std::array<real,4>;
			//vals shapefcts3d::dNdxi (bary_[i]);  
			//astd::array<real,4> dxi { shapefcts3d::dNdxi (bary_[i]) };
			/*
			std::array<point,4> nodes__ 
			{ 
				{it->get_node(0) },
				{it->get_node(1) },
				{it->get_node(2) },
				{it->get_node(3) }
			};*/
			//shapefcts3d::coeffs<1,point> nodes_{};
			//point dxi = shapefcts3d::apply<1,point>( nodes__, shapefcts3d::dNdxi  <1>(bary_[i]) );
		//tetrahedron<1> &t { *it };  // same as { mg.get_cell(it->id()) };
		
		std::cout<<std::fixed<<
			"elem " << std::distance( mg.grid_begin(lvl), it )+1 <<	", pos : " << it->pos() <<
			", node " << numb(*it)[i] <<" : "<< numb(numb(*it)[i]) <<
			", Chi () : " << it->Chi(bary_[i]) << ", u.eval(pos) : " << u.eval(numb(numb(*it)[i])) << 
			", u.eval(it,bary) : " << u.eval(*it,bary_[i]) <<"\n" // std::endl;
			", dChi () : " << it->dChi(bary_[i]) << 			// Jacobian
			", dChi.det() : "<< it->dChi(bary_[i]).det() << 		// determinant of Jacobian
			", dChi.inv()(0,0) : "<< it->dChi(bary_[i]).inv()(0,0) << //std::endl;//<< 		// inverse of Jacobian matrix
			", dNdxi : " << dNdxi<1>(bary_[i])[1] <<//std::endl; //dNdxi<1>(bary_[i]) << std::endl;//apply<1,point>( nodes_, shapefcts3d::dNdxi<1>(bary_[i]) )a
			", dhdx : " << it->dChi(bary_[i]).inv()(0,0)*dNdxi<1>(bary_[i])[0]<<std::endl;

		//std::cout<<
		//		std::endl;
		
		
		
		}
		/*
		std::cout << 
			"it->get_node : " << it->get_node(i) << ", " <<
			"numb.nodes() : " << numb.nodes()[ std::next( numb.cell_nodes().begin(), std::distance( mg.grid_begin(lvl), it ) )->second[i] ] << std::endl; // replaced with numb.nodes()[ nums->second[i] ].x
		std::cout << "get_node_num : " << std::next( numb.cell_nodes().begin(), std::distance( mg.grid_begin(lvl), it ) )->second[i] << std::endl; // replaced with nums->second[i]
		std::cout << "numb.cell_nodes().at(id) : "<< numb.cell_nodes().at( it->id() )[i] << std::endl;
		*/
	}
	
	//size_t order=1;
	geometry::tetrahedral_quadrule rule { get_refined_tetrahedral_quadrule(1,lvl) };
	auto begin = quadrule_iterator<1>::begin( mg.grid_begin(0), mg.grid_end(0), rule.begin(), rule.end() );
	auto   end = quadrule_iterator<1>::end  ( mg.grid_begin(0), mg.grid_end(0), rule.begin(), rule.end() );

	for ( auto it = begin; it != end; ++it )
	{
		std::cout<<//numb((*it))[1]<<", "<<
			"pos() : "<<std::next( mg.grid_begin(lvl), std::distance(begin, it) )->pos()<<
			", it.x : "<<(*it).x<<
			", it.w : "<<(*it).w<<			
			", rule.size() : "<<rule.size()<<
			", rule w : "<<rule.begin()->w<<
			//", rule b : "<<rule.begin()<<
			//
			std::endl;
	}


	for ( auto it = begin; it != end; ++it )
	{
		//it->t;
		int i = std::distance(begin,it);
		//math::vector dhdr { dNdxi<1>(bary_[i]) }; 
		math::vector dhdx { std::distance(begin,end), arma::fill::zeros };
		
		for ( int j = 0; j < nnel; ++j )
		{
		//	dhdx = it->dChi(bary_[i]).inv()(0,0)*dNdxi<1>(bary_[i])[0]
		}
		

		//for dNdxi<1>(bary[i])[0]
		//dhdx =
		std::cout<<nnel<<std::endl;
		std::cout<<//numb((*it))[1]<<", "<<
			"pos() : "<<std::next( mg.grid_begin(lvl), std::distance(begin, it) )->pos()<<
			", it.x : "<<(*it).x<<
			", it.w : "<<(*it).w<<			std::endl;
	}

	geoid id = mg.grid_begin(lvl)->id();
	std::cout<<"id.dim : "<<id.dim<<", id.lvl : "<<id.lvl<<", id.pos : "<<id.pos<<std::endl;
	std::cout<<"numb.cell_nodes().at(id) : "<<numb.cell_nodes().at(id)[1]<<std::endl;
	std::cout<<"numb.cell_nodes().at(id) : "<<numb.cell_nodes().at( std::next( mg.grid_begin(lvl),2)->id() )[1]<<std::endl;
	//std::cout<<numb.cell_nodes()[id][1]<<std::endl;
	std::cout<<mg.find_cell(id)/*->node_numbers*/<<std::endl;
	//std::cout<<mg.get_cell(id)/*->node_numbers*/<<std::endl;
	tetrahedron<1> &t { mg.get_cell(id) };
	tetrahedron<1> &t_3 { *mg.grid_begin(lvl) };
	
	std::cout<<"t : "<<t.pos()<<", t_3 : "<<t_3.pos()<<std::endl;
	std::cout<<numb.cell_nodes().begin()->second[1]<<std::endl;
	geoid id_2 = numb.cell_nodes().begin()->first;
	std::cout<<id_2.pos<<std::endl;
	std::cout<<numb.cell_nodes().begin()->first.pos<<std::endl; // trying to access first ( which is geoid id )
	std::cout<<numb(t_3).size()<<std::endl;
	std::cout<<numb(t_3)[1]<<std::endl;

	// aceess bary3d
	bary3d bary3d_1(1,1,1,1);
	bary3d bary3d_2(1,0,0);
	bary3d bary3d_3(-0.5,0.5,-0.5,0.5);
	bary3d bary3d_4( t_3.pos().x, t_3.pos().y, t_3.pos().z );
	bary3d bary3d_5(-0.5,-0.5,0.5);
	bary3d bary3d_6(numb(3).x,numb(3).y,numb(3).z);
	
	std::cout<<"Chi : "<<mg.grid_begin(lvl)->Chi(bary3d_2)<<", xi, eta, zeta = "<<bary3d_2.xi()<<", "<<bary3d_2.eta()<<", "<<bary3d_2.zeta()<<std::endl;
	std::cout<<"Chi : "<<mg.grid_begin(lvl)->Chi(bary3d_3)<<", xi, eta, zeta = "<<bary3d_3.xi()<<", "<<bary3d_3.eta()<<", "<<bary3d_3.zeta()<<std::endl;
	std::cout<<"Chi : "<<mg.grid_begin(lvl)->Chi(bary3d_4)<<", xi, eta, zeta = "<<bary3d_4.xi()<<", "<<bary3d_4.eta()<<", "<<bary3d_4.zeta()<<std::endl;
	std::cout<<"Chi : "<<mg.grid_begin(lvl)->Chi(bary3d_5)<<", xi, eta, zeta = "<<bary3d_5.xi()<<", "<<bary3d_5.eta()<<", "<<bary3d_5.zeta()<<std::endl;
	std::cout<<"Chi : "<<mg.grid_begin(lvl)->Chi(bary3d_6)<<", xi, eta, zeta = "<<bary3d_6.xi()<<", "<<bary3d_6.eta()<<", "<<bary3d_6.zeta()<<std::endl;
	std::cout<<bary3d_1.z0<<std::endl;
	std::cout<<numb(1).x<<std::endl;
	std::cout<<
		"Chi_bary(0,0,0) : "<<mg.grid_begin(lvl)->Chi(bary3d(0,0,0))<<", "<<
		"Chi_bary(1,0,0) : "<<mg.grid_begin(lvl)->Chi(bary3d(1,0,0))<<", "<<
		"Chi_bary(0,1,0) : "<<mg.grid_begin(lvl)->Chi(bary3d(0,1,0))<<", "<<
		"Chi_bary(0,0,1) : "<<mg.grid_begin(lvl)->Chi(bary3d(0,0,1))<<", "<<std::endl;
	for ( int i=0; i<numb(*mg.grid_begin(lvl)).size(); ++i )
		std::cout<<
			"node "<<i<<" : "<<numb(numb(*mg.grid_begin(lvl))[i])<<", ";
	std::cout<<std::endl;

	grid_function<point,1,1> grid_test( numb ); // use this one for solving vector Poisson equation
	std::cout<<"grid_test : "<<grid_test(numb(3))<<", eval : "<<grid_test.eval(numb(3))<<", eval_ : "<<grid_test.eval({0.5,0.5,0.5})<<std::endl; 
	

	int test_133_ = test_133( kk, nnel );std::cout<<test_133_<<std::endl;

	



	}
}

