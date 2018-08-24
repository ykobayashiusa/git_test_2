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


namespace shapefcts1d
{

////////////////////////
// General functions. //
////////////////////////

template <size_t order, typename T, size_t upto>
struct apply_helper
{
	static T eval( const coeffs<order,T> &c, const vals<order>& v )
	{
		return apply_helper<order,T,upto-1>::eval(c,v) + c[upto-1]*v[upto-1];
	}
};

template <size_t order, typename T>
struct apply_helper<order,T,0>
{
	static T eval( const coeffs<order,T>&, const vals<order>& )
	{
		return T {};
	}
};

template <size_t order, typename T> inline
T apply( const coeffs<order,T> &c, const vals<order> &v )
{
	return apply_helper<order,T,num<order>()>::eval(c,v);
}

//////////////
// Order 1. //
//////////////

template <> 
constexpr size_t num<1>()
{
    return 2;
}

template <>
constexpr coeffs<1,real> lattice<1>()
{
	return {{
	            0,
                1,   
           }};
}

template <> constexpr
vals<1> N<1>( real xi )
{
    return {{
                1. - xi,
                xi
           }};
}

template <> constexpr
vals<1> dNdxi<1>( real )
{
    return {{
                -1.,
                 1.,
           }};
}

//////////////
// Order 2. //
//////////////

template <> 
constexpr size_t num<2>()
{
    return 3;
}

template <>
constexpr coeffs<2,real> lattice<2>()
{
	return {{
	            0,
               .5,
                1,   
           }};
}

template <> constexpr
vals<2> N<2>( real xi )
{
    return {{
                2.*(xi-.5)*(xi-1.),
               -4.*xi*(xi-1.),
                2.*xi*(xi-.5)
           }};
}

template <> constexpr
vals<2> dNdxi<2>( real xi )
{
    return {{
                4*xi - 3,
               -8*xi + 4, 
                4*xi - 1,
           }};
}

}



namespace shapefcts2d
{

////////////////////////
// General functions. //
////////////////////////

template <size_t order, typename T, size_t upto>
struct apply_helper
{
	static T eval( const coeffs<order,T> &c, const vals<order>& v )
	{
		return apply_helper<order,T,upto-1>::eval(c,v) + c[upto-1]*v[upto-1];
	}
};

template <size_t order, typename T>
struct apply_helper<order,T,0>
{
	static T eval( const coeffs<order,T>&, const vals<order>& )
	{
		return T {};
	}
};


template <size_t order, typename T> inline
T apply( const coeffs<order,T> &c, const vals<order> &v )
{
	return apply_helper<order,T,num<order>()>::eval(c,v);
}

//////////////
// Order 0. //
//////////////

template <>
constexpr size_t num<0>()
{
	return 1;
}

template <>
constexpr coeffs<0,bary2d> lattice<0>()
{
	return {{ bary2d { 1./3., 1./3., 1./3. } }};
}

template <>
constexpr vals<0> N<0>( bary2d )
{
	return {{ 1 }};
}


template <>
constexpr vals<0> dNdxi<0>( bary2d )
{
	return {{ 0 }};
}

template <>
constexpr vals<0> dNdeta<0>( bary2d )
{
	return {{ 0 }};
}

//////////////
// Order 1. //
//////////////

template <>
constexpr size_t num<1>()
{
	return 3;
}

template <>
constexpr coeffs<1,bary2d> lattice<1>()
{
	return {{
	          bary2d( 1, 0, 0 ), // 0
	          bary2d( 0, 1, 0 ), // 1
	          bary2d( 0, 0, 1 )  // 2
	       }};
}

template <>
constexpr vals<1> N<1>( bary2d p )
{
	return {{ p.z0, p.z1, p.z2 }};
}


template <>
constexpr vals<1> dNdxi<1>( bary2d )
{
	return {{ -1, 1, 0 }};
}

template <>
constexpr vals<1> dNdeta<1>( bary2d )
{
	return {{ 0, -1, 1 }};
}


//////////////
// Order 2. //
//////////////

template <>
constexpr size_t num<2>()
{
	return 6;
}

template <>
constexpr coeffs<2,bary2d> lattice<2>()
{
	return {{
	          bary2d( 1, 0, 0 ),             // 0
	          bary2d( 0, 1, 0 ),             // 1
	          bary2d( 0, 0, 1 ),             // 2
	          bary2d( 1./2., 1./2.,    0. ), // 3
	          bary2d(    0., 1./2., 1./2. ), // 4
	          bary2d( 1./2.,    0., 1./2. )  // 5
	       }};
}

template <>
constexpr vals<2> N<2>( bary2d p )
{
	return {{ p.z0*(2*p.z0-1),
	          p.z1*(2*p.z1-1),
	          p.z2*(2*p.z2-1),
	          4*p.z0*p.z1,
	          4*p.z1*p.z2,
	          4*p.z2*p.z0 }};
}


template <>
constexpr vals<2> dNdxi<2>( bary2d p )
{
	return {{ 1-4*p.z0,
	          4*p.z1-1,
	          0,
	          4*(p.z0-p.z1),
	          4*p.z2,
	         -4*p.z2 }};
} 

template <>
constexpr vals<2> dNdeta<2>( bary2d p )
{
	return {{ 0,
	          1-4*p.z1,
	          4*p.z2-1,
	         -4*p.z0,
	          4*(p.z1-p.z2),
	          4*p.z0 }};
}

}

namespace shapefcts3d
{

////////////////////////
// General functions. //
////////////////////////

template <size_t order, typename T, size_t upto>
struct apply_helper
{
	static T eval( const coeffs<order,T> &c, const vals<order>& v )
	{
		return apply_helper<order,T,upto-1>::eval(c,v) + c[upto-1]*v[upto-1];
	}
};

template <size_t order, typename T>
struct apply_helper<order,T,0>
{
	static T eval( const coeffs<order,T>&, const vals<order>& )
	{
		return T {};
	}
};


template <size_t order, typename T> inline
T apply( const coeffs<order,T> &c, const vals<order> &v )
{
	return apply_helper<order,T,num<order>()>::eval(c,v);
}


//////////////
// Order 0. //
//////////////

template <>
constexpr size_t num<0>()
{
	return 1;
}

template <>
constexpr coeffs<0,bary3d> lattice<0>()
{
	return {{ bary3d(.25,.25,.25,.25) }};
}

template <>
constexpr vals<0> N<0>( bary3d )
{
	return {{ 1 }};
}


template <>
constexpr vals<0> dNdxi<0>( bary3d )
{
	return {{ 0 }};
}

template <>
constexpr vals<0> dNdeta<0>( bary3d )
{
	return {{ 0 }};
}

template <>
constexpr vals<0> dNdzeta<0>( bary3d )
{
	return {{ 0 }};
}

template <>
constexpr coeffs<0,point> dN<0>( bary3d )
{
    return {{ point { 0, 0, 0 } }};
}

//////////////
// Order 1. //
//////////////

template <>
constexpr size_t num<1>()
{
	return 4;
}

template <>
constexpr coeffs<1,bary3d> lattice<1>()
{
	return {{
	         bary3d( 1, 0, 0, 0 ), // 0
	         bary3d( 0, 1, 0, 0 ), // 1
	         bary3d( 0, 0, 1, 0 ), // 2
             bary3d( 0, 0, 0, 1 ), // 3
	       }};
}

template <>
constexpr vals<1> N<1>( bary3d p )
{
	return {{ p.z0, p.z1, p.z2, p.z3 }};
}


template <>
constexpr vals<1> dNdxi<1>( bary3d )
{
	return {{ -1, 1, 0, 0 }};
}

template <>
constexpr vals<1> dNdeta<1>( bary3d )
{
	return {{ -1, 0, 1, 0 }};
}

template <>
constexpr vals<1> dNdzeta<1>( bary3d )
{
	return {{ -1, 0, 0, 1 }};
}

template <>
constexpr coeffs<1,point> dN<1>( bary3d )
{
    return {{
                point { -1, -1, -1 },
                point {  1,  0,  0 },
                point {  0,  1,  0 },
                point {  0,  0,  1 }
           }};
}

//////////////
// Order 2. //
//////////////

template <>
constexpr size_t num<2>()
{
	return 10;
}

template <>
constexpr coeffs<2,bary3d> lattice<2>()
{
	return {{
	         bary3d(  1,  0,  0,  0 ),  // 0
	         bary3d(  0,  1,  0,  0 ),  // 1
	         bary3d(  0,  0,  1,  0 ),  // 2
             bary3d(  0,  0,  0,  1 ),  // 3
             bary3d( .5, .5,  0,  0 ),  // 4
             bary3d(  0, .5, .5,  0 ),  // 5
             bary3d( .5,  0, .5,  0 ),  // 6
             bary3d( .5,  0,  0, .5 ),  // 7
             bary3d(  0, .5,  0, .5 ),  // 8
             bary3d(  0,  0, .5, .5 ),  // 9  
	       }};
}

template <>
constexpr vals<2> N<2>( bary3d p )
{
	return {{ p.z0*(2*p.z0 - 1), // 0
              p.z1*(2*p.z1 - 1), // 1
              p.z2*(2*p.z2 - 1), // 2
              p.z3*(2*p.z3 - 1), // 3
              4*p.z0*p.z1,       // 4
              4*p.z1*p.z2,       // 5
              4*p.z2*p.z0,       // 6
              4*p.z0*p.z3,       // 7
              4*p.z1*p.z3,       // 8
              4*p.z2*p.z3        // 9
            }};
}


template <>
constexpr vals<2> dNdxi<2>( bary3d p )
{
	return {{
                -(4*p.z0 - 1), // 0
                  4*p.z1 - 1 , // 1
                           0 , // 2
                           0 , // 3
              4*(p.z0 - p.z1), // 4
                      4*p.z2 , // 5
                     -4*p.z2 , // 6
                     -4*p.z3 , // 7
                      4*p.z3 , // 8
                           0   // 9
           }};
} 

template <>
constexpr vals<2> dNdeta<2>( bary3d p )
{
	return {{
                -(4*p.z0 - 1), // 0
                           0 , // 1
                  4*p.z2 - 1 , // 2
                           0 , // 3
                     -4*p.z1 , // 4
                      4*p.z1 , // 5
              4*(p.z0 - p.z2), // 6
                     -4*p.z3 , // 7
                           0 , // 8
                      4*p.z3   // 9
           }};
}

template <>
constexpr vals<2> dNdzeta<2>( bary3d p )
{
	return {{ 
                -(4*p.z0 - 1), // 0
                           0 , // 1
                           0 , // 2
                  4*p.z3 - 1 , // 3
                     -4*p.z1 , // 4
                           0 , // 5 
                     -4*p.z2 , // 6
              4*(p.z0 - p.z3), // 7
                      4*p.z1 , // 8
                      4*p.z2   // 9
           }};
}

template <>
constexpr coeffs<2,point> dN<2>( bary3d p )
{
    return
    {{
        point {   -(4*p.z0 - 1),   -(4*p.z0 - 1),    -(4*p.z0 - 1) }, // 0
        point {     4*p.z1 - 1 ,              0 ,               0  }, // 1
        point {              0 ,     4*p.z2 - 1 ,               0  }, // 2
        point {              0 ,              0 ,      4*p.z3 - 1  }, // 3
        point { 4*(p.z0 - p.z1),        -4*p.z1 ,         -4*p.z1  }, // 4
        point {         4*p.z2 ,         4*p.z1 ,               0  }, // 5
        point {        -4*p.z2 , 4*(p.z0 - p.z2),         -4*p.z2  }, // 6
        point {        -4*p.z3 ,        -4*p.z3 ,  4*(p.z0 - p.z3) }, // 7
        point {         4*p.z3 ,              0 ,          4*p.z1  }, // 8
        point {              0 ,         4*p.z3 ,          4*p.z2  }  // 9
    }};
}

}

}

