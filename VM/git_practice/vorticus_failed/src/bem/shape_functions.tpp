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

namespace shapefcts1d
{

////////////////////////
// General functions. //
////////////////////////

template <uint order, typename T, uint upto>
struct apply_helper
{
	static T eval( const coeffs<order,T> &c, const vals<order>& v )
	{
		return apply_helper<order,T,upto-1>::eval(c,v) + c[upto-1]*v[upto-1];
	}
};

template <uint order, typename T>
struct apply_helper<order,T,0>
{
	static T eval( const coeffs<order,T>&, const vals<order>& )
	{
		return T {};
	}
};

template <uint order, typename T> inline
T apply( const coeffs<order,T> &c, const vals<order> &v )
{
	return apply_helper<order,T,num<order>()>::eval(c,v);
}

template <uint order> constexpr
uint num()
{
    return order + 1;
}

//////////////
// Order 1. //
//////////////

template <> constexpr
vals<1> N<1>( real xi )
{
    return {{
                1. - xi,
                xi
           }};
}

template <> constexpr
vals<1> dN<1>( real )
{
    return {{
                -1.,
                 1.,
           }};
}

//////////////
// Order 2. //
//////////////

template <> constexpr
vals<2> N<2>( real xi )
{
    return {{
                2.*(xi-.5)*(xi-1.),
               -4.*xi*(xi-1.),
                2.*xi*(xi-.5)
           }};
}

//////////////
// Order 3. //
//////////////

template <> constexpr
vals<3> N<3>( real xi )
{
    return {{
                -4.5*(xi-1./3.)*(xi-2./3.)*(xi-1.),
                13.5*xi*(xi-2./3.)*(xi-1.),
               -13.5*xi*(xi-1./3.)*(xi-1.),
                 4.5*xi*(xi-1./3.)*(xi-2./3.)
           }};
}

}

namespace shapefcts2d
{

////////////////////////
// General functions. //
////////////////////////

template <uint order, typename T, uint upto>
struct apply_helper
{
	static T eval( const coeffs<order,T> &c, const vals<order>& v )
	{
		return apply_helper<order,T,upto-1>::eval(c,v) + c[upto-1]*v[upto-1];
	}
};

template <uint order, typename T>
struct apply_helper<order,T,0>
{
	static T eval( const coeffs<order,T>&, const vals<order>& )
	{
		return T {};
	}
};


template <uint order, typename T> inline
T apply( const coeffs<order,T> &c, const vals<order> &v )
{
	return apply_helper<order,T,num<order>()>::eval(c,v);
}


//////////////
// Order 1. //
//////////////

template <>
constexpr uint num<1>()
{
	return 3;
}

template <>
constexpr vals<1> N<1>( bary p )
{
	return {{ p.z0, p.z1, p.z2 }};
}


template <>
constexpr vals<1> dNdxi<1>( bary )
{
	return {{ -1, 1, 0 }};
}

template <>
constexpr vals<1> dNdeta<1>( bary )
{
	return {{ 0, -1, 1 }};
}


//////////////
// Order 2. //
//////////////

template <>
constexpr uint num<2>()
{
	return 6;
}

template <>
constexpr vals<2> N<2>( bary p )
{
	return {{ p.z0*(2*p.z0-1),
	          p.z1*(2*p.z1-1),
	          p.z2*(2*p.z2-1),
	          4*p.z0*p.z1,
	          4*p.z1*p.z2,
	          4*p.z2*p.z0 }};
}


template <>
constexpr vals<2> dNdxi<2>( bary p )
{
	return {{ 1-4*p.z0,
	          4*p.z1-1,
	          0,
	          4*(p.z0-p.z1),
	          4*p.z2,
	         -4*p.z2 }};
} 

template <>
constexpr vals<2> dNdeta<2>( bary p )
{
	return {{ 0,
	          1-4*p.z1,
	          4*p.z2-1,
	         -4*p.z0,
	          4*(p.z1-p.z2),
	          4*p.z0 }};
}

//////////////
// Order 3. //
//////////////

template <>
constexpr uint num<3>()
{
	return 10;
}

template <>
constexpr vals<3> N<3>( bary p )
{
	return {{
	           4.5*p.z0*(p.z0-(1.0/3.0))*(p.z0-(2.0/3.0)),
	           4.5*p.z1*(p.z1-(1.0/3.0))*(p.z1-(2.0/3.0)),
	           4.5*p.z2*(p.z2-(1.0/3.0))*(p.z2-(2.0/3.0)),
	           4.5*p.z0*p.z1*(3*p.z0-1),
	           4.5*p.z0*p.z1*(3*p.z1-1),
	           4.5*p.z1*p.z2*(3*p.z1-1),
	           4.5*p.z1*p.z2*(3*p.z2-1),
	           4.5*p.z2*p.z0*(3*p.z2-1),
	           4.5*p.z2*p.z0*(3*p.z0-1),
	           27.*p.z0*p.z1*p.z2
	       }};
}	

}

}

