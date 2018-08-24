namespace bem
{

namespace lattice
{

//////////////
// Order 1. //
//////////////

template <>
constexpr uint size<1>()
{
	return 3;
}

template <>
constexpr positions<1> lattice<1>()
{
	return {{
	         bary( 1, 0, 0 ), // 0
	         bary( 0, 1, 0 ), // 1
	         bary( 0, 0, 1 )  // 2
	       }};
}

//////////////
// Order 2. //
//////////////

template <>
constexpr uint size<2>()
{
	return 6;
}

template <>
constexpr positions<2> lattice<2>()
{
	return {{
	         bary( 1, 0, 0 ),             // 0
	         bary( 0, 1, 0 ),             // 1
	         bary( 0, 0, 1 ),             // 2
	         bary( 1./2., 1./2.,    0. ), // 3
	         bary(    0., 1./2., 1./2. ), // 4
	         bary( 1./2.,    0., 1./2. )  // 5
	       }};
}


//////////////
// Order 3. //
//////////////

template <>
constexpr uint size<3>()
{
	return 10;
}

template <>
constexpr positions<3> lattice<3>()
{
	return {{
	         bary( 1., 0., 0. ),          // 0
	         bary( 0., 1., 0. ),          // 1
	         bary( 0., 0., 1. ),          // 2
	         bary( 2./3., 1./3.,    0. ), // 3
	         bary( 1./3., 2./3.,    0. ), // 4
	         bary(    0., 2./3., 1./3. ), // 5
	         bary(    0., 1./3., 2./3. ), // 6
	         bary( 1./3.,    0., 2./3. ), // 7
	         bary( 2./3.,    0., 1./3. ), // 8
	         bary( 1./3., 1./3., 1./3. )  // 9
	       }};
}


//////////////
// Order 4. //
//////////////

template <>
constexpr uint size<4>()
{
	return 15;
}

template <>
constexpr positions<4> lattice<4>()
{
	return {{
	         bary( 1., 0., 0. ),          // 0
	         bary( 0., 1., 0. ),          // 1
	         bary( 0., 0., 1. ),          // 2
	         bary( 3./4., 1./4.,    0. ), // 3
	         bary( 2./4., 2./4.,    0. ), // 4
	         bary( 1./4., 3./4.,    0. ), // 5
	         bary(    0., 3./4., 1./4. ), // 6
	         bary(    0., 2./4., 2./4. ), // 7
	         bary(    0., 1./4., 3./4. ), // 8
	         bary( 1./4.,    0., 3./4. ), // 9
	         bary( 2./4.,    0., 2./4. ), // 10
	         bary( 3./4.,    0., 1./4. ), // 11
	         bary( 2./4., 1./4., 1./4. ), // 12
	         bary( 1./4., 2./4., 1./4. ), // 13
	         bary( 1./4., 1./4., 2./4. )  // 14
	       }};
}

} // End of namespace lattice.

} // End of namespace bem.

