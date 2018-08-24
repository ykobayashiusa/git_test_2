#ifndef MISC_RUNTIME_INFO_H
#define MISC_RUNTIME_INFO_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace misc
{


inline
int thread_num()
{
	#ifdef _OPENMP
	return omp_get_thread_num();
	#else
	return 0;
	#endif
}

inline
int num_threads()
{
	#ifdef _OPENMP
	return omp_get_num_threads();
	#else
	return 1;
	#endif
}

inline
int max_threads()
{
	#ifdef _OPENMP
	return omp_get_max_threads();
	#else
	return 1;
	#endif
}

inline
int thead_limit()
{
	#ifdef _OPENMP
	return omp_get_thread_limit();
	#else
	return 1;
	#endif
}

}

#endif

