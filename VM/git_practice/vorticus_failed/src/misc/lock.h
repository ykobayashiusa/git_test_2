#ifndef MISC_LOCK_H
#define MISC_LOCK_H

#ifdef _OPENMP

#include <omp.h>

namespace misc
{

class mutex
{
public:
	mutex();
	~mutex();

	void   lock();
	void unlock();

	mutex( const mutex& rhs ) = delete;
	mutex& operator=( const mutex& rhs ) = delete;

private:
	omp_lock_t l;
};

inline
mutex::mutex()
{
	omp_init_lock( &l );
}

inline
mutex::~mutex()
{
	omp_destroy_lock( &l );
}

inline
void mutex::lock()
{
	omp_set_lock( &l );
}

inline
void mutex::unlock()
{
	omp_unset_lock( &l );
}

}

#else

namespace misc
{

class mutex
{
public:
	void lock()   {}
	void unlock() {}
};

}

#endif


namespace misc
{

class lock
{
public:
	 lock( mutex& mut );
	~lock();

	void unlock();
	void relock();

	lock( const lock& rhs ) = delete;
	lock& operator=( const lock& rhs ) = delete;

private:
	mutex& m;
	bool locked;
};

inline
lock::lock( mutex& mut ):
 m(mut), locked(true)
{
	m.lock();
}

inline
lock::~lock()
{
	unlock();
}

inline
void lock::unlock()
{
	if ( locked )
	{
		m.unlock();
		locked = false;
	}
}

inline
void lock::relock()
{
	if ( ! locked )
	{
		m.lock();
		locked = true;
	}
}

}

#endif

