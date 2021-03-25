#ifndef KEVIN_ASSERTIONS
#define KEVIN_ASSERTIONS

#include <assert.h>
#include <iostream.h>
//#include <stdio.h>
//#include "Log.h"   // for now, I have the output going to the logger... TODO: develop a better output/logging/iostream/console paradigm
/*
namespace kevnAsserts
{
	//: Sanity check a pointer.
	//#define assertPointer(pointer)\
	//    assert( pointer != NULL                  );\
	//    assert( (int)pointer != (int)0xcdcdcdcd  );\
	//    assert( (int)pointer != (int)0xdddddddd );\
	//    assert( (int)pointer > 1000             )\

   //: Alert
	//  basically pops an informational message to the user if "verify" is false
	//  if expression "verify" is false, then "text" is logged
	inline void Alert( const bool& verify, const char* const text )
	{
		if (verify == false)
		{
			cerr<<text<<"\n"<<flush;
		}
	}
      
	//: Assert
	//  if expression "verify" is false, then "text" is logged, and code asserts only in DEBUG builds
	#if defined(_DEBUG) || defined(DEBUG)
		#define Assert( verify, text )					\
			if (verify == false)						\
			{											\
				cerr<<"ASSERT: "<<text<<"\n"<<flush;          \
			}											\
			assert( verify );							\
			((void)0)
	#else
		#define  Assert( verify, text ) ((void)0)
	#endif

	//: Verify
	//  if expression "verify" is false, then "text" is logged, and code asserts
	#define Verify( verify, text )						\
		if (verify == false)							\
		{												\
			cerr<<"VERIFY: "<<text<<"\n"<<flush;           \
		}												\
		assert( verify );								\
		((void)0)

	//: Fatal
	//  if expression "verify" is false, then "text" is logged, and current thread aborts
	#define Fatal( verify, text )						\
		if (verify == false)							\
		{												\
			cerr<<"FATAL: "<<text<<"\n"<<flush;  \
			exit(0);									\
		}												\
		((void)0)

	//#define assert(arg) \
	//	Use Assert, Verify or Fatal instead
};

using namespace kevnAsserts;
*/

#endif
