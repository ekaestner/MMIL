
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _CONSTANTS_AND_TYPEDEFS_H_
#define _CONSTANTS_AND_TYPEDEFS_H_

#include <stdio.h>
#include <math.h>
#include <iostream.h>

// Put everything in the kev namespace
namespace defines
{

	/********/
	/* Misc */
	/********/
	// Converts a number in radians to a number in degrees
	const double    TO_DEG_D         = 57.2957795131;

	// Converts a number in radians to a number in degrees
	const float     TO_DEG_F         = 57.2957795131f;

	// Converts a number in degrees to a number in radians
	const double    TO_RAD_D         = 0.01745329252;

	// Converts a number in degrees to a number in radians
	const float     TO_RAD_F         = 0.01745329252f;

	// Zero - Why?
	const double    ZERO_D           = 1e-6;

	// Zero - Why?
	const float     ZERO_F           = 1e-6f;

	// A huge number - not nessesarily the largest for a 64 bit number
	const double    INFINITY_D       = 1e20;

	// A huge number - not nessesarily the largest for a 32 bit number
	const float     INFINITY_F       = 1e20f;

	// Powell tau d - in 64bit precision
	const double    POWELL_TAU_D     = 0.6180339887499;

	// Powell tau d - in 32bit precision
	const float     POWELL_TAU_F     = 0.6180339887499f;

	// log of Powell tau d - in 64bit precision
	const double    LOG_POWELL_TAU_D	= -0.20898764025;

	// log of Powell tau d - in 32bit precision
	const float     LOG_POWELL_TAU_F	= -0.20898764025f;

	/*************/
	/* Log stuff */
	/*************/
	// e
	const float     E_F           = 2.7182818284590452354f;

	// e
	const double    E_D           = 2.7182818284590452354;

	// Log (base 10) 2E
	const double    LOG2E_D       = 1.4426950408889634074;

	// Log (base 10) 2E
	const float     LOG2E_F       = 1.4426950408889634074f;

	// Log (base 10) 10E
	const double    LOG10E_D      = 0.43429448190325182765;

	// Log (base 10) 10E
	const float     LOG10E_F      = 0.43429448190325182765f;

	// Natural Log of 2
	const float     LN2_F         = 0.69314718055994530942f;

	// Natural Log of 2
	const double    LN2_D         = 0.69314718055994530942;

	// Natural Log of 10
	const float     LN10_F        = 2.30258509299404568402f;

	// Natural Log of 10
	const double    LN10_D        = 2.30258509299404568402;

	/*****************/
	/* Trig/pi stuff */
	/*****************/
	// PI - 64 bit precision
	const double    PI_D          = 3.14159265358979323846264338327950288419716939937510;

	// PI - 32 bit precision
	const float     PI_F          = 3.14159265358979323846f;

	// 2 * PI - 32 bit precision
	const float     TWO_PI_F        = 6.28318530718f;

	// 2 * PI - 64 bit precision
	const double    TWO_PI_D        = 6.28318530718;

	// 1 / PI - 64 bit precision
	const double    INV_PI_D      = 0.31830988618379067154;

	// 1 / PI - 32 bit precision
	const float     INV_PI_F      = 0.31830988618379067154f;

	// 1 / (2*PI) - 64 bit precision
	const double    INV_TWO_PI_D    = 0.159154943092;

	// 1 / (2*PI) - 32 bit precision
	const float     INV_TWO_PI_F    = 0.159154943092f;

	// PI / 2 - 64 bit precision
	const double    PI_OVER_TWO_D   = 1.57079632679489661923;

	// PI / 2 - 32 bit precision
	const float     PI_OVER_TWO_F   = 1.57079632679489661923f;

	// PI / 4  - 64 bit precision
	const double    PI_OVER_FOUR_D   = 0.78539816339744830962;

	// PI / 4  - 32 bit precision
	const float     PI_OVER_FOUR_F   = 0.78539816339744830962f;

	// 2 / PI  - 64 bit precision
	const double    TWO_OVER_PI_D   = 0.63661977236758134308;

	// 2 / PI  - 32 bit precision
	const float     TWO_OVER_PI_F   = 0.63661977236758134308f;

	// 2 / sqrt(PI) - 64 bit precision
	const double    TWO_OVER_SQRT_PI_D = 1.12837916709551257390;

	// 2 / sqrtf(PI) - 32 bit precision
	const float     TWO_OVER_SQRT_PI_F = 1.12837916709551257390f;

	// (2 * PI) / 3 - 64 bit precision
	const double    TWO_PI_OVER_THREE_D = 2.0943951023933334;

	// (2 * PI) / 3 - 32 bit precision
	const float     TWO_PI_OVER_THREE_F = 2.0943951023933334f;

	/**********/
	/* Square */
	/* Root   */
	/* Stuff  */
	/**********/

	// Same as sqrt(2.0)
	const double    SQRT_TWO_D      = 1.41421356237309504880;

	// Same as sqrtf(2.0f)
	const float     SQRT_TWO_F      = 1.41421356237309504880f;

	// Same as 1.0 / sqrt(2.0)
	const double    INV_SQRT_TWO_D  = 0.70710678118654752440;

	// Same as 1.0f / sqrtf(2.0f)
	const float     INV_SQRT_TWO_F  = 0.70710678118654752440f;

	// Same as sqrt(3.0)
	const double    SQRT_THREE_D      = 1.73205080757;

	// Same as sqrtf(3.0f)
	const float     SQRT_THREE_F      = 1.73205080757f;

	// Same as 1.0 / sqrt(3.0)
	const double    INV_SQRT_THREE_D	 = 0.577350269189;

	// Same as 1.0f / sqrtf(3.0f)
	const float     INV_SQRT_THREE_F	 = 0.577350269189f;

	/********/
	/* Bool */
	/********/
	// True
	const bool XzTrue               = true;

	// False
	const bool XzFalse	             = false;

	/*********/
	/* Color */
	/*********/
	// should move this someday to the color class
	const float UNDEFINED_HUE_COLOR = 6.9696969696969f;

	/*******************************/
	/* Standard Data Type Typedefs */
	/*******************************/
	//: Bool
	typedef bool					XzBool;

	//: char
	typedef char					XzChar;

	//: unsigned char
	typedef unsigned char			XzUchar;

	//: unsigned short
	typedef unsigned short			XzUshort;

	//: short
	typedef	short					XzShort;

	//: integer
	typedef int						XzInt;

	//: unsigned integer
	typedef unsigned int			XzUint;

	//: long
	typedef	long					XzLong;

	//: unsigned long
	typedef unsigned long			XzUlong;

	//: float
	typedef float					XzFloat;

	//: double
	typedef double					XzDouble;

	/***********************************/
	/* Unambiguous Data Type Typedefs  */
	/***********************************/
	//: an 8bit value
	typedef char					XzInt8;

	//: an unsigned 8bit integer
	typedef unsigned char			XzUint8;

	//: a 16bit value
	typedef short					XzInt16;

	//: an unsigned 16bit value
	typedef unsigned short			XzUint16;

	//: a 32bit integer
	typedef int						XzInt32;

	//: an unsigned 32bit integer
	typedef unsigned int			XzUint32;

	//: 32bit floating point number
	typedef float					XzFloat32;

	//: a 64bit floating point number
	typedef double					XzFloat64;

	//: a 64bit integer
	typedef long					XzInt64;

	//: an unsigned 64bit integer
	typedef unsigned long			XzUint64;


	/*****************************/
	/* Other Data Type Typedefs  */
	/*****************************/
	//: One byte - an 8 bit unsigned number
	typedef	unsigned char			XzByte;

	//: One word - a 16 bit unsigned number
	typedef	unsigned short			XzWord;

	//: a double word - a 32 bit unsigned number
	typedef	unsigned int			XzDword;
};
using namespace defines;

namespace kev
{
	/*******************************/
	/* math and assorted functions */
	/*******************************/
	inline double SIN( const double& arg )
	{
		return ::sin( arg );
	}
	
	inline float SIN( const float& arg )
	{
		return ::sinf( arg );
	}
	
	inline double SQRT( const double& arg )
	{
		return ::sqrt( arg );
	}
	
	inline float SQRT( const float& arg )
	{
		return ::sqrtf( arg );
	}
	
	inline double COS( const double& arg)
	{
		return ::cos( arg );
	}
	
	inline float COS( const float& arg )
	{
		return ::cosf( arg );
	}

   inline double TAN( const double& arg )
	{
		return ::tan( arg );
	}
	
	inline float TAN( const float& arg )
	{
		return ::tanf( arg );
	}
	
   
	// Compute the factorial
	// give - an object who's type has operator++, operator=, operator<=, and operator*= defined.
	//        it should be a single valued scalar type such as an int, float, double etc....
	// NOTE: This could be faster with a lookup table, but then wouldn't work templated : kevin
	template<class T>
	inline T FACTORIAL(T rhs)
	{
		T lhs = 1;
    
		for( T x = 1; x <= rhs; ++x )
		{
			lhs *= x;
		}
    
		return lhs;
	}

	//: Absolute value
	// give - any single value type that can be 
	//        compared to a scalar such as int, float, double, complex, short, etc...
	// i.e: 1, 1.0, or 1.0f are all valid
	template<class T>
	inline T ABS(const T& num)
		 { return (num < 0) ? -num : num; }
		 

   // test for equality with tolerance...
   template <class dataType>
   inline bool isEqual( const dataType& a, const dataType& b, const dataType& tolerance )
   {
      return ABS( a - b ) <= tolerance;
   }

	/*************/
	/* Max & min */
	/*************/

	// returns the maximum of a, b, and c.
	template<class T> 
	inline const T& max(const T& a, const T& b, const T& c)
		 { return (a>=b) ? ((a>=c)?a:c) : ((b>=c)?b:c); }

	// returns the minimum of a, b, and c.
	template<class T> 
	inline const T& min(const T& a, const T& b, const T& c)
		 { return (a<=b) ? ((a<=c)?a:c) : ((b<=c)?b:c); }

	// returns the maximum of a, b
	template<class T> 
	inline const T& max(const T& a, const T& b)
		 { return (a > b) ? a : b; }

	// returns the minimum of a, b
	template<class T> 
	inline const T& min(const T& a, const T& b)
		 { return (a < b) ? a : b; }

	//: [something] to the power of 2
	template<class T> 
	inline const T& POW2(const T& num)
		{ return num * num; }

	//: [something] to the power of 3
	template<class T> 
	inline T POW3(const T& num)
		{ return POW2(num) * num; }

	//: [something] to the power of 4
	template<class T> 
	inline T POW4(const T& num)
		{ return POW3(num) * num; }

	//: [something] to the power of 5
	template<class T> 
	inline T POW5(const T& num) 
		{ return POW4(num) * num; }

	// Round the floating point value up or down 
	// result - rounds up if value is >= x.5
	//          rounds down if value is < x.5
	template<class T> 
	inline T round( T value )
	{
		return static_cast<T>( static_cast<int>( value + static_cast<T>(0.5) ) );
	}

   // Linear Interpolation
   template<class DataType>
   inline void Lerp(const DataType& from, 
		   const DataType& to, 
		   const float& lerp, 
		   DataType& result )
   {
       DataType offset = to - from;
       result = from + offset*lerp;
   }
   
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
}; // end of namespace kev



#endif
