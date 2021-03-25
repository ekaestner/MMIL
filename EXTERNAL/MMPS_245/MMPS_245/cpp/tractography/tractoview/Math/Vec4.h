
//////////////////////////////////////////////////////////////////
//
////
#ifndef VECTOR_3DHOMOGENEOUS_32BIT_INCLUDED
#define VECTOR_3DHOMOGENEOUS_32BIT_INCLUDED

#include "Vec3.h"

template <class Type>
class Vec4
{
protected:

	//union 
	//{
		Type _v[4];
       		//struct _Vec4 { Type _x, _y, _z, _w; };
       		struct _Vec4  
		{
			static Type _x; 
			static Type _y; 
			static Type _z; 
			static Type _w;
		} _Vec4; 
	//} ;

public:

	Vec4(){}
	Vec4 ( const Vec4<Type>& vec );
	Vec4 ( const Vec3<Type>& vec );
	Vec4 ( const Type vec[4] );
	Vec4 ( Type v0, Type v1, Type v2, Type v3 );


	//: get the absolute value of the vector
	void				absolute( Vec4<Type>& absVec );

	Type				dot ( const Vec4 &vec ) const;

	bool					equals ( const Vec4 &vec, Type tolerence ) const;

	Type				getDistance ( const Vec4<Type>& vec ) const;
	Type				getDistanceSquared ( const Vec4<Type>& vec ) const;
	Vec3<Type>					getRealPart() const;
	void					getRealPart ( Vec3<Type>& vec ) const;
	Type				getRealDistance ( const Vec4<Type>& vec ) const;
	Type				getRealDistance ( const Vec3<Type>& vec ) const;
	Type				getRealDistanceSquared ( const Vec4<Type>& vec ) const;
	Type				getRealDistanceSquared ( const Vec3<Type>& vec ) const;
	void				get( Type &v0, Type &v1, Type &v2, Type &v3 ) const;

	Type				length() const;

	Type				normalize();

	inline Type* 			data() { return _v; }
	inline const Type* 		data() const { return _v; }
	
	inline Type&		operator[]( int i ) { return _v[i]; }
	inline const Type&	operator[]( int i ) const { return _v[i]; }
   
   Vec4<Type>& 				operator*=( Type scalar );
	Vec4<Type>& 				operator/=( Type scalar );
	Vec4<Type>& 				operator+=( const Vec4<Type>& vec );
	Vec4<Type>& 				operator+=( const Vec3<Type>& vec );
	Vec4<Type>& 				operator+=( const Type arrayOfScalars[4] );
	Vec4<Type>& 				operator-=( const Vec4<Type>& vec );
	Vec4<Type>& 				operator-=( const Vec3<Type>& vec );
	Vec4<Type>& 				operator=( const Vec4<Type>& vec );
	Vec4<Type>& 				operator=( const Vec3<Type>& vec );
	Vec4<Type>& 				operator=( const Type arrayOfScalars[4] );
	Vec4<Type>				operator-() const;
	Vec4<Type>				operator*( Type scalar ) const;
	Vec4<Type>				operator/( Type scalar ) const;
	Vec4<Type>				operator+( const Vec4<Type>& vec ) const;
 	Vec4<Type>				operator+( const Vec3<Type>& vec ) const;
	Vec4<Type>				operator-( const Vec4<Type>& vec ) const;
	Vec4<Type>				operator-( const Vec3<Type>& vec ) const;
	bool					operator==( const Vec4<Type>& vec ) const;
	bool					operator!=( const Vec4<Type>& vec ) const;
	
	//template <class Type>
	//friend ostream &			operator<<( ostream &out, const Vec4<Type>& vec );
	
	//template <class Type>
	//friend istream &			operator>>( istream &in, Vec4<Type>& vec );

	void					setLength( Type value );

	void					set( const Vec4<Type>& vec );
	void					set( const Vec3<Type>& vec );
	void					set( const Type vec[4] );
	void					set( Type v0, Type v1, Type v2, Type v3 );
};

//  Constructor.
template <class Type>
inline Vec4<Type>::Vec4( const Vec4<Type>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
	_v[3] = vec[3];
}

//  Constructor.
template <class Type>
inline Vec4<Type>::Vec4( const Vec3<Type>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
	_v[3] = 1.0f;
}

//  Constructor.
template <class Type>
inline Vec4<Type>::Vec4( const Type v[4] )
{
	_v[0] = v[0];
	_v[1] = v[1];
	_v[2] = v[2];
	_v[3] = v[3];
}

//  Constructor.
template <class Type>
inline Vec4<Type>::Vec4( Type v0, Type v1, Type v2, Type v3 )
{
	_v[0] = v0;
	_v[1] = v1;
	_v[2] = v2;
	_v[3] = v3;
}

//  Set the value.
template <class Type>
inline void Vec4<Type>::set( const Vec4 &vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
	_v[3] = vec[3];
}

//  Set the value.
template <class Type>
inline void Vec4<Type>::set( const Type v[4] )
{
	_v[0] = v[0];
	_v[1] = v[1];
	_v[2] = v[2];
	_v[3] = v[3];
}

//  Set the value.
template <class Type>
inline void Vec4<Type>::set( Type v0, Type v1, Type v2, Type v3 )
{
	_v[0] = v0;
	_v[1] = v1;
	_v[2] = v2;
	_v[3] = v3;
}

//: get the absolute value of the vector
//  Make all the vector's components positive.
template <class Type>
inline void Vec4<Type>::absolute( Vec4<Type>& absVec )
{
	absVec._v[0] = ABS ( _v[0] );
	absVec._v[1] = ABS ( _v[1] );
	absVec._v[2] = ABS ( _v[2] );
	absVec._v[3] = ABS ( _v[3] );
}


//  Return the dot product.
template <class Type>
inline Type Vec4<Type>::dot( const Vec4<Type>& vec ) const
{
	return 
		_v[0]* vec[0] +
		_v[1]* vec[1] +
		_v[2]* vec[2] +
		_v[3]* vec[3];
}

//  See if the vectors are equal within the given tolerance.
template <class Type>
inline bool Vec4<Type>::equals( const Vec4<Type>& vec, Type tolerance ) const
{
	return ( ABS ( _v[0] - vec._v[0] ) <= tolerance &&
		  ABS ( _v[1] - vec._v[1] ) <= tolerance &&
		  ABS ( _v[2] - vec._v[2] ) <= tolerance &&
		  ABS ( _v[3] - vec._v[3] ) <= tolerance );
}
		

//  Get the value.
template <class Type>
inline void Vec4<Type>::get( Type &v0, Type &v1, Type &v2, Type &v3 ) const
{
	v0 = _v[0];
	v1 = _v[1];
	v2 = _v[2];
	v3 = _v[3];
}

//  Get the distance.
template <class Type>
inline Type Vec4<Type>::getDistance( const Vec4<Type>& vec ) const
{
	return SQRT( ( _v[0] - vec[0] )* ( _v[0] - vec[0] ) +
				   ( _v[1] - vec[1] )* ( _v[1] - vec[1] ) +
				   ( _v[2] - vec[2] )* ( _v[2] - vec[2] ) +
				   ( _v[3] - vec[3] )* ( _v[3] - vec[3] ) );
}

//  Get the squared distance.
template <class Type>
inline Type Vec4<Type>::getDistanceSquared( const Vec4<Type>& vec ) const
{
	return ( _v[0] - vec[0] )* ( _v[0] - vec[0] ) +
		   ( _v[1] - vec[1] )* ( _v[1] - vec[1] ) +
		   ( _v[2] - vec[2] )* ( _v[2] - vec[2] ) +
		   ( _v[3] - vec[3] )* ( _v[3] - vec[3] );
}

//  Get the real value from the homogeneous one.
template <class Type>
inline Vec3<Type> Vec4<Type>::getRealPart() const
{
	Type invW = 1.0f / _v[3];

	return Vec3<Type> ( _v[0]* invW, 
					 _v[1]* invW,
					 _v[2]* invW );
}

//  Get the real value from the homogeneous one.
template <class Type>
inline void Vec4<Type>::getRealPart( Vec3<Type>& vec ) const
{
	Type invW = 1.0f / _v[3];

	vec[0] = _v[0]* invW;
	vec[1] = _v[1]* invW;
	vec[2] = _v[2]* invW;
}

//  Get the real distance from this vector to the given one.
template <class Type>
inline Type Vec4<Type>::getRealDistance( const Vec4<Type>& vec ) const
{
	// Get the real part of this vector.

	Type vA[3];
	Type invW = 1.0f / _v[3];

	vA[0] = _v[0]* invW;
	vA[1] = _v[1]* invW;
	vA[2] = _v[2]* invW;

	// Get the real part of the given vector.

	Type vB[3];
	invW = 1.0f / vec[3];

	vB[0] = vec[0]* invW;
	vB[1] = vec[1]* invW;
	vB[2] = vec[2]* invW;

	// Return the distance between them.

	return SQRT( ( vA[0] - vB[0] )* ( vA[0] - vB[0] ) +
				   ( vA[1] - vB[1] )* ( vA[1] - vB[1] ) +
				   ( vA[2] - vB[2] )* ( vA[2] - vB[2] ) );
}

//  Get the real distance from this vector to the given one.
template <class Type>
inline Type Vec4<Type>::getRealDistance( const Vec3<Type>& vec ) const
{
	// Get the real part of this vector.

	Type v[3];
	Type invW = 1.0f / _v[3];

	v[0] = _v[0]* invW;
	v[1] = _v[1]* invW;
	v[2] = _v[2]* invW;

	// Return the distance between them.

	return SQRT( ( v[0] - vec[0] ) * ( v[0] - vec[0] ) +
				   ( v[1] - vec[1] ) * ( v[1] - vec[1] ) +
				   ( v[2] - vec[2] ) * ( v[2] - vec[2] ) );
}

//  Return the length.
template <class Type>
inline Type Vec4<Type>::length() const
{
	return SQRT( _v[0]* _v[0] + 
				   _v[1]* _v[1] + 
				   _v[2]* _v[2] + 
				   _v[3]* _v[3] );
}

//  Normalize, return the length prior to normalization.
template <class Type>
inline Type Vec4<Type>::normalize()
{
	Type length = SQRT( _v[0]* _v[0] + 
						   _v[1]* _v[1] +
						   _v[2]* _v[2] +
						   _v[3]* _v[3] );

	Type invLength = 1.0f / length;

	_v[0]*= invLength;
	_v[1]*= invLength;
	_v[2]*= invLength;
	_v[3]*= invLength;

	return length;
}

//  Multiply all the components by the value.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator*=( Type value )
{
	_v[0]*= value;
	_v[1]*= value;
	_v[2]*= value;
	_v[3]*= value;

	return*this;
}

//  Divide all the components by the value.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator/=( Type value )
{
	Type invValue = 1.0f / value;

	_v[0]*= invValue;
	_v[1]*= invValue;
	_v[2]*= invValue;
	_v[3]*= invValue;

	return*this;
}

//  Add the vector to this one.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator+=( const Vec4<Type>& vec )
{
	_v[0] += vec[0];
	_v[1] += vec[1];
	_v[2] += vec[2];
	_v[3] += vec[3];

	return*this;
}

//  Add the vector to this vector. It becomes this vector.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator+=( const Vec3<Type>& vec )
{
	// Get the real part of this vector.

	Type v[3];
	Type invW = 1.0f / _v[3];

	v[0] = _v[0]* invW;
	v[1] = _v[1]* invW;
	v[2] = _v[2]* invW;

	// Add the vector to it.
	
	v[0] += vec[0];
	v[1] += vec[1];
	v[2] += vec[2];

	// Put the weight back.

	_v[0] = _v[3]* v[0];
	_v[1] = _v[3]* v[1];
	_v[2] = _v[3]* v[2];

	return*this;
}

//  Subtract the vector from this one.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator-=( const Vec4<Type>& vec )
{
	_v[0] -= vec[0];
	_v[1] -= vec[1];
	_v[2] -= vec[2];
	_v[3] -= vec[3];

	return*this;
}

//  Subtract the vector from this vector. It becomes this vector.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator-=( const Vec3<Type>& vec )
{
	// Get the real part of this vector.

	Type v[3];
	Type invW = 1.0f / _v[3];

	v[0] = _v[0]* invW;
	v[1] = _v[1]* invW;
	v[2] = _v[2]* invW;

	// Subtract the vector from it.
	
	v[0] -= vec[0];
	v[1] -= vec[1];
	v[2] -= vec[2];

	// Assign it to this vector.

	_v[0] = _v[3]* v[0];
	_v[1] = _v[3]* v[1];
	_v[2] = _v[3]* v[2];

	return*this;
}

//  Assign this vector.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator=( const Vec4<Type>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
	_v[3] = vec[3];

	return*this;
}

//  Assign this vector.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator=( const Vec3<Type>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
	_v[3] = 1.0f;

	return*this;
}

//	Assign this vector.
template <class Type>
inline Vec4<Type>& Vec4<Type>::operator=( const Type v[4] )
{
	_v[0] = v[0];
	_v[1] = v[1];
	_v[2] = v[2];
	_v[3] = v[3];

	return*this;
}

//  Return the negative of this vector.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator-() const
{
	return Vec4<Type> ( -_v[0], 
					 -_v[1], 
					 -_v[2], 
					 -_v[3] );
}

//  Return the component-wise product with the given value.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator*( Type scalar ) const
{
	return Vec4<Type> ( _v[0] * scalar, 
					 _v[1]* scalar, 
					 _v[2]* scalar, 
					 _v[3]* scalar );
}

//  Return the component-wise division with the given value.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator/( Type value ) const
{
	Type invValue = 1.0f / value;

	return Vec4<Type> ( _v[0]* invValue, 
					 _v[1]* invValue, 
					 _v[2]* invValue, 
					 _v[3]* invValue );
}

//  Return the vector sum.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator+( const Vec4<Type>& vec ) const
{
	return Vec4<Type> ( _v[0] + vec._v[0], 
				_v[1] + vec._v[1], 
				_v[2] + vec._v[2], 
				_v[3] + vec._v[3] );
}

//  4D Point + 3D vector = 4D point.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator+( const Vec3<Type>& vec ) const
{
	// Get the real part of the point.

	Type p[3];
	Type invW = 1.0f / _v[3];

	p[0] = _v[0]* invW;
	p[1] = _v[1]* invW;
	p[2] = _v[2]* invW;

	// Add the vector to it.
	
	p[0] += vec[0];
	p[1] += vec[1];
	p[2] += vec[2];

	// Put the weight back.

	return Vec4<Type>( _v[3] * p[0], 
			   _v[3] * p[1], 
			   _v[3] * p[2], 
			   _v[3] );
}

//  Return the vector difference.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator-( const Vec4<Type>& vec ) const
{
	return Vec4<Type>( _v[0] - vec._v[0], 
			   _v[1] - vec._v[1], 
			   _v[2] - vec._v[2], 
			   _v[3] - vec._v[3] );
}

//  4D Point - 3D vector = 4D point.
template <class Type>
inline Vec4<Type> Vec4<Type>::operator-( const Vec3<Type>& vec ) const
{
	// Get the real part of the point.

	Type p[3];
	Type invW = 1.0f / _v[3];

	p[0] = _v[0]* invW;
	p[1] = _v[1]* invW;
	p[2] = _v[2]* invW;

	// Subtract the vector from it.
	
	p[0] -= vec[0];
	p[1] -= vec[1];
	p[2] -= vec[2];

	// Put the weight back.

	return Vec4<Type>(	_v[3]* p[0], 
				_v[3]* p[1], 
				_v[3]* p[2], 
				_v[3] );
}

//  See if they're equal.
template <class Type>
inline bool Vec4<Type>::operator==( const Vec4<Type>& vec ) const
{
	return ( _v[0] == vec._v[0] && 
		_v[1] == vec._v[1] &&
		_v[2] == vec._v[2] &&
		_v[3] == vec._v[3] );
}

//  See if they're not equal.
template <class Type>
inline bool Vec4<Type>::operator!=( const Vec4<Type>& vec ) const
{
	return ( _v[0] != vec._v[0] || 
		_v[1] != vec._v[1] ||
		_v[2] != vec._v[2] ||
		_v[3] != vec._v[3] );
}
/*
//  Output operator.
template <class Type>
inline ostream &operator<<( ostream &out, const Vec4<Type>& vec )
{
	out << vec[0] << " " 
		<< vec[1] << " " 
		<< vec[2] << " " 
		<< vec[3];

	return out;
}

//  Input operator.
template <class Type>
inline istream &operator>>( istream &in, Vec4<Type>& vec )
{
	in >> vec[0] 
	   >> vec[1] 
	   >> vec[2] 
	   >> vec[3];

	return in;
}
*/
//  Set the length.
template <class Type>
inline void Vec4<Type>::setLength( Type newLength )
{
	Type oldLength = SQRT( _v[0]* _v[0] + 
				_v[1]* _v[1] + 
				_v[2]* _v[2] + 
				_v[3]* _v[3] );
	
	Type factor = newLength / oldLength;

	_v[0] *= factor;
	_v[1] *= factor;
	_v[2] *= factor;
	_v[3] *= factor;
}





#endif
