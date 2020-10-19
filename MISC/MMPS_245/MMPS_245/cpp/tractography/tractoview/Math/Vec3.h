
//////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////
#ifndef VEC3__INCLUDED
#define VEC3__INCLUDED

#include "Utils/Defines.h"

template <class DataType>
class Vec3
{
protected:

	DataType _v[3];

public:

	Vec3(){}
	Vec3( const Vec3<DataType>& vec );
	Vec3( const DataType vec[3] );
	Vec3( const DataType& v0, const DataType& v1, const DataType& v2 );

	void				      makeAbsolute();

	Vec3<DataType>			cross( const Vec3<DataType>& vec ) const;

	DataType				   dot( const Vec3<DataType>& vec ) const;

	bool				      isEqual( const Vec3<DataType>& vec, const DataType& tolerence ) const;

	Vec3<DataType>			getAbs() const;
	void				      getAbs( Vec3<DataType>& vec ) const;
	DataType				   getAngle( const Vec3<DataType>& vec ) const;
	DataType				   getDistance( const Vec3<DataType>& vec ) const;
	DataType				   getDistanceSquared( const Vec3<DataType>& vec ) const;
	DataType				   getRealDistance( const Vec3<DataType>& vec ) const;
	DataType				   getRealDistanceSquared( const Vec3<DataType>& vec ) const;
	void				      get( DataType &v0, DataType &v1, DataType &v2 ) const;

	DataType				   length() const;
	DataType				   lengthSquared() const;
	
	// Linear Interpolation between two vectors.
	static void			   Lerp(const Vec3<DataType>& from, 
					             const Vec3<DataType>& to, 
					             const DataType &lerp, 
					             Vec3<DataType>& lerpedValue );
				    
	DataType			      normalize();

	inline DataType* 			   data() { return _v; }
	inline const DataType* 		data() const { return _v; }
	
	inline DataType &		      operator[]( int i ) { return _v[i]; }
	inline const DataType &		operator[]( int i ) const { return _v[i]; }
	Vec3<DataType>& 		      operator*=( const DataType& value );
	Vec3<DataType>& 		      operator/=( const DataType& value );
	Vec3<DataType>& 		      operator+=( const Vec3<DataType>& vec );
	Vec3<DataType>& 		      operator+=( const DataType v[4] );
	Vec3<DataType>& 		      operator-=( const Vec3<DataType>& vec );
	Vec3<DataType>& 		      operator=( const Vec3<DataType>& vec );
	Vec3<DataType>& 		      operator=( const DataType v[4] );
	Vec3<DataType>& 		      operator=( const DataType& value );
	Vec3<DataType>		         operator-() const;
	Vec3<DataType>		         operator*( const DataType& value ) const;
	Vec3<DataType>		         operator/( const DataType& value ) const;
	Vec3<DataType>		         operator+( const Vec3<DataType>& vecB ) const;
	Vec3<DataType>		         operator-( const Vec3<DataType>& vecB ) const;
	bool				            operator==( const Vec3<DataType>& vecB ) const;
	bool				            operator!=( const Vec3<DataType>& vecB ) const;
	
	//template<class DataType>
	//friend ostream &		operator<<( ostream& out, const Vec3<DataType>& vec );
	
	//template<class DataType>
	//friend istream &	    operator>>( istream& in, Vec3<DataType>& vec );

	void				setLength( const DataType& value );

	void				set( const Vec3<DataType>& vec );
	void				set( const DataType vec[4] );
	void				set( const DataType& v0, const DataType& v1, const DataType& v2 );
};

//  Constructor.
template<class DataType>
inline Vec3<DataType>::Vec3( const Vec3<DataType>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
}

//  Constructor.
template<class DataType>
inline Vec3<DataType>::Vec3( const DataType v[3] )
{
	_v[0] = v[0];
	_v[1] = v[1];
	_v[2] = v[2];
}

//  Constructor.
template<class DataType>
inline Vec3<DataType>::Vec3( const DataType& v0, const DataType& v1, const DataType& v2 )
{
	_v[0] = v0;
	_v[1] = v1;
	_v[2] = v2;
}

//  Set the value.
template<class DataType>
inline void Vec3<DataType>::set( const Vec3<DataType>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];
}

//  Set the value.
template<class DataType>
inline void Vec3<DataType>::set( const DataType v[3] )
{
	_v[0] = v[0];
	_v[1] = v[1];
	_v[2] = v[2];
}

//  Set the value.
template<class DataType>
inline void Vec3<DataType>::set( const DataType& v0, const DataType& v1, const DataType& v2 )
{
	_v[0] = v0;
	_v[1] = v1;
	_v[2] = v2;
}


//  Make all the vector's components positive.
template<class DataType>
inline void Vec3<DataType>::makeAbsolute()
{
	_v[0] = ABS ( _v[0] );
	_v[1] = ABS ( _v[1] );
	_v[2] = ABS ( _v[2] );
}

//  Return the cross product.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::cross( const Vec3<DataType>& vec ) const
{
	return Vec3<DataType> ( _v[1] * vec[2] - _v[2] * vec[1],
					 _v[2] * vec[0] - _v[0] * vec[2],
					 _v[0] * vec[1] - _v[1] * vec[0] );
}

//  Return the dot product.
template<class DataType>
inline DataType Vec3<DataType>::dot ( const Vec3<DataType>& vec ) const
{
	return _v[0] * vec[0] +
		    _v[1] * vec[1] +
		    _v[2] * vec[2];
}

//  See if the vectors are equal within the given tolerance.
template<class DataType>
inline bool Vec3<DataType>::isEqual( const Vec3<DataType>& vec, const DataType& tolerance ) const
{
   return ( kev::isEqual( _v[0], vec._v[0], tolerance ) &&
		      kev::isEqual( _v[1], vec._v[1], tolerance ) &&
		      kev::isEqual( _v[2], vec._v[2], tolerance ) );
}

//  Return the absolute value.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::getAbs() const
{
	return Vec3<DataType>( ABS(_v[0]), ABS(_v[1]), ABS(_v[2]) );
}

//  Return the absolute value.
template<class DataType>
inline void Vec3<DataType>::getAbs ( Vec3<DataType>& vec ) const
{
	vec.set( ABS( _v[0] ), ABS( _v[1] ), ABS( _v[2] ) );
}

//  Get the value.
template<class DataType>
inline void Vec3<DataType>::get( DataType &v0, DataType &v1, DataType &v2 ) const
{
	v0 = _v[0];
	v1 = _v[1];
	v2 = _v[2];
}

//  Return the angle between the two vectors.
template<class DataType>
inline DataType Vec3<DataType>::getAngle ( const Vec3<DataType>& vec ) const
{
    // This is: theta = acosf ( A dot B / |A||B| ).
    return acosf( 
       ( _v[0] * vec[0] + _v[1] * vec[1] + _v[2] * vec[2] ) / 
		 ( kev::SQRT ( _v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2] ) * 
		   kev::SQRT ( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] ) ) 
    );
}

//  Get the distance.
template<class DataType>
inline DataType Vec3<DataType>::getDistance ( const Vec3<DataType>& vec ) const
{
	return kev::SQRT( ( _v[0] - vec[0] ) * ( _v[0] - vec[0] ) +
				    ( _v[1] - vec[1] ) * ( _v[1] - vec[1] ) +
				    ( _v[2] - vec[2] ) * ( _v[2] - vec[2] )   );
}

//  Get the squared distance.
template<class DataType>
inline DataType Vec3<DataType>::getDistanceSquared ( const Vec3<DataType>& vec ) const
{
	return ( _v[0] - vec[0] ) * ( _v[0] - vec[0] ) +
		   ( _v[1] - vec[1] ) * ( _v[1] - vec[1] ) +
		   ( _v[2] - vec[2] ) * ( _v[2] - vec[2] );
}

//  Get the real distance from this vector to the given one.
template<class DataType>
inline DataType Vec3<DataType>::getRealDistance ( const Vec3<DataType>& vec ) const
{
	// Get the real part of this vector.

	DataType vA[2];
	DataType invW = 1.0f / _v[2];

	vA[0] = _v[0] * invW;
	vA[1] = _v[1] * invW;

	// Get the real part of the given vector.

	DataType vB[2];
	invW = 1.0f / vec[2];

	vB[0] = vec[0] * invW;
	vB[1] = vec[1] * invW;

	// Return the distance between them.

	return kev::SQRT ( ( vA[0] - vB[0] ) * ( vA[0] - vB[0] ) +
				   ( vA[1] - vB[1] ) * ( vA[1] - vB[1] ) );
}

//  Return the length. (kev::SQRT of dot product)
template<class DataType>
inline DataType Vec3<DataType>::length() const
{
	return kev::SQRT ( _v[0] * _v[0] + 
                 _v[1] * _v[1] + 
                 _v[2] * _v[2] );
}

//  Return the length squared. (same as dot product)
template<class DataType>
inline DataType Vec3<DataType>::lengthSquared() const
{
	return _v[0] * _v[0] + 
          _v[1] * _v[1] + 
          _v[2] * _v[2];
}

//  Normalize, return the length prior to normalization.
template<class DataType>
inline DataType Vec3<DataType>::normalize()
{
	DataType length = kev::SQRT ( _v[0] * _v[0] + 
				_v[1] * _v[1] +
				_v[2] * _v[2] );

	DataType invLength = 1.0f / length;

	_v[0] *= invLength;
	_v[1] *= invLength;
	_v[2] *= invLength;

	return length;
}

//  Multiply all the components by the value.
template<class DataType>
inline Vec3<DataType>& Vec3<DataType>::operator*=( const DataType& value )
{
	_v[0] *= value;
	_v[1] *= value;
	_v[2] *= value;

	return *this;
}

//  Divide all the components by the value.
template<class DataType>
inline Vec3<DataType>& Vec3<DataType>::operator/=( const DataType& value )
{
	DataType invValue = 1.0f / value;

	_v[0] *= invValue;
	_v[1] *= invValue;
	_v[2] *= invValue;

	return *this;
}

//  Add the vector to this one.
template<class DataType>
inline Vec3<DataType>& Vec3<DataType>::operator+=( const Vec3<DataType>& vec )
{
	_v[0] += vec[0];
	_v[1] += vec[1];
	_v[2] += vec[2];

	return *this;
}

//  Subtract the vector from this one.
template<class DataType>
inline Vec3<DataType>& Vec3<DataType>::operator-=( const Vec3<DataType>& vec )
{
	_v[0] -= vec[0];
	_v[1] -= vec[1];
	_v[2] -= vec[2];

	return *this;
}

//  Assign this vector.
template<class DataType>
inline Vec3<DataType>& Vec3<DataType>::operator=( const Vec3<DataType>& vec )
{
	_v[0] = vec[0];
	_v[1] = vec[1];
	_v[2] = vec[2];

	return *this;
}

//	Assign this vector.
template<class DataType>
inline Vec3<DataType>& Vec3<DataType>::operator=( const DataType v[3] )
{
	_v[0] = v[0];
	_v[1] = v[1];
	_v[2] = v[2];

	return *this;
}

//  Return the negative of this vector.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::operator-() const
{
	return Vec3<DataType>( -_v[0], -_v[1], -_v[2] );
}

//  Return the component-wise product with the given value.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::operator*( const DataType& value ) const
{
    return Vec3<DataType>( _v[0]*value, _v[1]*value, _v[2]*value );
}

//  Return the component-wise division with the given value.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::operator/( const DataType& value ) const
{
	DataType invValue = static_cast<DataType>(1.0) / value;
	return Vec3<DataType> ( _v[0] * invValue, 
					 _v[1] * invValue, 
					 _v[2] * invValue );
}

//  Return the vector sum.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::operator+( const Vec3<DataType>& vecB ) const
{
	return Vec3<DataType> ( _v[0] + vecB[0], 
				_v[1] + vecB[1], 
				_v[2] + vecB[2] );
}

//  Return the vector difference.
template<class DataType>
inline Vec3<DataType> Vec3<DataType>::operator-( const Vec3<DataType>& vecB ) const
{
	return Vec3<DataType>( _v[0] - vecB._v[0], 
				_v[1] - vecB._v[1], 
				_v[2] - vecB._v[2] );
}

//  See if they're equal.
template<class DataType>
inline bool Vec3<DataType>::operator==( const Vec3<DataType>& vecB ) const
{
	return ( _v[0] == vecB._v[0] && 
		 _v[1] == vecB._v[1] &&
		_v[2] == vecB._v[2] );
}

//  See if they're not equal.
template<class DataType>
inline bool Vec3<DataType>::operator!=( const Vec3<DataType>& vecB ) const
{
	return ( _v[0] != vecB._v[0] || 
		_v[1] != vecB._v[1] ||
		_v[2] != vecB._v[2] );
}

/*
//  Output operator.
template<class DataType>
inline ostream& operator<<( ostream &out, const Vec3<DataType>& vec )
{
	out << vec[0] << " " 
		<< vec[1] << " " 
		<< vec[2];

	return out;
}

//  Input operator.
template<class DataType>
inline istream& operator>>( istream &in, Vec3<DataType>& vec )
{
	in >> vec[0] 
	   >> vec[1] 
	   >> vec[2];

	return in;
}
*/

//  Set the length.
template<class DataType>
inline void Vec3<DataType>::setLength( const DataType& newLength )
{
	DataType oldLength = kev::SQRT( _v[0] * _v[0] + 
				_v[1] * _v[1] + 
				_v[2] * _v[2] );
	
	DataType factor = newLength / oldLength;

	_v[0] *= factor;
	_v[1] *= factor;
	_v[2] *= factor;
}

// Linear Interpolation between two vectors.
template<class DataType>
inline void Vec3<DataType>::Lerp(const Vec3<DataType>& from, 
		const Vec3<DataType>& to, 
		const DataType &lerp, 
		Vec3<DataType>& lerpedValue )
{
    Vec3<DataType> offset = to - from;
    lerpedValue = from + offset*lerp;
}

#endif
