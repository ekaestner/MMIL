
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#ifndef _32BIT_4X4_MATRIX
#define _32BIT_4X4_MATRIX

class Quatf;
#include "Vec3.h"
#include "Vec4.h"
#include <iostream.h>

//: 4x4 matrix oriented in columnar order
class  Matrix4f
{
// Constructors
public:
   //: Default Constructor
   //  NOTE: does no initialization, call makeIdentity() to init to identity matrix
   //        or use the copy constructor instead with Matrix4f::identity()
   Matrix4f();

	//: Copy constructor
	Matrix4f( const Matrix4f& M );
	
	//: Construct 4x4 matrix from 16 floats
	Matrix4f( float a0, float a4, float a8,  float a12,
		 float a1, float a5, float a9,  float a13,
		 float a2, float a6, float a10, float a14,
		 float a3, float a7, float a11, float a15 );

// Matrix methods
public:
	//: perform an abs (absolute value) function to each matrix cell
	void			absolute();
	
	//: copy translation from a matrix
	//  copy translation only
	void                    copyTranslation( const Matrix4f& matrix );
	
	//: copy rotation from a matrix
	//  copy the upper 3x3 only
	void                    copyRotation( const Matrix4f& matrix );
	
	//: get a pointer to the matrix data
	float *			data() { return _m; }
	
	//: get a const pointer to the matrix data
	const float *		data() const { return _m; }
	
	//: Get the determinate of the matrix.
	float			det() const;
	
	//: Get the determinate of the upper left 3 x 3 submatrix.
	float			det3() const;
	
	//: Get the determinate of the 3 x 3 submatrix.
	float			det3( int r1, int r2, int r3, int c1, int c2, int c3 ) const;

	//: check if 'this' == matrix, within a certain tolerence
	bool			equals( const Matrix4f& matrix, float tolerence ) const;
	
	//: check if 'this' == matrix, within a certain tolerence
	// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
	bool			equals( const float* m, float tolerence ) const;
	
	//: get the inverse of this matrix
	bool			getInverse( Matrix4f& M ) const;
	
	//: get the rotation and scale part of this matrix
	void			getRotScale( Matrix4f& R ) const;

	void			getRotationXYZ( float& xRot, float& yRot, float& zRot ) const;

	//: set the scale part of this matrix
	//  sets the 3 diagonal cells
	void			getScale( float& sx, float& sy, float& sz ) const;

	//: get the translation component of this matrix
	void			getTranslation( float& tx, float& ty, float& tz ) const;

	//: get the transpose of this matrix
	void			getTranspose( Matrix4f& M ) const;

	//: get a copy of the matrix.
	inline void             get( Matrix4f& M ) const { M = *this; }
	
	//: get a copy of the matrix data.
	// Each matrix entry is copied to memory location m
	// NOTE: you allocate this memory that m points to
	// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
	void			get( float* m ) const;
	
	//: get a copy of the matrix
	void			get( float& a0, float& a4, float& a8,  float& a12,
					float& a1, float& a5, float& a9,  float& a13,
					float& a2, float& a6, float& a10, float& a14,
					float& a3, float& a7, float& a11, float& a15 ) const;
	
	//: make matrix an identity matrix
	void			makeIdentity();

        //: returns an identity matrix.
        inline static const Matrix4f& identity() { static Matrix4f im( 1.0f,0.0f,0.0f,0.0f, 0.0f,1.0f,0.0f,0.0f,
                                                                       0.0f,0.0f,1.0f,0.0f, 0.0f,0.0f,0.0f,1.0f);
                                                   return im; }

	//: invert the matrix.
	bool			invert();
	
	//: invert matrix using alternative algorithm
	bool			invertF();

	//: is the matrix an identity matrix?
	bool			isIdentity() const;
	
	//: is the matrix an identity matrix within some tolerence?
	bool			isIdentity( float tol ) const;
	
	//: is the matrix the inverse of the matrix M?
	bool			isInverse( const Matrix4f& M ) const;
	
	//: c = a * b
	// required: c, a, and b must each point to 16 floats
	// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
	static void             multiply( float* c, const float* a, const float* b );
	
	//: c = a * b
	// required: c, a, and b must each point to 16 floats
	static void             multiply( Matrix4f& c, const Matrix4f& a, const Matrix4f& b );
	
	//: this = M * this
	void                    multLeft( const Matrix4f& M );
	
	//: this = this * M
	void                    multRight( const Matrix4f& M );
	
	//: negate every entry
	// NOTE: not very efficient, it returns a copy.
	Matrix4f                operator-() const;

	//: returns memory element i(out of [0..15])
	inline float&           operator[]( int i ) { return _m[i]; }
	
	//: returns memory element i(out of [0..15])
	inline const float&     operator[]( int i ) const { return _m[i]; }
	
	//: returns element i, j
	inline float&           operator()( const int& i, const int& j ) { return _m[i*4+j]; }
	
	//: returns element i, j
	inline const float&     operator()( const int& i, const int& j ) const { return _m[i*4+j]; }

	//: this = M
	Matrix4f&               operator=( const Matrix4f& M );
	
	//: this = m[16]
	// required: m must point to 16 floats
	// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
	//Matrix4f&               operator=( const float* m );
	
	//: this = this * M
	Matrix4f&               operator*=( const Matrix4f& M );
	
	//: divide every element of 'this' by some scalar
	Matrix4f&               operator/=( float value );

	//: copyOfResult = this * M
	// NOTE: for efficiency, try to use *= or multiply...  This function is returning a copy (slow)
	Matrix4f                operator*( const Matrix4f& MB ) const;
	
	//: test for equality
	bool					operator==( const Matrix4f& MB ) const;
	
	//: test for equality
	// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
	//bool					operator==( const float* m ) const;
	
	//: test for un-equality
	bool					operator!=( const Matrix4f& MB ) const;
	
	//: Rotate this matrix about an axis
	void					rotate( const float& rad, const float& x, const float& y, const float& z );
	
	//: Rotate the matrix about the x-axis.
	void					rotateX( float angle );
	
	//: Rotate the matrix about the y-axis.
	void					rotateY( float angle );
	
	//: Rotate the matrix about the z-axis.
	void					rotateZ( float angle );

	//: Scale the matrix evenly by some scalar.
	void					scale( float x, float y, float z );
	void					setScale( float x, float y, float z );
	void					makeScale( float x, float y, float z );

	//: Set to a Look-at Matrix
	// like a gluLookAt, basically sets your camera position.
	void makeLookAt( const float& eyex, 
			const float& eyey, 
			const float& eyez, 
			const float& centerx, 
			const float& centery, 
			const float& centerz, 
			const float& upx, 
			const float& upy, 
			const float& upz );
	
   
 
	//: Set to a projection matrix
	//  Like a glFrustum, set a projection matrix.
	void makeFrustum( const float& left, const float& right, 
			    const float& bottom, const float& top, 
			    const float& Near, const float& Far);
	
   // fovy 
   //   The field of view angle, in degrees, in the y-direction. 
   // aspect 
   //   The aspect ratio that determines the field of view in the x-direction. 
   //   The aspect ratio is the ratio of x (width) to y (height). 
   // zNear 
   //   The distance from the viewer to the near clipping plane (always positive). 
   // zFar 
   //   The distance from the viewer to the far clipping plane (always positive). 
	void makePerspective( float fovy, float aspect, float zNear, float zFar );

	//: Set the matrix to a rotation matrix defined by the rotation part of M
	void					setRotation( const Matrix4f& M );
	// NOTE: erases any translation
	void					makeRotation( const Matrix4f& M );
	
	//: set the twist about an arbitrary axis.
	// NOTE: this erases any translation in this matrix
	void					makeRotation( const float& rad, const float& x, const float& y, const float& z );
	
	void					setRotationX( float angle );
	void					setRotationY( float angle );
	void					setRotationZ( float angle );

	void					makeRotationX( float angle );
	void					makeRotationY( float angle );
	void					makeRotationZ( float angle );

	void					makeRotationXYZ( const float& xRot, const float& yRot, const float& zRot );
	

	void					setTranslation( float x, float y, float z );
	void					makeTranslation( float x, float y, float z );

	//: set the matrix
	void				        set( const Matrix4f& M );
	
	//: set the matrix with a float pointer
	// required: float pointer must point to a user-allocated array of 16 floats
	// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
	void					set( const float* mat );
	
	//: set the matrix with 16 floats
	void					set( float a0, float a4, float a8,  float a12,
							float a1, float a5, float a9,  float a13,
							float a2, float a6, float a10, float a14,
							float a3, float a7, float a11, float a15 );

	void					translate( float x, float y, float z );
    
	//: make this matrix equal to it's transpose
	void					transpose();
	
	//: get the transpose of m
	static void				transpose( const float* m, float* result );
	
// These all depend on Vec
public:
	//: Set to a Look-at Matrix
	// like a gluLookAt, basically sets your camera position.
	void makeLookAt( const Vec3<float>& eye, const Vec3<float>& center, const Vec3<float>& up);

	void scale( const Vec3<float>& s );
	void setScale( const Vec3<float>& s );
	void makeScale( const Vec3<float>& s );
	void getScale( Vec3<float>& s ) const;

	//void setRotation( const float& rad, const Vec3<float>& axis);
	void makeRotation( const float& rad, const Vec3<float>& axis);
	void rotate( const float& rad, const Vec3<float>& axis);
	//void getRotation( float& rad, Vec3<float>& axis) const;

	void setTranslation( const Vec3<float>& t );
	void makeTranslation( const Vec3<float>& t );
	void translate( const Vec3<float>& t );
	void getTranslation( Vec3<float>& t ) const;

	Vec4<float> operator*( const Vec4<float>& b ) const;
	Vec3<float> operator*( const Vec3<float>& b ) const;
	
	//: matrix = quat
	Matrix4f& operator=( const Quatf& q );
	
//: friends
public:
	friend ostream&  operator<<( ostream& out, const Matrix4f& M );
	friend istream&  operator>>( istream& in, Matrix4f& M );
	friend Vec4<float> operator/( const Vec4<float>& v, const Matrix4f& M );
	friend Vec3<float> operator/( const Vec3<float>& v, const Matrix4f& M );

//: Data members
public:
	// the actual matrix data.
	// NOTE: it's a better idea to use the functions in this class, unless you're writing an external method
	float _m[16];
};

	

#endif
