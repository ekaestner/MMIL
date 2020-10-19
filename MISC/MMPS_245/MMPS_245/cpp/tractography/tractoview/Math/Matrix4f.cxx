
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include <assert.h>
#include "Matrix4f.h"

//: Default Constructor
//  NOTE: does no initialization, call identity() to init to identity matrix
Matrix4f::Matrix4f()
{
}

      
//: Copy constructor
Matrix4f::Matrix4f( const Matrix4f& M )
{
	_m[0] = M._m[0]; _m[4] = M._m[4]; _m[8]  = M._m[8];  _m[12] = M._m[12];
	_m[1] = M._m[1]; _m[5] = M._m[5]; _m[9]  = M._m[9];  _m[13] = M._m[13];
	_m[2] = M._m[2]; _m[6] = M._m[6]; _m[10] = M._m[10]; _m[14] = M._m[14];
	_m[3] = M._m[3]; _m[7] = M._m[7]; _m[11] = M._m[11]; _m[15] = M._m[15];
}


//: Construct 4x4 matrix with an array of floats
//  required: you must pass a pointer to an array of 16 floats or more
//  WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
/*
Matrix4f::Matrix4f( const float* m )
{
   assert( false );
	_m[0] = m[0]; _m[4] = m[4]; _m[8]  = m[8];  _m[12] = m[12];
	_m[1] = m[1]; _m[5] = m[5]; _m[9]  = m[9];  _m[13] = m[13];
	_m[2] = m[2]; _m[6] = m[6]; _m[10] = m[10]; _m[14] = m[14];
	_m[3] = m[3]; _m[7] = m[7]; _m[11] = m[11]; _m[15] = m[15];
}
*/
//: Construct 4x4 matrix from 16 floats
Matrix4f::Matrix4f( float a0, float a4, float a8,  float a12,
			float a1, float a5, float a9,  float a13,
			float a2, float a6, float a10, float a14,
			float a3, float a7, float a11, float a15 )
{
	_m[0] = a0; _m[4] = a4; _m[8]  = a8;  _m[12] = a12;
	_m[1] = a1; _m[5] = a5; _m[9]  = a9;  _m[13] = a13;
	_m[2] = a2; _m[6] = a6; _m[10] = a10; _m[14] = a14;
	_m[3] = a3; _m[7] = a7; _m[11] = a11; _m[15] = a15;
}



//: perform an abs (absolute value) function to each matrix cell
void Matrix4f::absolute()
{
	kev::ABS( _m[0] );
	kev::ABS( _m[1] );
	kev::ABS( _m[2] );
	kev::ABS( _m[3] );
	
	kev::ABS( _m[4] );
	kev::ABS( _m[5] );
	kev::ABS( _m[6] );
	kev::ABS( _m[7] );
	
	kev::ABS( _m[8] );
	kev::ABS( _m[9] );
	kev::ABS( _m[10] );
	kev::ABS( _m[11] );
	
	kev::ABS( _m[12] );
	kev::ABS( _m[13] );
	kev::ABS( _m[14] );
	kev::ABS( _m[15] );
}

//: copy translation from a matrix
//  copy translation only
void Matrix4f::copyTranslation( const Matrix4f& mat )
{
	_m[12] = mat[12];
	_m[13] = mat[13];
	_m[14] = mat[14];
}
	
//: copy rotation from a matrix
//  copy the upper 3x3 only
void Matrix4f::copyRotation( const Matrix4f& mat )
{
	_m[0] = mat[0];
	_m[1] = mat[1];
	_m[2] = mat[2];

	_m[4] = mat[4];
	_m[5] = mat[5];
	_m[6] = mat[6];

	_m[8] = mat[8];
	_m[9] = mat[9];
	_m[10] = mat[10];
}

//: Get the determinate of the matrix.
float Matrix4f::det() const
{
	return
		_m[0]  * this->det3( 1, 2, 3, 1, 2, 3 ) -
		_m[4]  * this->det3( 1, 2, 3, 0, 2, 3 ) +
		_m[8]  * this->det3( 1, 2, 3, 0, 1, 3 ) -
		_m[12] * this->det3( 1, 2, 3, 0, 1, 2 );
}

//: Get the determinate of the upper left 3 x 3 submatrix.
float Matrix4f::det3() const
{
	return this->det3( 0, 1, 2, 0, 1, 2 );
}

//: Get the determinate of the 3 x 3 submatrix.
float Matrix4f::det3( int r1, int r2, int r3, int c1, int c2, int c3 ) const
{
	// Convert to the correct indices.
	c1 *= 4;
	c2 *= 4;
	c3 *= 4;

	return  _m[r1+c1] *( _m[r2+c2] * _m[r3+c3] - _m[r2+c3] * _m[r3+c2] ) -
		_m[r1+c2] *( _m[r2+c1] * _m[r3+c3] - _m[r2+c3] * _m[r3+c1] ) +
		_m[r1+c3] *( _m[r2+c1] * _m[r3+c2] - _m[r2+c2] * _m[r3+c1] );
}

//: check if 'this' == matrix, within a certain tolerence
bool Matrix4f::equals( const Matrix4f& M, float tolerance ) const
{
	return 
		( fabsf( _m[0]  - M._m[0] )  <= tolerance &&
		  fabsf( _m[1]  - M._m[1] )  <= tolerance &&
		  fabsf( _m[2]  - M._m[2] )  <= tolerance &&
		  fabsf( _m[3]  - M._m[3] )  <= tolerance &&

		  fabsf( _m[4]  - M._m[4] )  <= tolerance &&
		  fabsf( _m[5]  - M._m[5] )  <= tolerance &&
		  fabsf( _m[6]  - M._m[6] )  <= tolerance &&
		  fabsf( _m[7]  - M._m[7] )  <= tolerance &&

		  fabsf( _m[8]  - M._m[8] )  <= tolerance &&
		  fabsf( _m[9]  - M._m[9] )  <= tolerance &&
		  fabsf( _m[10] - M._m[10] ) <= tolerance &&
		  fabsf( _m[11] - M._m[11] ) <= tolerance &&

		  fabsf( _m[12] - M._m[12] ) <= tolerance &&
		  fabsf( _m[13] - M._m[13] ) <= tolerance &&
		  fabsf( _m[14] - M._m[14] ) <= tolerance &&
		  fabsf( _m[15] - M._m[15] ) <= tolerance );
}


//: check if 'this' == matrix, within a certain tolerence
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
bool Matrix4f::equals( const float* m, float tolerance ) const
{
	return 
		( fabsf( _m[0]  - m[0] )  <= tolerance &&
		  fabsf( _m[1]  - m[1] )  <= tolerance &&
		  fabsf( _m[2]  - m[2] )  <= tolerance &&
		  fabsf( _m[3]  - m[3] )  <= tolerance &&

		  fabsf( _m[4]  - m[4] )  <= tolerance &&
		  fabsf( _m[5]  - m[5] )  <= tolerance &&
		  fabsf( _m[6]  - m[6] )  <= tolerance &&
		  fabsf( _m[7]  - m[7] )  <= tolerance &&

		  fabsf( _m[8]  - m[8] )  <= tolerance &&
		  fabsf( _m[9]  - m[9] )  <= tolerance &&
		  fabsf( _m[10] - m[10] ) <= tolerance &&
		  fabsf( _m[11] - m[11] ) <= tolerance &&

		  fabsf( _m[12] - m[12] ) <= tolerance &&
		  fabsf( _m[13] - m[13] ) <= tolerance &&
		  fabsf( _m[14] - m[14] ) <= tolerance &&
		  fabsf( _m[15] - m[15] ) <= tolerance );
}

//: get the inverse of this matrix
bool Matrix4f::getInverse( Matrix4f& M ) const
{
	M.set( _m );
	return M.invert();
}

//: set the scale part of this matrix
//  sets the 3 diagonal cells
void Matrix4f::getScale( float& sx, float& sy, float& sz ) const
{
	sx = _m[0];
	sy = _m[5];
	sz = _m[10];
}

//: get the translation component of this matrix
void Matrix4f::getTranslation( float& tx, float& ty, float& tz ) const
{
	tx = _m[12];
	ty = _m[13];
	tz = _m[14];
}

//: get the transpose of this matrix
void Matrix4f::getTranspose( Matrix4f& M ) const
{
	M.set( _m[0],  _m[1],  _m[2],  _m[3],
				 _m[4],  _m[5],  _m[6],  _m[7],
				 _m[8],  _m[9],  _m[10], _m[11],
				 _m[12], _m[13], _m[14], _m[15] );
}

//: get the rotation and scale part of this matrix
void Matrix4f::getRotScale( Matrix4f& R ) const
{
	R.set(	_m[0], _m[4], _m[8],  0.0f,
			_m[1], _m[5], _m[9],  0.0f,
			_m[2], _m[6], _m[10], 0.0f,
			0.0f,  0.0f,  0.0f, _m[15] );
}

void Matrix4f::getRotationXYZ( float& xRot, float& yRot, float& zRot ) const
{
    float sx;
    const Matrix4f& mat = *this;
    
    zRot = atan2f(mat(0, 1), mat(0, 0));      //(cy*sz)/(cy*cz) - cy falls out
    xRot = atan2f(mat(1, 2), mat(2, 2));      //(sx*cy)/(cx*cy) - cy falls out
    sx   = kev::SIN( xRot );
    yRot = atan2f(-mat(0, 2),(mat(1, 2)/sx) );   // -(-sy)/((sx*cy)/sx)
    
    xRot = TO_DEG_F * xRot;
    yRot = TO_DEG_F * yRot;
    zRot = TO_DEG_F * zRot;
}

//: get a copy of the matrix data.
// Each matrix entry is copied to memory location m
// NOTE: you allocate this memory that m points to
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
void Matrix4f::get( float* m ) const
{
    m[0] = _m[0]; m[4] = _m[4]; m[8]  = _m[8];  m[12] = _m[12];
    m[1] = _m[1]; m[5] = _m[5]; m[9]  = _m[9];  m[13] = _m[13];
    m[2] = _m[2]; m[6] = _m[6]; m[10] = _m[10]; m[14] = _m[14];
    m[3] = _m[3]; m[7] = _m[7]; m[11] = _m[11]; m[15] = _m[15];
}

//: get a copy of the matrix
void Matrix4f::get( float& a0, float& a4, float& a8,  float& a12,
							float& a1, float& a5, float& a9,  float& a13,
							float& a2, float& a6, float& a10, float& a14,
							float& a3, float& a7, float& a11, float& a15 ) const
{
	a0 = _m[0]; a4 = _m[4]; a8  = _m[8];  a12 = _m[12];
	a1 = _m[1]; a5 = _m[5]; a9  = _m[9];  a13 = _m[13];
	a2 = _m[2]; a6 = _m[6]; a10 = _m[10]; a14 = _m[14];
	a3 = _m[3]; a7 = _m[7]; a11 = _m[11]; a15 = _m[15];
}

//: make matrix an identity matrix
void Matrix4f::makeIdentity()
{
	_m[0] = 1.0f; _m[4] = 0.0f; _m[8]  = 0.0f; _m[12] = 0.0f;
	_m[1] = 0.0f; _m[5] = 1.0f; _m[9]  = 0.0f; _m[13] = 0.0f;
	_m[2] = 0.0f; _m[6] = 0.0f; _m[10] = 1.0f; _m[14] = 0.0f;
	_m[3] = 0.0f; _m[7] = 0.0f; _m[11] = 0.0f; _m[15] = 1.0f;
}

//: invert the matrix.
bool Matrix4f::invert()
{
        Matrix4f  a = *this;
	Matrix4f&  b = *this;
	
	int	n = 4;
	int     i, j, k;
        int     r[ 4], c[ 4], row[ 4], col[ 4];
        float  m[ 4][ 4*2], pivot, max_m, tmp_m, fac;

	
        // Initialization
        for (i = 0; i < n; i ++ )
        {
                r[ i] = c[ i] = 0;
                row[ i] = col[ i] = 0;
        }

        // Set working matrix
        for (i = 0; i < n; i++ )
        {
                for (j = 0; j < n; j++ )
                {
                        m[ i][ j] = a[ i * n + j];
                        m[ i][ j + n] =( i == j ) ? 1.0f : 0.0f ;
                }
        }

        // Begin of loop
        for (k = 0; k < n; k++ )
        {
                // Choosing the pivot
                for (i = 0, max_m = 0; i < n; i++ )
                {
                        if (row[ i]  )      continue;
                        for (j = 0; j < n; j++ )
                        {
                                if (col[ j] )          continue;
                                tmp_m = kev::ABS( m[ i][ j]);
                                if (tmp_m > max_m)
                                {
                                        max_m = tmp_m;
                                        r[ k] = i;
                                        c[ k] = j;
                                }
                        }
                }
                row[ r[k] ] = col[ c[k] ] = 1;
                pivot = m[ r[ k] ][ c[ k] ];


                if (fabs( pivot) <= 1e-20)
                {
                        cerr << "*** pivot = %f in mat_inv. ***\n";
                        //exit( 0);
			return false;
                }

                // Normalization
                for (j = 0; j < 2*n; j++ )
                {
                        if (j == c[ k] )
                                m[ r[ k]][ j] = 1.0;
                        else
                                m[ r[ k]][ j] /= pivot;
                }

                // Reduction
                for (i = 0; i < n; i++ )
                {
                        if (i == r[ k] )
                                continue;

                        for (j=0, fac = m[ i][ c[k]]; j < 2*n; j++ )
                        {
                                if (j == c[ k] )
                                        m[ i][ j] = 0.0;
                                else
                                        m[ i][ j] -= fac * m[ r[k]][ j];
                        }
                }
        }

        // Assign inverse to a matrix
        for (i = 0; i < n; i++ )
                for (j = 0; j < n; j++ )
                        row[ i] =( c[ j] == i ) ? r[ j] : row[ i];

        for (i = 0; i < n; i++ )
                for (j = 0; j < n; j++ )
                        b[ i * n +  j] = m[ row[ i]][ j + n];
			
	return true;   // It worked
}

//: invert matrix using alternative algorithm
bool Matrix4f::invertF()
{   
	float d = this->det();

	if (d == 0.0f ) 
	{
	    return false;
	}

	Matrix4f M( Matrix4f::identity() );

	M[0]  =   this->det3( 1, 2, 3, 1, 2, 3 );
	M[4]  = - this->det3( 1, 2, 3, 0, 2, 3 );
	M[8]  =   this->det3( 1, 2, 3, 0, 1, 3 );
	M[12] = - this->det3( 1, 2, 3, 0, 1, 2 );

	M[1]  = - this->det3( 0, 2, 3, 1, 2, 3 );
	M[5]  =   this->det3( 0, 2, 3, 0, 2, 3 );
	M[9]  = - this->det3( 0, 2, 3, 0, 1, 3 );
	M[13] =   this->det3( 0, 2, 3, 0, 1, 2 );

	M[2]  =   this->det3( 0, 1, 3, 1, 2, 3 );
	M[6]  = - this->det3( 0, 1, 3, 0, 2, 3 );
	M[10] =   this->det3( 0, 1, 3, 0, 1, 3 );
	M[14] = - this->det3( 0, 1, 3, 0, 1, 2 );

	M[3]  = - this->det3( 0, 1, 2, 1, 2, 3 );
	M[7]  =   this->det3( 0, 1, 2, 0, 2, 3 );
	M[11] = - this->det3( 0, 1, 2, 0, 1, 3 );
	M[15] =   this->det3( 0, 1, 2, 0, 1, 2 );
	
	M.transpose();
	M /= d;
	this->set( M );

	return true;
}


//: is the matrix an identity matrix?
bool Matrix4f::isIdentity() const
{
	return( _m[0] == 1.0f && _m[4] == 0.0f && _m[8]  == 0.0f && _m[12] == 0.0f &&
			 _m[1] == 0.0f && _m[5] == 1.0f && _m[9]  == 0.0f && _m[13] == 0.0f &&
			 _m[2] == 0.0f && _m[6] == 0.0f && _m[10] == 1.0f && _m[14] == 0.0f &&
			 _m[3] == 0.0f && _m[7] == 0.0f && _m[11] == 0.0f && _m[15] == 1.0f );
}


//: is the matrix an identity matrix within some tolerence?
bool Matrix4f::isIdentity( float tol ) const
{
	return( kev::ABS( _m[0] - 1.0f ) < tol && 
			 kev::ABS( _m[1] ) < tol && 
			 kev::ABS( _m[2] ) < tol && 
			 kev::ABS( _m[3] ) < tol && 
			 kev::ABS( _m[4] ) < tol && 
			 kev::ABS( _m[5] - 1.0f ) < tol && 
			 kev::ABS( _m[6] ) < tol && 
			 kev::ABS( _m[7] ) < tol && 
			 kev::ABS( _m[8] ) < tol && 
			 kev::ABS( _m[9] ) < tol && 
			 kev::ABS( _m[10] - 1.0f ) < tol && 
			 kev::ABS( _m[11] ) < tol && 
			 kev::ABS( _m[12] ) < tol && 
			 kev::ABS( _m[13] ) < tol && 
			 kev::ABS( _m[14] ) < tol && 
			 kev::ABS( _m[15] - 1.0f ) < tol );
}


//: is the matrix the inverse of the matrix M?
bool Matrix4f::isInverse( const Matrix4f& M1 ) const
{
	Matrix4f M2 = *this * M1;
	return M2.isIdentity();
}


//: c = a * b
// required: c, a, and b must each point to 16 floats
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
void Matrix4f::multiply( float* c, const float* a, const float* b )
{
	c[0]  = a[0] * b[0]  + a[4] * b[1]  + a[8]  * b[2]  + a[12] * b[3];
	c[1]  = a[1] * b[0]  + a[5] * b[1]  + a[9]  * b[2]  + a[13] * b[3];
	c[2]  = a[2] * b[0]  + a[6] * b[1]  + a[10] * b[2]  + a[14] * b[3];
	c[3]  = a[3] * b[0]  + a[7] * b[1]  + a[11] * b[2]  + a[15] * b[3];

	c[4]  = a[0] * b[4]  + a[4] * b[5]  + a[8]  * b[6]  + a[12] * b[7];
	c[5]  = a[1] * b[4]  + a[5] * b[5]  + a[9]  * b[6]  + a[13] * b[7];
	c[6]  = a[2] * b[4]  + a[6] * b[5]  + a[10] * b[6]  + a[14] * b[7];
	c[7]  = a[3] * b[4]  + a[7] * b[5]  + a[11] * b[6]  + a[15] * b[7];

	c[8]  = a[0] * b[8]  + a[4] * b[9]  + a[8]  * b[10] + a[12] * b[11];
	c[9]  = a[1] * b[8]  + a[5] * b[9]  + a[9]  * b[10] + a[13] * b[11];
	c[10] = a[2] * b[8]  + a[6] * b[9]  + a[10] * b[10] + a[14] * b[11];
	c[11] = a[3] * b[8]  + a[7] * b[9]  + a[11] * b[10] + a[15] * b[11];

	c[12] = a[0] * b[12] + a[4] * b[13] + a[8]  * b[14] + a[12] * b[15];
	c[13] = a[1] * b[12] + a[5] * b[13] + a[9]  * b[14] + a[13] * b[15];
	c[14] = a[2] * b[12] + a[6] * b[13] + a[10] * b[14] + a[14] * b[15];
	c[15] = a[3] * b[12] + a[7] * b[13] + a[11] * b[14] + a[15] * b[15];
}

//: c = a * b
void Matrix4f::multiply( Matrix4f& c, const Matrix4f& a, const Matrix4f& b )
{
   Matrix4f::multiply( c.data(), a.data(), b.data() );
}

//: this = M * this
void Matrix4f::multLeft( const Matrix4f& M )
{
	float a[16];
	Matrix4f::multiply( a, M._m, _m );
	this->set( a );
}

//: this = this * M
void Matrix4f::multRight( const Matrix4f& M )
{
	float a[16];
	Matrix4f::multiply( a, _m, M._m );
	this->set( a );
}

//: negate every entry
// NOTE: not very efficient, it returns a copy.
Matrix4f Matrix4f::operator-() const
{
	return Matrix4f( -_m[0], -_m[4], -_m[8],  -_m[12],
						-_m[1], -_m[5], -_m[9],  -_m[13],
						-_m[2], -_m[6], -_m[10], -_m[14],
						-_m[3], -_m[7], -_m[11], -_m[15] );
}

//: this = m[16]
// required: m must point to 16 floats
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
/*
Matrix4f& Matrix4f::operator=( const float* m )
{
	_m[0] = m[0]; _m[4] = m[4]; _m[8]  = m[8];  _m[12] = m[12];
	_m[1] = m[1]; _m[5] = m[5]; _m[9]  = m[9];  _m[13] = m[13];
	_m[2] = m[2]; _m[6] = m[6]; _m[10] = m[10]; _m[14] = m[14];
	_m[3] = m[3]; _m[7] = m[7]; _m[11] = m[11]; _m[15] = m[15];

    return *this;
}
*/
//: this = M
Matrix4f& Matrix4f::operator=( const Matrix4f& M )
{
	_m[0] = M._m[0]; _m[4] = M._m[4]; _m[8]  = M._m[8];  _m[12] = M._m[12];
	_m[1] = M._m[1]; _m[5] = M._m[5]; _m[9]  = M._m[9];  _m[13] = M._m[13];
	_m[2] = M._m[2]; _m[6] = M._m[6]; _m[10] = M._m[10]; _m[14] = M._m[14];
	_m[3] = M._m[3]; _m[7] = M._m[7]; _m[11] = M._m[11]; _m[15] = M._m[15];

    return *this;
}

//: this = this * M
Matrix4f& Matrix4f::operator*=( const Matrix4f& M )
{
	float a[16];
	Matrix4f::multiply( a, _m, M._m );
	this->set( a );
	return *this;
}

//: this = this * m[16]
// required: m must point to 16 floats
// WARNING: This operator is dangerous since you could pass a bad pointer, use at your own risk
/*
Matrix4f& Matrix4f::operator*=( const float* m )
{
	float a[16];
	Matrix4f::multiply( _m, m, a );
	this->set( a );
	return *this;
}
*/
      
//: divide every element of 'this' by some scalar
Matrix4f& Matrix4f::operator/=( float value )
{
	float inv = 1.0f / value;

	_m[0] /= inv; _m[4] /= inv; _m[8]  /= inv; _m[12] /= inv;
	_m[1] /= inv; _m[5] /= inv; _m[9]  /= inv; _m[13] /= inv;
	_m[2] /= inv; _m[6] /= inv; _m[10] /= inv; _m[14] /= inv;
	_m[3] /= inv; _m[7] /= inv; _m[11] /= inv; _m[15] /= inv;

	return *this;
}


//: copyOfResult = this * M
// NOTE: for efficiency, try to use *= or multiply...  This function is returning a copy (slow)
Matrix4f Matrix4f::operator*( const Matrix4f& M ) const
{
	Matrix4f M3( Matrix4f::identity() );
	Matrix4f::multiply( M3._m, this->_m, M._m );
	return M3;
}

//: copyOfResult = this * m[16]
// required: m must point to 16 floats
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
// NOTE: for efficiency, try to use *= or multiply...  This function is returning a copy (slow)
/*
Matrix4f Matrix4f::operator*( const float* m ) const
{
	Matrix4f M3( Matrix4f::identity() );
	Matrix4f::multiply( M3._m, this->_m, m );
	return M3;
}
*/
//: test for equality
bool Matrix4f::operator==( const Matrix4f& M2 ) const
{
	return( this->_m[0] == M2[0] && this->_m[4] == M2[4] && this->_m[8]  == M2[8]  && this->_m[12] == M2[12] && 
		this->_m[1] == M2[1] && this->_m[5] == M2[5] && this->_m[9]  == M2[9]  && this->_m[13] == M2[13] && 
		this->_m[2] == M2[2] && this->_m[6] == M2[6] && this->_m[10] == M2[10] && this->_m[14] == M2[14] && 
		this->_m[3] == M2[3] && this->_m[7] == M2[7] && this->_m[11] == M2[11] && this->_m[15] == M2[15] );
}
//: test for equality
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
/*
bool Matrix4f::operator==( const float* m2 ) const
{
	return( this->_m[0] == m2[0] && this->_m[4] == m2[4] && this->_m[8]  == m2[8]  && this->_m[12] == m2[12] && 
		 this->_m[1] == m2[1] && this->_m[5] == m2[5] && this->_m[9]  == m2[9]  && this->_m[13] == m2[13] && 
		this->_m[2] == m2[2] && this->_m[6] == m2[6] && this->_m[10] == m2[10] && this->_m[14] == m2[14] && 
	 	this->_m[3] == m2[3] && this->_m[7] == m2[7] && this->_m[11] == m2[11] && this->_m[15] == m2[15] );
}
*/
      
//: test for un-equality
bool Matrix4f::operator!=( const Matrix4f& M2 ) const
{
	return( this->_m[0] != M2[0] || this->_m[4] != M2[4] || this->_m[8]  != M2[8]  || this->_m[12] != M2[12] || 
		this->_m[1] != M2[1] || this->_m[5] != M2[5] || this->_m[9]  != M2[9]  || this->_m[13] != M2[13] || 
		this->_m[2] != M2[2] || this->_m[6] != M2[6] || this->_m[10] != M2[10] || this->_m[14] != M2[14] || 
		this->_m[3] != M2[3] || this->_m[7] != M2[7] || this->_m[11] != M2[11] || this->_m[15] != M2[15] );
}

//: test for un-equality
// WARNING: This function is dangerous since you could pass a bad pointer, use at your own risk
/*
bool Matrix4f::operator!=( const float* m2 ) const
{
	return( this->_m[0] != m2[0] || this->_m[4] != m2[4] || this->_m[8]  != m2[8]  || this->_m[12] != m2[12] || 
		this->_m[1] != m2[1] || this->_m[5] != m2[5] || this->_m[9]  != m2[9]  || this->_m[13] != m2[13] || 
		this->_m[2] != m2[2] || this->_m[6] != m2[6] || this->_m[10] != m2[10] || this->_m[14] != m2[14] || 
		this->_m[3] != m2[3] || this->_m[7] != m2[7] || this->_m[11] != m2[11] || this->_m[15] != m2[15] );
}
*/
//: Rotate this matrix about an axis
void Matrix4f::rotate( const float& rad, const float& x, const float& y, const float& z )
{
	Matrix4f M;
	M.makeRotation( rad, x, y, z );
	this->multRight( M );
}

//: Rotate the matrix about the x-axis.
void Matrix4f::rotateX( float angle )
{
	Matrix4f R;
	R.makeRotationX( angle );
	this->multRight( R );
}

//: Rotate the matrix about the y-axis.
void Matrix4f::rotateY( float angle )
{
	Matrix4f R;
	R.makeRotationY( angle );
	this->multRight( R );
}

//: Rotate the matrix about the z-axis.
void Matrix4f::rotateZ( float angle )
{
	Matrix4f R;
	R.makeRotationZ( angle );
	this->multRight( R );
}

void Matrix4f::scale( float sx, float sy, float sz )
{
	Matrix4f S;
	S.makeScale( sx, sy, sz );
	this->multRight( S );
}

void Matrix4f::setScale( float sx, float sy, float sz )
{
	_m[0]  = sx;
	_m[5]  = sy;
	_m[10] = sz;
}

void Matrix4f::makeScale( float sx, float sy, float sz )
{
	this->makeIdentity();
	this->setScale( sx, sy, sz );
}



//: Create a Look-at Matrix
// like a gluLookAt, basically sets your camera position.
void Matrix4f::makeLookAt(   const float& eyex, 
			    const float& eyey, 
			    const float& eyez, 
			    const float& centerx, 
			    const float& centery, 
			    const float& centerz, 
			    const float& upx, 
			    const float& upy, 
			    const float& upz )
{
    // use Vec3s, for code reuse...
    Vec3<float> eye(eyex, eyey, eyez);
    Vec3<float> center(centerx, centery, centerz);
    Vec3<float> up(upx, upy, upz);
    up.normalize();
    
    this->identity();
    
    Vec3<float> u, v, n;
    
    n = eye - center;
    n.normalize();
    
    v = up - center;
    v -= n * v.dot(n);
    v.normalize();
    
    v = up - n * up.dot(n);
    v.normalize();
    
    u = v.cross(n);
    
    _m[0] = u[0];   //0,0
    _m[4] = u[1];   //1,0
    _m[8] = u[2];   //2,0
    _m[1] = v[0];   //0,1
    _m[5] = v[1];   //1,1
    _m[9] = v[2];   //2,1
    _m[2] = n[0];   //0,2
    _m[6] = n[1];   //1,2
    _m[10] = n[2];  //2,2
    
    _m[12] = -u.dot(eye);   //3,0
    _m[13] = -v.dot(eye);   //3,1
    _m[14] = -n.dot(eye);   //3,2
}

//: Set to a projection matrix
//  Like a glFrustum, set a projection matrix.
void Matrix4f::makeFrustum( const float& left, const float& right, 
			    const float& bottom, const float& top, 
			    const float& Near, const float& Far)
{
    _m[0] =( 2.0f * Near ) /( right - left );
    _m[1] = 0.0f;
    _m[2] = 0.0f;
    _m[3] = 0.0f;
    
    _m[4] = 0.0f;
    _m[5] =( 2.0f * Near ) /( top - bottom );
    _m[6] = 0.0f;
    _m[7] = 0.0f;
    
    _m[8] =( right + left ) /( right - left );
    _m[9] =( top + bottom ) /( top - bottom );
    _m[10] = -( Far + Near ) /( Far - Near );
    _m[11] = -1.0f;
    
    _m[12] = 0.0f;
    _m[13] = 0.0f;
    _m[14] = -( 2.0f * Far * Near ) /( Far - Near );
    _m[15] = 0.0f;
}

// fovy 
//   The field of view angle, in degrees, in the y-direction. 
// aspect 
//   The aspect ratio that determines the field of view in the x-direction. 
//   The aspect ratio is the ratio of x (width) to y (height). 
// zNear 
//   The distance from the viewer to the near clipping plane (always positive). 
// zFar 
//   The distance from the viewer to the far clipping plane (always positive). 
void Matrix4f::makePerspective( float fovy, float aspect, float zNear, float zFar )
{
   assert( zNear > 0 && zFar > zNear && "invalid near and far values" );
   float right, top;
   float theta = fovy * 0.5 * TO_RAD_F;
   float tangentTheta = kev::TAN( theta );

   // tan(theta) = right / zNear
   // top = tan(theta) * zNear
   // right = tan(theta) * zNear * aspect

   top = tangentTheta * zNear;
   right = top * aspect; // aspect determines the fieald of view in the x-axis

   makeFrustum( -right, right, -top, top, zNear, zFar );
}

//: set the twist about an arbitrary axis.
// NOTE: this erases any translation in this matrix
void Matrix4f::makeRotation( const float& rad, const float& x, const float& y, const float& z )
{
    float cosine = cosf( rad );
    float sine = sinf( rad );
    float one_minus_cosine = 1 - cosine;

    // rotation part    
    _m[0] = x * x +( 1.0f - x*x ) * cosine;
    _m[1] = x * y * one_minus_cosine + z * sine;
    _m[2] = x * z * one_minus_cosine - y * sine;
    
    _m[4] = y * x * one_minus_cosine - z * sine;
    _m[5] = y * y +( 1.0f - y*y ) * cosine;
    _m[6] = y * z * one_minus_cosine + x * sine;
    
    _m[8] = z * x * one_minus_cosine + y * sine;
    _m[9] = z * y * one_minus_cosine - x * sine;
    _m[10] = z * z +( 1.0f - z*z ) * cosine;
    
    // translation part
    _m[12] = 0;
    _m[13] = 0;
    _m[14] = 0;

    // bottom row
    _m[3] = 0;
    _m[7] = 0;
    _m[11] = 0;
    _m[15] = 1;
}

//: copy the rotation part of the matrix mat
void Matrix4f::setRotation( const Matrix4f& mat )
{
	// rotation part
	_m[0] = mat[0];
	_m[1] = mat[1];
	_m[2] = mat[2];

	_m[4] = mat[4];
	_m[5] = mat[5];
	_m[6] = mat[6];

	_m[8] = mat[8];
	_m[9] = mat[9];
	_m[10] = mat[10];
}

//: Set the matrix to a rotation matrix defined by the rotation part of M
// NOTE: erases any translation
void Matrix4f::makeRotation( const Matrix4f& mat )
{
	this->setRotation( mat );
	
	// translation part
	_m[12] = 0;
	_m[13] = 0;
	_m[14] = 0;
	
	// bottom row
	_m[3] = 0;
	_m[7] = 0;
	_m[11] = 0;
	_m[14] = 1;
}


void Matrix4f::setRotationX( float angle )
{
	float value = cosf( angle );

	_m[5]  = value;
	_m[10] = value;

	value = sinf( angle );

	_m[9] = - value;
	_m[6] =   value;
}

void Matrix4f::setRotationY( float angle )
{
	float value = cosf( angle );

	_m[0]  = value;
	_m[10] = value;

	value = sinf( angle );

	_m[2] = - value;	// Is this one right?
	_m[8] =   value;
}

void Matrix4f::setRotationZ( float angle )
{
	float value = cosf( angle );

	_m[0] = value;
	_m[5] = value;

	value = sinf( angle );

	_m[4] = - value;
	_m[1] =   value;
}

void Matrix4f::makeRotationX( float angle )
{
        this->makeIdentity();
	this->setRotationX( angle );
}
void Matrix4f::makeRotationY( float angle )
{
        this->makeIdentity();
	this->setRotationY( angle );
}
void Matrix4f::makeRotationZ( float angle )
{
        this->makeIdentity();
	this->setRotationZ( angle );
}

void Matrix4f::makeRotationXYZ( const float& xRot, const float& yRot, const float& zRot )
{
    Matrix4f& mat = *this;
    
    float sx = kev::SIN( TO_RAD_F * xRot );  
    float cx = kev::SIN( TO_RAD_F * xRot );
    
    float sy = kev::SIN( TO_RAD_F * yRot );  
    float cy = kev::SIN( TO_RAD_F * yRot );
    
    float sz = kev::SIN( TO_RAD_F * zRot );  
    float cz = kev::SIN( TO_RAD_F * zRot );

   // Derived by multiplying the matrices by hand
   // Z*X*Y
   mat(0, 0) = cy*cz - sx*sy*sz;    mat(1, 0) = -cx*sz;     mat(2, 0) = sy*cz + sx*cy*sz;    mat(3, 0) = 0.0f;
   mat(0, 1) = cy*sz + sx*sy*cz;    mat(1, 1) = cx*cz;      mat(2, 1) = sy*sz - sx*cy*cz;    mat(3, 1) = 0.0f;
   mat(0, 2) = -cx*sy;              mat(1, 2) = sx;         mat(2, 2) = cx*cy;               mat(3, 2) = 0.0f;
   mat(0, 3) = 0.0f;                mat(1, 3) = 0.0f;       mat(2, 3) = 0.0f;               mat(3, 3) = 1.0f;

   //zeroClamp();
}

void Matrix4f::setTranslation( float tx, float ty, float tz )
{
	_m[12] = tx;
	_m[13] = ty;
	_m[14] = tz;
}

void Matrix4f::makeTranslation( float tx, float ty, float tz )
{
	this->makeIdentity();
	this->setTranslation( tx, ty, tz );
}


void Matrix4f::set( const Matrix4f& M )
{
	_m[0] = M._m[0]; _m[4] = M._m[4]; _m[8]  = M._m[8];  _m[12] = M._m[12];
	_m[1] = M._m[1]; _m[5] = M._m[5]; _m[9]  = M._m[9];  _m[13] = M._m[13];
	_m[2] = M._m[2]; _m[6] = M._m[6]; _m[10] = M._m[10]; _m[14] = M._m[14];
	_m[3] = M._m[3]; _m[7] = M._m[7]; _m[11] = M._m[11]; _m[15] = M._m[15];
}

//: set the matrix
void Matrix4f::set( const float* m )
{
	_m[0] = m[0]; _m[4] = m[4];	_m[8]  = m[8];	_m[12] = m[12];
	_m[1] = m[1]; _m[5] = m[5];	_m[9]  = m[9];	_m[13] = m[13];
	_m[2] = m[2]; _m[6] = m[6];	_m[10] = m[10];	_m[14] = m[14];
	_m[3] = m[3]; _m[7] = m[7];	_m[11] = m[11];	_m[15] = m[15];
}

//: set the matrix with 16 floats
void Matrix4f::set( float a0, float a4, float a8,  float a12,
			float a1, float a5, float a9,  float a13,
			float a2, float a6, float a10, float a14,
			float a3, float a7, float a11, float a15 )
{
	_m[0] = a0; _m[4] = a4; _m[8]  = a8;  _m[12] = a12;
	_m[1] = a1; _m[5] = a5; _m[9]  = a9;  _m[13] = a13;
	_m[2] = a2; _m[6] = a6; _m[10] = a10; _m[14] = a14;
	_m[3] = a3; _m[7] = a7; _m[11] = a11; _m[15] = a15;
}

void Matrix4f::translate( float tx, float ty, float tz )
{
	Matrix4f T;
	T.makeTranslation( tx, ty, tz );
	this->multRight( T );
}



//: get the transpose of m
void Matrix4f::transpose( const float* m, float* a )
{
	a[0] = m[0];  a[4] = m[1];  a[8]  = m[2];  a[12] = m[3];
	a[1] = m[4];  a[5] = m[5];  a[9]  = m[6];  a[13] = m[7];
	a[2] = m[8];  a[6] = m[9];  a[10] = m[10]; a[14] = m[11];
	a[3] = m[12]; a[7] = m[13]; a[11] = m[14]; a[15] = m[15];
}


//: make this matrix equal to it's transpose
void Matrix4f::transpose()
{
	// Make a transposed copy(don't do the diagonals).

	float a[16];
	                a[4]  = _m[1];	a[8]  = _m[2];	a[12] = _m[3];
	a[1]  = _m[4];                  a[9]  = _m[6];	a[13] = _m[7];
	a[2]  = _m[8];  a[6]  = _m[9];                  a[14] = _m[11];
	a[3]  = _m[12];	a[7]  = _m[13];	a[11] = _m[14];

	// Copy back into this matrix.

	                _m[4] = a[4];   _m[8] = a[8];	_m[12] = a[12];
	_m[1] = a[1];	                _m[9] = a[9];   _m[13] = a[13];
	_m[2] = a[2];	_m[6] = a[6];                   _m[14] = a[14];
	_m[3] = a[3];	_m[7] = a[7];	_m[11] = a[11];
}

//////////////////////////////////////////////////////
// These all depend on Vec3 or Vec4
//////////////////////////////////////////////////////
#include "Vec3.h"
#include "Vec4.h"

//: Set to a Look-at Matrix
// like a gluLookAt, basically sets your camera position.
void Matrix4f::makeLookAt( const Vec3<float>& eye, const Vec3<float>& center, const Vec3<float>& up)
{
    this->makeLookAt( eye[0], eye[1], eye[2], center[0], center[1], center[2], up[0], up[1], up[2] );
}

void Matrix4f::scale( const Vec3<float>& s )
{
	Matrix4f scale_matrix;
	scale_matrix.makeScale( s );
	this->multRight( scale_matrix );
}

void Matrix4f::setScale( const Vec3<float>& s )
{
	_m[0]  = s[0];
	_m[5]  = s[1];
	_m[10] = s[2];
}
void Matrix4f::makeScale( const Vec3<float>& s )
{
	this->makeIdentity();
	this->setScale( s );
}

void Matrix4f::makeRotation( const float& rad, const Vec3<float>& axis)
{
    this->makeRotation( rad, axis[0], axis[1], axis[2] );
}

void Matrix4f::rotate( const float& rad, const Vec3<float>& axis)
{
	Matrix4f M;
	M.makeRotation( rad, axis );
	this->multRight( M );
}

void Matrix4f::setTranslation( const Vec3<float>& t )
{
	_m[12] = t[0];
	_m[13] = t[1];
	_m[14] = t[2];
}
void Matrix4f::makeTranslation( const Vec3<float>& t )
{
	this->makeIdentity();
	this->setTranslation( t );
}

void Matrix4f::translate( const Vec3<float>& t )
{
	Matrix4f T;
	T.makeTranslation( t );
	this->multRight( T );
}

void Matrix4f::getTranslation( Vec3<float>& t ) const
{
	t[0] = _m[12];
	t[1] = _m[13];
	t[2] = _m[14];
}

void Matrix4f::getScale( Vec3<float>& s ) const
{
	s[0] = _m[0];
	s[1] = _m[5];
	s[2] = _m[10];
}

Vec4<float> Matrix4f::operator*( const Vec4<float>& b ) const
{
	return Vec4<float>( _m[0] * b[0] + _m[4] * b[1] + _m[8]  * b[2] + _m[12] * b[3],
			_m[1] * b[0] + _m[5] * b[1] + _m[9]  * b[2] + _m[13] * b[3],
			_m[2] * b[0] + _m[6] * b[1] + _m[10] * b[2] + _m[14] * b[3],
			_m[3] * b[0] + _m[7] * b[1] + _m[11] * b[2] + _m[15] * b[3] );
}
Vec3<float> Matrix4f::operator*( const Vec3<float>& b ) const
{
	float a[4];

	a[3] = 1.0f /( _m[3] * b[0] + _m[7] * b[1] + _m[11] * b[2] + _m[15] );
	a[2] = a[3] *( _m[2] * b[0] + _m[6] * b[1] + _m[10] * b[2] + _m[14] );
	a[1] = a[3] *( _m[1] * b[0] + _m[5] * b[1] + _m[9]  * b[2] + _m[13] );
	a[0] = a[3] *( _m[0] * b[0] + _m[4] * b[1] + _m[8]  * b[2] + _m[12] );

	return Vec3<float>( a );
}

/////////////////////////////////////////////////
//  Members that depend on Quatf
/////////////////////////////////////////////////
#include "Quatf.h"

//: matrix = quat
Matrix4f& Matrix4f::operator=( const Quatf& quat )
{ 
    const int W = 0;
    const int X = 1;
    const int Y = 2;
    const int Z = 3;

    this->_m[0]  = 1.0f  - 2.0f * (quat[Y] * quat[Y] + quat[Z] * quat[Z]);
    this->_m[1]  = 2.0f         * (quat[X] * quat[Y] + quat[W] * quat[Z]);
    this->_m[2]  = 2.0f         * (quat[X] * quat[Z] - quat[W] * quat[Y]);
    this->_m[3]  = 0.0f;

    this->_m[4]  = 2.0f         * (quat[X] * quat[Y] - quat[W] * quat[Z]);
    this->_m[5]  = 1.0f  - 2.0f * (quat[X] * quat[X] + quat[Z] * quat[Z]);
    this->_m[6]  = 2.0f         * (quat[Y] * quat[Z] + quat[W] * quat[X]);
    this->_m[7]  = 0.0f;

    this->_m[8]  = 2.0f         * (quat[X] * quat[Z] + quat[W] * quat[Y]);
    this->_m[9]  = 2.0f         * (quat[Y] * quat[Z] - quat[W] * quat[X]);
    this->_m[10] = 1.0f  - 2.0f * (quat[X] * quat[X] + quat[Y] * quat[Y]);
    this->_m[11] = 0.0f;

    this->_m[12] = 0.0f;
    this->_m[13] = 0.0f;
    this->_m[14] = 0.0f;
    this->_m[15] = 1.0f;
    
    return *this;
}


////////////////////////////////////
//: Friends
////////////////////////////////////
ostream&  operator<<( ostream& out, const Matrix4f& M )
{
	out << M[0] << " " << M[4] << " " << M[8]  << " " << M[12] << endl;
	out << M[1] << " " << M[5] << " " << M[9]  << " " << M[13] << endl;
	out << M[2] << " " << M[6] << " " << M[10] << " " << M[14] << endl;
	out << M[3] << " " << M[7] << " " << M[11] << " " << M[15];

	return out;
}

istream&  operator>>( istream& in, Matrix4f& M )
{
	in >> M[0] >> M[4] >> M[8]  >> M[12];
	in >> M[1] >> M[5] >> M[9]  >> M[13];
	in >> M[2] >> M[6] >> M[10] >> M[14];
	in >> M[3] >> M[7] >> M[11] >> M[15];

	return in;
}

Vec4<float> operator/( const Vec4<float>& v, const Matrix4f& M )
{
	Matrix4f invM( M );
	invM.invert();
	return invM * v;
}

Vec3<float> operator/( const Vec3<float>& v, const Matrix4f& M )
{
	Matrix4f invM( M );
	invM.invert();
	return invM * v;
}


