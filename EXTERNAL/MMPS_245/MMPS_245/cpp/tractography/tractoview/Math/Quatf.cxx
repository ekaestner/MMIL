
//////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Defines.h"

#include <iostream.h>
#include <string.h> // strlen
#include "Vec3.h"

#include "Quatf.h"


////////////////////
// Constructors
////////////////////
Quatf::Quatf()
{
	_s = 1;
	_V[0] = 1;
	_V[1] = _V[2] = 0;
	
}

Quatf::Quatf( const float twistvec[4] )
{
	_s    = cosf( twistvec[0] * .5f );
	_V[0] = sinf( twistvec[0] * .5f ) * twistvec[1];
	_V[1] = sinf( twistvec[0] * .5f ) * twistvec[2];
	_V[2] = sinf( twistvec[0] * .5f ) * twistvec[3];
}

//quat = quat(w,x,y,z)
Quatf::Quatf( const float &w, const float &x, 
                  const float &y, const float &z )

{
	_s    = w;
	_V[0] = x;
	_V[1] = y;
	_V[2] = z;
}

Quatf::Quatf( const Quatf &quat )
{
	_s = quat._s;
	_V = quat._V;
}



//takes a scalar (twist) and a vector for w and [x,y,z];
Quatf::Quatf( const float &s, const Vec3<float> &vec )
{
	_s = s;
	_V = vec;
}
//takes a vector (vec), and makes a 'pure' quaternion.
Quatf::Quatf( const Vec3<float> &vec )
{
	_s    = 0;
	_V[0] = vec[0];
	_V[1] = vec[1];
	_V[2] = vec[2];
}
////////////////////
// Public functions
////////////////////

Quatf Quatf::makeRoll(const float &s)
{
	return Quatf(cosf(s*0.5f), sinf(s*0.5f), 0, 0 );
}

Quatf  Quatf::makePitch(const float & t)
{
	return Quatf(cosf(t*0.5f), 0, sinf(t*0.5f), 0 );
}

Quatf  Quatf::makeYaw(const float & p)
{
	return Quatf(cosf(p*0.5f), 0, 0, sinf(p*0.5f) ); 
}

//takes an angle (twist) in radians and a vector (vec).
Quatf Quatf::setRotation(const float &rad, const Vec3<float> &vec)
{
	float halfRad = rad * .5f;
	float sinHalfRad = sinf(halfRad);
	Vec3<float> vecNormalized = vec;
	vecNormalized.normalize();

	_s    = cosf(halfRad);
	_V[0] = sinHalfRad * vecNormalized[0];
	_V[1] = sinHalfRad * vecNormalized[1];
	_V[2] = sinHalfRad * vecNormalized[2];
	this->unit();
	return *this;
}

//returns twist and dir in the form of rad, vec
void Quatf::getRotation(float &rad, Vec3<float> &vec) const
{
	float oneOverASin2 = 1 / ( asinf(_s) * 2 );
	rad = acosf(_s) * 2;
	vec.set( _V[0] * ( oneOverASin2 ),
	              _V[1] * ( oneOverASin2 ),
		      _V[2] * ( oneOverASin2 ));
	vec.normalize();
}	

////////////////////
// Operators
////////////////////

//quat = quat
Quatf & Quatf::operator= ( const Quatf &q )
{
	_s    = q._s;
	_V[0] = q._V[0];
	_V[1] = q._V[1];
	_V[2] = q._V[2];
	return *this;
}

//quat =  float[4] quaternion
Quatf & Quatf::operator= ( const float quat[4] )
{
	_s    = quat[0];
	_V[0] = quat[1];
	_V[1] = quat[2];
	_V[2] = quat[3];
	return *this;
}

//quat = 'pure' quat
Quatf & Quatf::operator= ( const Vec3<float> &vec )
{
	//construct a 'pure' quaternion
	_s = 0;
	_V = vec;
	return *this;
}
//////////////////////////////////////////////////////////////


Quatf	Quatf::operator* ( const Quatf &q2 )
{
	return Quatf( _s * q2._s - _V.dot(q2._V),  
		        q2._V * _s + _V * q2._s + 
			_V.cross(q2._V) );
}
Vec3<float>	Quatf::operator* ( const Vec3<float> &vec )
{
	Quatf pure(vec);
	Quatf result = this->getInverse() * pure * (*this);
	return Vec3<float>(result[1],result[2],result[3]);
	
}

Quatf	Quatf::operator/ ( const Quatf &q2 )
{
	return  (((Quatf)q2).getInverse()) * (*this);
}

Quatf	Quatf::operator* ( float val )
{
	return Quatf(_s * val, _V * val);
}
Quatf	operator* ( float val, const Quatf &q )
{
	return Quatf(q._s * val, q._V * val);
}
Quatf Quatf::operator/ ( float val )
{
	val = 1 / val;
	return val * (*this);
}

Quatf	Quatf::operator+( const Quatf &q2 )
{
	return Quatf( (_s + q2._s), (_V + q2._V) );
}
Quatf	Quatf::operator-( const Quatf &q2 )
{
	return Quatf( (_s - q2._s), (_V - q2._V) );
}

int Quatf::operator==( const Quatf &q2 )
{
	if( _s == q2._s &&
		_V == q2._V   ) return 1;
	else return 0;
}
int Quatf::operator!= ( const Quatf &q2 )
{
	if( _s == q2._s &&
		_V == q2._V   ) return 0;
	else return 1;
}

ostream & operator<<( ostream &out, Quatf q)
{
	float rad;
	Vec3<float> vec;
	q.getRotation(rad,vec);
	
	out << rad*TO_DEG_F << " deg, " << vec[0] << ", " << vec[1] << ", " << vec[2];
	return out;
}


istream & operator>>( istream &in, Quatf &q )
{
	float deg;
	Vec3<float> vec;
	char *throwaway = (char *) new char[::strlen(" deg, ")];
	
	in >> deg;
	in >> throwaway;
	in >> vec[0] >> vec[1] >> vec[2];

	q.setRotation(deg * TO_RAD_F, vec);
	delete throwaway;
	return in;
}
	
////////////////////////////////////////////////////////////////
//  Linear interpolation between two quaternion positions
//
//  Arguments: Quatf  (first and second quaternion)
//             GLfloat  (interpolation parameter [0..1])
//             Quatf  (resulting quaternion, inbetween)
//  Comments:  Fast but not nearly as smooth as Slerp
////////////////////////////////////////////////////////////////
void Quatf::Lerp(const Quatf &from, const Quatf &to, const float &t, Quatf &res)
{
        Quatf		    to1;
        float           cosom;
        float           scale0, scale1;

        // calc cosine
        cosom = from._V[0] * to._V[0] + 
	        from._V[1] * to._V[1] + 
		from._V[2] * to._V[2] + 
		from._s    * to._s;

        // adjust signs (if necessary)
        if ( cosom < 0.0f )
	{
			//could probably do a to1 = -to;
			to1._V[0] = - to._V[0];
			to1._V[1] = - to._V[1];
			to1._V[2] = - to._V[2];
			to1._s    = - to._s;
        } else  {

			to1._V[0] = to._V[0];
			to1._V[1] = to._V[1];
			to1._V[2] = to._V[2];
			to1._s    = to._s;
        }

	// interpolate linearly
        scale0 = 1.0f - t;
        scale1 = t;

	// calculate final values
	res._V[0] = scale0 * from._V[0] + scale1 * to1._V[0];
	res._V[1] = scale0 * from._V[1] + scale1 * to1._V[1];
	res._V[2] = scale0 * from._V[2] + scale1 * to1._V[2];
	res._s    = scale0 * from._s    + scale1 * to1._s;
}

///////////////////////////////////////////////////////////////////////////////
// Spherical Linear Interpolation Between two Quaternions
// Arguments:	Two Quaternions, blend factor, result quaternion
// Notes:	The comments explain the handling of the special cases.
//		The comments in the magazine were wrong.  Adjust the
//		DELTA constant to play with the LERP vs. SLERP level
///////////////////////////////////////////////////////////////////////////////
void Quatf::Slerp(const Quatf &q1, 
                    const Quatf &q2,
                    const float &slerp, 
                    Quatf &result, 
                    const float &delta )
{
    float omega, cosom, sinom, s0, s1;
    
    // do a dot product to get cos(theta) between the two quats
    cosom = q1._V[0] * q2._V[0] + 
	    q1._V[1] * q2._V[1] + 
	    q1._V[2] * q2._V[2] + 
	    q1._s    * q2._s; 

    // CHECK A COUPLE OF SPECIAL CASES.
    // MAKE SURE THE TWO QUATERNIONS ARE NOT EXACTLY OPPOSITE?
    // (WITHIN A LITTLE SLOP)
    if ((1.0f + cosom) > delta)
    {
	// ARE THEY MORE THAN A LITTLE BIT DIFFERENT? 
	// AVOID A DIVIDED BY ZERO AND LERP IF NOT
	if ( (1.0f - cosom) > delta) 
	{
	    // YES, DO A SLERP
	    omega = acosf(cosom);
	    sinom = sinf(omega);
	    s0 = sinf((1.0f - slerp) * omega) / sinom;
	    s1 = sinf(slerp * omega) / sinom;
	} else {
	    // NOT A VERY BIG DIFFERENCE, DO A LERP
	    s0 = 1.0f - slerp;
	    s1 = slerp;
	}
	result._V[0] = s0 * q1._V[0] + s1 * q2._V[0];
	result._V[1] = s0 * q1._V[1] + s1 * q2._V[1];
	result._V[2] = s0 * q1._V[2] + s1 * q2._V[2];
	result._s    = s0 * q1._s    + s1 * q2._s;
    } else {
	// THE QUATERNIONS ARE NEARLY OPPOSITE SO TO 
	// AVOID A DIVIDED BY ZERO ERROR CALCULATE A
	// PERPENDICULAR QUATERNION AND SLERP THAT
	// DIRECTION
	result._V[0] = -q2._V[1];
	result._V[1] =  q2._V[0];
	result._V[2] = -q2._s;
	result._s    =  q2._V[2];
	
	s0 = sinf((1.0f - slerp) * PI_OVER_TWO_F);
	s1 = sinf(slerp * PI_OVER_TWO_F);
	result._V[0] = s0 * q1._V[0] + s1 * result._V[0];
	result._V[1] = s0 * q1._V[1] + s1 * result._V[1];
	result._V[2] = s0 * q1._V[2] + s1 * result._V[2];
	result._s    = s0 * q1._s    + s1 * result._s;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Spherical Linear Interpolation Between two Quaternions
// Arguments:	Two Quaternions, 
//              blend factor, 
//              result quaternion
// Note: This one will work nicely across angle boundarys
//       such as (170, 0,1,0) and (160, 0,-1,0) which are 
//       actually only 30 degrees apart.
///////////////////////////////////////////////////////////////////////////////
void Quatf::Slerp2(const Quatf &q1, 
		     const Quatf &q2,
		     const float &slerp, 
		     Quatf &result, 
		     const float &delta )
{
    Quatf   q1b;
    float       omega, cosom, sinom, s0, s1;
    
    // USE THE DOT PRODUCT TO GET THE COSINE 
    // OF THE ANGLE BETWEEN THE QUATERNIONS
    cosom = q1._V[0] * q2._V[0] + 
	    q1._V[1] * q2._V[1] + 
	    q1._V[2] * q2._V[2] + 
	    q1._s * q2._s; 

    // MAKE SURE WE ARE TRAVELING ALONG THE SHORTER PATH
    if (cosom < 0.0f)
    {
	// IF WE ARE NOT, REVERSE ONE OF THE QUATERNIONS
	cosom = -cosom;
	q1b._V[0] = - q1._V[0];
	q1b._V[1] = - q1._V[1];
	q1b._V[2] = - q1._V[2];
	q1b._s = - q1._s;
    } else {
	q1b._V[0] = q1._V[0];
	q1b._V[1] = q1._V[1];
	q1b._V[2] = q1._V[2];
	q1b._s = q1._s;
    }

    
    if ((1.0f - cosom) > delta) 
    {
	omega = acosf(cosom);
	sinom = sinf(omega);
	s0 = sinf((1.0f - slerp) * omega) / sinom;
	s1 = sinf(slerp * omega) / sinom;
    } else {
	s0 = 1.0f - slerp;
	s1 = slerp;
    }

    result._V[0] = s0 * q1b._V[0] + s1 * q2._V[0];
    result._V[1] = s0 * q1b._V[1] + s1 * q2._V[1];
    result._V[2] = s0 * q1b._V[2] + s1 * q2._V[2];
    result._s    = s0 * q1b._s    + s1 * q2._s;
}



//////////////////////////////////////////////////////////////////
//
// Convert a set of Euler angles to a Quaternion
// Arguments:	A rotation set of 3 angles, a quaternion to set
// Discussion:  As the order of rotations is important I am
//		using the Quantum Mechanics convention of (p,y,r)
//		a Yaw-Pitch-Roll (y,p,r) system would have to be
//		adjusted.  It is more efficient this way though.
//
//////////////////////////////////////////////////////////////////
void Quatf::setPyrEuler(const float &pitch, 
                          const float &yaw, 
			  const float &roll)
{
    float pOver2, yOver2, rOver2, 
	  cosp, cosy, cosr, 
	  sinp, siny, sinr, 
	  cc, cs, sc, ss;

    // precompute half angles
    pOver2 = pitch * 0.5f;
    yOver2 = yaw * 0.5f;
    rOver2 = roll * 0.5f;
    
    cosp = cosf( pOver2 );
    cosy = cosf( yOver2 );
    cosr = cosf( rOver2 );
    
    sinp = sinf( pOver2 );
    siny = sinf( yOver2 );
    sinr = sinf( rOver2 );

    cc = cosp * cosr;
    cs = cosp * sinr;
    sc = sinp * cosr;
    ss = sinp * sinr;

    _s    = (cosy * cc) + (siny * ss);
    _V[0] = (cosy * sc) - (siny * cs);
    _V[1] = (cosy * ss) + (siny * cc);
    _V[2] = (cosy * cs) - (siny * sc);
    
    // insure the quaternion is normalized
    this->normalize();
}

//////////////////////////////////////////////////////////////////
//
// Convert a set of Euler angles (Y*P*R) to a Quaternion
// Arguments:	A rotation set of 3 angles, a quaternion to set
// Discussion:  a Yaw-Pitch-Roll (y,p,r) system.
//
//////////////////////////////////////////////////////////////////
void Quatf::setYprEuler(const float &y, 
                          const float &p, 
			  const float &r)
{
    float cosr, cosp, cosy, 
          sinr, sinp, siny, 
	  cosp_cosy, sinp_siny;
    
    cosr = cosf(r * 0.5f);
    cosp = cosf(p * 0.5f);
    cosy = cosf(y * 0.5f);
    
    sinr = sinf(r * 0.5f);
    sinp = sinf(p * 0.5f);
    siny = sinf(y * 0.5f);
    
    cosp_cosy = cosp * cosy;
    sinp_siny = sinp * siny;
    
    _s    = (cosr * cosp_cosy)   + (sinr * sinp_siny);
    _V[0] = (sinr * cosp_cosy)   - (cosr * sinp_siny);
    _V[1] = (cosr * sinp * cosy) + (sinr * cosp * siny);
    _V[2] = (cosr * cosp * siny) - (sinr * sinp * cosy);
    
    // insure the quaternion is normalized
    this->normalize();
}

///////////////////////////////////////////////////////////////////////////////
// Convert a set of Euler angles to a Quaternion
// Arguments:	A rotation set of 3 angles, a quaternion to set
// Discussion:  This is a second variation.  It creates a
//		Series of quaternions and multiplies them together
//		It would be easier to extend this for other rotation orders
///////////////////////////////////////////////////////////////////////////////
void Quatf::setPyrEuler2(const float &p, 
                           const float &y, 
			   const float &r)
{
    float   pOver2, yOver2, rOver2;
    Quatf qx, qy, qz;
    
    // precompute half angles
    pOver2 = p * 0.5f;
    yOver2 = y * 0.5f;
    rOver2 = r * 0.5f;

    //make the pitch quat
    qx._s    = cosf(pOver2);
    qx._V[0] = sinf(pOver2); 
    qx._V[1] = 0.0f; 
    qx._V[2] = 0.0f; 
    
    //make the yaw quat
    qy._s    = cosf(yOver2);
    qy._V[0] = 0.0f; 
    qy._V[1] = sinf(yOver2); 
    qy._V[2] = 0.0f; 
    
    //make the roll quat
    qz._s    = cosf(rOver2);
    qz._V[0] = 0.0f; 
    qz._V[1] = 0.0f; 
    qz._V[2] = sinf(rOver2); 
    
    //compose the three in pyr order...
    (*this) = qx * qy * qz;

    // insure the quaternion is normalized
    this->normalize();
}
		 
////////////////////////////////////////
//: members who depend on Matrix4f
////////////////////////////////////////
#include "Matrix4f.h"
 
//: quat = matrix
Quatf& Quatf::operator=( const Matrix4f& m )
{
    float    tr, s;
    float    q[4];
    int      i, j, k;
    int      nxt[3] = {1, 2, 0};

    tr = m[0] + m[5] + m[10];

    // check the diagonal
    if (tr > 0.0f) 
    {
	s = sqrtf (tr + 1.0f);
	this->_s = s * 0.5f;
	s = 0.5f / s;
	
	this->_V[0] = (m[6] - m[9]) * s;
	this->_V[1] = (m[8] - m[2]) * s;
	this->_V[2] = (m[1] - m[4]) * s;
    }
    // diagonal is negative
    else 
    {
	i = 0;
	
	if (m[5]  > m[0]) i = 1;
	if (m[10] > m(i, i)) i = 2;
	
	j = nxt[i];
	k = nxt[j];
	
	s = sqrtf(( m(i, i) - (m(j, j)+m(k, k)) ) + 1.0f );
	
	q[i] = s * 0.5f;
	
	if ( s != 0.0f ) 
	    s = 0.5f / s;

	q[3] = (m(j, k) - m(k, j)) * s;
	q[j] = (m(i, j) + m(j, i)) * s;
	q[k] = (m(i, k) + m(k, i)) * s;

	this->_V[0] = q[0];
	this->_V[1] = q[1];
	this->_V[2] = q[2];
	this->_s    = q[3];
    }
    
    return *this;
}

///////////////////////////////////////////////////////////////
// UNIMPLEMENTED
///////////////////////////////////////////////////////////////
/*

  Name:		gluQuatSetFromAx_EXT

  Action:   Constructs quaternion to rotate from one direction vector to
			another

  Params:   GLfloat (x1, y1, z1 - from vector),
			GLfloat (x2, y2, z2 - to vector), GL_QUAT* (resulting quaternion)

  Returns:  nothing

  Comments: Two vectors have to be UNIT vectors (so make sure you normalize
			them before calling this function
			Some parts are heavily optimized so readability is not so great :(

void gluQuatSetFromAx_EXT(GLfloat x1,GLfloat y1, GLfloat z1,
			GLfloat x2,GLfloat y2, GLfloat z2, GL_QUAT *quat)

{
    GLfloat tx, ty, tz, temp, dist;

    GLfloat cost, len, ss;

    // get dot product of two vectors
    cost = x1 * x2 + y1 * y2 + z1 * z2;

    // check if parallel
    if (cost > 0.99999f) {
	quat->x = quat->y = quat->z = 0.0f;
	quat->w = 1.0f;
	return;
    }
    else if (cost < -0.99999f) {		// check if opposite

	// check if we can use cross product of from vector with [1, 0, 0]
	tx = 0.0;
	ty = x1;
	tz = -y1;

	len = sqrt(ty*ty + tz*tz);

	if (len < DELTA)
	{
		// nope! we need cross product of from vector with [0, 1, 0]
		tx = -z1;
		ty = 0.0;
		tz = x1;
	}

	// normalize
	temp = tx*tx + ty*ty + tz*tz;

    dist = (GLfloat)(1.0 / sqrt(temp));

    tx *= dist;
    ty *= dist;
    tz *= dist;
	
	quat->x = tx;
	quat->y = ty;
	quat->z = tz;
	quat->w = 0.0;

	return;
    }

	// ... else we can just cross two vectors

	tx = y1 * z2 - z1 * y2;
	ty = z1 * x2 - x1 * z2;
	tz = x1 * y2 - y1 * x2;

	temp = tx*tx + ty*ty + tz*tz;

    dist = (GLfloat)(1.0 / sqrt(temp));

    tx *= dist;
    ty *= dist;
    tz *= dist;


    // we have to use half-angle formulae (sin^2 t = ( 1 - cos (2t) ) /2)
	
	ss = (float)sqrt(0.5f * (1.0f - cost));

    tx *= ss;
    ty *= ss;
    tz *= ss;

    // scale the axis to get the normalized quaternion
    quat->x = tx;
    quat->y = ty;
    quat->z = tz;

    // cos^2 t = ( 1 + cos (2t) ) / 2
    // w part is cosine of half the rotation angle
    quat->w = (float)sqrt(0.5f * (1.0f + cost));

}
//

  Name:		gluQuatSquare_EXT

  Action:   Square quaternion

  Params:   GL_QUAT* (q1 * q1 = res)

  Returns:  nothing

  Comments: none

void APIENTRY gluQuatSquare_EXT(GL_QUAT* q1, GL_QUAT* res)
{
	GLfloat  tt;


	tt = 2 * q1->w;
	res->x = tt * q1->x;
	res->y = tt * q1->y;
	res->z = tt * q1->z;
	res->w = (q1->w * q1->w - q1->x * q1->x - q1->y * q1->y - q1->z * q1->z);
}

//

  Name:		gluQuatSqrt_EXT

  Action:   Find square root of a quaternion

  Params:   GL_QUAT* (sqrt(q1) = res)

  Returns:  nothing

  Comments: none

void APIENTRY gluQuatSqrt_EXT(GL_QUAT* q1, GL_QUAT* res)
{
	GLfloat  length, m, r1, r2;
	GL_QUAT r;

	length = sqrt (q1->w * q1->w + q1->x * q1->x + q1->y * q1->y);
	if (length != 0.0)
		length = 1.0 / length;
	else length = 1.0;

	r.x = q1->x * length;
	r.y = q1->z * length;
	r.z = 0.0f;
	r.w = q1->w * length;

	m = 1.0 / sqrt (r.w * r.w + r.x * r.x);
	r1 = sqrt ((1.0 + r.y) * 0.5);
	r2 = sqrt ((1.0 - r.y) * 0.5);

	res->x = sqrt (length) * r2 * r.x * m;
	res->y = sqrt (length) * r1;
	res->z = q1->z;
	res->w = sqrt (length) * r1 * r.w * m;

}

//

  Name:		gluQuatDot_EXT

  Action:   Computes the dot product of two unit quaternions

  Params:   GL_QUAT (first and second quaternion)

  Returns:  (GLfloat) Dot product

  Comments: Quaternion has to be normalized (i.e. it's a unit quaternion)

GLfloat APIENTRY gluQuatDot_EXT(GL_QUAT* q1, GL_QUAT* q2)
{
  return (GLfloat)(q1->w * q2->w + q1->x * q2->x + q1->y * q2->y+q1->z*q2->z);
}

//

  Name:		gluQuatLength_EXT

  Action:   Calculates the length of a quaternion

  Params:   GL_QUAT* (quaternion)

  Returns:  GLfloat (length)

  Comments: none

GLfloat APIENTRY gluQuatLength_EXT(GL_QUAT* q1)
{
  return sqrt (q1->w * q1->w + q1->x * q1->x + q1->y * q1->y + q1->z * q1->z);
}

//
  Name:		gluQuatExp_EXT

  Action:   Calculates exponent of a quaternion

  Params:   GL_QUAT* (Source and destination quaternion)

  Returns:  nothing

  Comments: none

void APIENTRY gluQuatExp_EXT(GL_QUAT* q1, GL_QUAT* q2)
{
	GLfloat  len1, len2;

	len1 = (GLfloat) sqrt (q1->x * q1->x + q1->y * q1->y + q1->z * q1->z);
	if (len1 > 0.0)
		len2 = (GLfloat)sin(len1) / len1;
	else
		len2 = 1.0;

	q2->x = q1->x * len2;
	q2->y = q1->y * len2;
	q2->z = q1->z * len2;
	q2->w = cos (len1);
}


//
  Name:		gluQuatLog_EXT

  Action:   Calculates natural logarithm of a quaternion

  Params:   GL_QUAT* (Source and destination quaternion)

  Returns:  nothing

  Comments: none

void APIENTRY gluQuatLog_EXT(GL_QUAT* q1, GL_QUAT* q2)
{
	GLfloat  length;

	length = sqrt (q1->x * q1->x + q1->y * q1->y + q1->z * q1->z);

	//make sure we do not divide by 0
	if (q1->w != 0.0)
		length = atan (length / q1->w);
	else length = (GLfloat)M_PI/2;

	q2->w = 0.0f;
	q2->x = q1->x * length;
	q2->y = q1->y * length;
	q2->z = q1->z * length;
}

//
  Name:		gluQuatLnDif_EXT

  Action:   Computes the "natural log difference" of two quaternions,
			q1 and q2 as  ln(qinv(q1)*q2)

  Params:   GL_QUAT* (Source quaternions  and a destination quaternion)

  Returns:  nothing

  Comments: none

void APIENTRY gluQuatLnDif_EXT(GL_QUAT *q1, GL_QUAT *q2, GL_QUAT *res)
{

	GL_QUAT inv, dif, temp;
	GLfloat  len, len1, s;

	qt_inverse (a, &inv);
	qt_mul (&inv, b, &dif);
	len = sqrt (dif.x*dif.x + dif.y*dif.y + dif.z*dif.z);
	s = qt_dot (a, b);
	if (s != 0.0) len1 = atan (len / s); else len1 = M_PI/2;
	if (len != 0.0) len1 /= len;
	temp.w = 0.0;
	temp.x = dif.x * len1;
	temp.y = dif.y * len1;
	temp.z = dif.z * len1;
	qt_copy (&temp, out);
}

 */
