
//////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#ifndef QUATERNION_32BIT
#define QUATERNION_32BIT

class Matrix4f;
class Quatf
{
//: Constructors
public:
	//: Default constructor
	Quatf();

	//: copy constructor
	Quatf( const Quatf& quat );

	//takes a scalar( twist in radians) and a vector for axis [x,y,z] of rotation
	Quatf( const float& w, const float& x, const float& y, const float& z );
	
	//takes a scalar( twist in radians) and a vector for axis [x,y,z] of rotation
	Quatf( const float& s, const Vec3<float>& vec ); 
	
	//pure quat, [0, vec]
	Quatf( const Vec3<float>& vec ); 

	//from degrees/direction to a quaternion, does conversion
	Quatf( const float twistvec[4] );


//: methods
private:
	inline float      norm(){ return _s*_s + _V.dot(_V); }
public:
	//:normalize it to a unit quaternion = 1.
	inline void       unit()      { float n = this->norm(); _s = _s / n; _V = _V / n; }
	inline void       normalize() { this->unit(); }

	//:only return the unit quaternion = 1
	inline Quatf    getUnit(){ return( *this) /( this->norm());}

	//:make it a conjugate.
	inline void       conjugate(){ _V[0] = -_V[0]; _V[1] = -_V[1]; _V[2] = -_V[2]; }

	//:only return the conjugate.
	inline Quatf    getConjugate(){ return Quatf(_s, -_V); }

	//:make it inversed
	inline void       inverse()   { this->conjugate(); this->unit(); }
	//:only return the inverse.
	inline Quatf    getInverse(){ return this->getConjugate().getUnit(); }

	//:return the magnitude( or absolute)
	inline float      mag()       { return sqrtf(this->norm()); }

	//:takes an angle( twist) in radians and a vector( vec).
	Quatf           setRotation( const float& rad, const Vec3<float>& vec );
	//:returns twist and dir in the form of rad, vec
	void              getRotation( float& rad, Vec3<float>& vec ) const;

	Quatf           makeRoll(const float& s);
	Quatf	          makePitch(const float& t);
	Quatf	          makeYaw(const float& p);
	
	
	//////////////////////////////////////////////////////////
	//:  Linear interpolation between two quaternion positions
	//
	//  Arguments: SlQuat( first and second quaternion), 
	//             float ( interpolation parameter [0..1]), 
	//             SlQuat( resulting quaternion, inbetween)
	//  Comments:  Fast but not nearly as smooth as Slerp
	//////////////////////////////////////////////////////////
	static void Lerp(const Quatf& from, const Quatf& to, 
	                 const float& t, Quatf& res);
	

	//////////////////////////////////////////////////////////
	//: Spherical Linear Interpolation Between two Quaternions
	// Arguments:	Two Quaternions, 
	//              blend factor, 
	//              result quaternion
	// Notes:	Adjust the DELTA constant to play with 
	//              the LERP vs. SLERP level
	//////////////////////////////////////////////////////////
	static void Slerp(const Quatf& q1, 
				const Quatf& q2,
				const float& slerp, 
				Quatf& result, 
				const float& delta = 0.0001f );
				
	//////////////////////////////////////////////////////////
	//: Spherical Linear Interpolation Between two Quaternions
	// Arguments:	Two Quaternions, 
	//              blend factor, 
	//              result quaternion,
	//              delta - lerp vs slerp level
	// Note: This one will work nicely across angle boundarys
	//       such as( 170, 0,1,0) and( 160, 0,-1,0) which are 
	//       actually only 30 degrees apart.
	//////////////////////////////////////////////////////////
	static void Slerp2(const Quatf& q1, 
				 const Quatf& q2,
				 const float& slerp, 
				 Quatf& result, 
				 const float& delta = 0.0001f );
	
	/////////////////////////////////////////////////////////
	//: Convert a set of Euler angles to a Quaternion
	// Arguments:	A rotation set of 3 angles, a quaternion to set
	// Discussion:  As the order of rotations is important I am
	//		using the Quantum Mechanics convention of( p,y,r)
	//		a Yaw-Pitch-Roll( y,p,r) system would have to be
	//		adjusted.  It is more efficient this way though.
	//////////////////////////////////////////////////////////
	void setPyrEuler(const float& p, 
		         const float& y, 
			 const float& r);
	
	//////////////////////////////////////////////////////////
	//: Convert a set of Euler angles( Y*P*R) to a Quaternion
	// Arguments:	A rotation set of 3 angles, a quaternion to set
	// Discussion:  a Yaw-Pitch-Roll( y,p,r) system.
	//////////////////////////////////////////////////////////
	void setYprEuler(const float& y, 
			 const float& p, 
			 const float& r);
	
	
	//////////////////////////////////////////////////////////
	//: Convert a set of Euler angles to a Quaternion
	// Arguments:	A rotation set of 3 angles, a quaternion to set
	// Discussion:  This is a second variation.  It creates a
	//		Series of quaternions and multiplies them together
	//		It would be easier to extend this for other rotation orders
	//////////////////////////////////////////////////////////
	void setPyrEuler2(const float& p, 
			  const float& y, 
			  const float& r);
	
//: Operators
public:
	inline float& 	      operator[]( int i )       { if(i==0) return _s; else return _V[i-1]; }
	inline const float&  operator[]( int i ) const { if(i==0) return _s; else return _V[i-1]; }
	
	inline float* getValue(float d[4]){ d[0] = _s; d[1] = _V[0]; d[2] = _V[1]; d[3] = _V[2]; return d; }	
	inline float* getwtkValue(float d[4]){ d[3] = _s; d[0] = _V[0]; d[1] = -_V[1]; d[2] = -_V[2]; return d; }	
	inline float* getpfValue(float d[4]){ return d; }
	
	//quat = quat
	Quatf& 		operator=( const Quatf& q );
    
	//quat =  float[4] quaternion
	Quatf& 		operator=( const float quat[4] );
	
	//quat = 'pure' quat
	Quatf& 		operator=( const Vec3<float>& vec );
	//////////////////////////////////////////////////////////////


	Quatf		operator*( const Quatf& q2 );
	Vec3<float>	operator*( const Vec3<float>& vec );
	Quatf		operator/( const Quatf& q2 );

	Quatf	        operator*( float val );
	friend Quatf 	operator*( float val, const Quatf& q );
	Quatf		operator/( float val );

	Quatf	operator+( const Quatf& q2 );
	Quatf	operator-( const Quatf& q2 );
	
	int		operator==( const Quatf& q2 );
	int		operator!=( const Quatf& q2 );

	friend ostream&  operator<<( ostream& out, Quatf q );
	friend istream&  operator>>( istream& in, Quatf& q );

//: Matrix dependent
public:
	//: quat = matrix
	Quatf& operator=( const Matrix4f& m );
	
//: Data members
public:
	float _s;
	Vec3<float> _V;
};
#endif

