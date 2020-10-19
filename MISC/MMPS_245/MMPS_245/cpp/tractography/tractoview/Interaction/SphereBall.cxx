
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#include "Defines.h"
#include "Vec3.h"
#include "Matrix4f.h"
#include "Quatf.h"

#include "SphereBall.h"

//NODE_IMPLEMENT(SphereBall, TrackBall)

SphereBall::SphereBall() : TrackBall()
{ 
    _size = 0.65f;
}

/////////////////////////////////////////////////////////////////////////////
//
//  Simulate a track-ball.  Project the points onto the virtual
//  trackball, then figure out the axis of rotation, which is the cross
//  product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
//
//  Note:  This is a deformed trackball-- is a trackball in the center,
//  but is deformed into a hyperbolic sheet of rotation away from the
//  center.  This particular function was chosen after trying out
//  several variations.
//
//  It is assumed that the arguments to this routine are in the range
//  (-1.0 ... 1.0)
//
/////////////////////////////////////////////////////////////////////////////

void SphereBall::_getRotation ( const float& prevX, const float& prevY,
				const float& currentX, const float& currentY,
				Quatf& quat )
{
    Vec3<float> axis;		// Axis of rotation.
    float angle;		// How much to rotate about axis.
    Vec3<float> p1, p2, d;
    float t;

    if ( prevX == currentX && prevY == currentY ) // Zero rotation.
    {
        quat.setRotation ( 0.0f, Vec3<float>(0, 1, 0) );
        return;
    }

    // First, figure out z-coordinates for projection of P1 and P2 to
    // deformed sphere
    p1.set( prevX, prevY, _projectToSphere( _size, prevX, prevY ) );
    p2.set( currentX, currentY, _projectToSphere( _size, currentX, currentY ) );

    
    // Find the angle of rotation about the axis
    d = p1 - p2;
    t = d.length() / ( 2.0f * _size );

    // Avoid problems with out-of-control values...
    //if ( t >  1.0f ) t =  1.0;
    //if ( t < -1.0f ) t = -1.0;
    
    angle = 2.0f * asinf ( t );
    
    // Find the axis of rotation: cross product of p1 and p2
    axis = p1.cross ( p2 );
    axis.normalize();

    quat.setRotation(angle, axis);
}

/////////////////////////////////////////////////////////////////////////////
//
// Project an x,y pair onto a sphere of radius r or a hyperbolic sheet
// if we are away from the center of the sphere.
// returns the Z distance.
/////////////////////////////////////////////////////////////////////////////

float SphereBall::_projectToSphere ( float r, float x, float y )
{
    float d, t, z;

    d = sqrtf ( x * x + y * y );

    if ( d < r * INV_SQRT_TWO_F ) // Inside sphere.
	{
        z = sqrtf ( r * r - d * d );
    } 
	
	else // On hyperbola.
	{
        t = r * INV_SQRT_TWO_F;
        z = t * t / d;
    }

    return z;
}
