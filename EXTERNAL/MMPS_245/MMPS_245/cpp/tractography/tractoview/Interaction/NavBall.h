
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#ifndef _NAVBALL_H_
#define _NAVBALL_H_

#include "TrackBall.h"
class NavBall : public TrackBall
{
	//NODE_DECLARE(NavBall, TrackBall)
public:
    NavBall();
    virtual ~NavBall() {;}
    
    // get the eye position
    virtual const Vec3<float>& getTranslation();
    
    // set translation
    virtual void setTranslation(const Vec3<float>& trans);    
    
    // set translation
    virtual void setTranslation(const float& x, const float& y, const float& z);    
    
    //set rotation
    virtual void setRotation(const Matrix4f& trans);
    
    //set rotation
    virtual void setRotation(const float& rad, const float& x, const float& y, const float& z);
    
    // get the up direction reletive to your eye
    const Vec3<float>& getUp();
    void setUp( const Vec3<float>& up );
    
    // get the direction you're looking
    const Vec3<float>& getDirection();
    void setDirection( const Vec3<float>& dir );
    
    // get the left direction in relation to your eye
    Vec3<float> getX();
    
    // modify the matrix such that x mouse offset will
    // translate trackball in the Y axis.
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyXtrans(const int& xPosition, 
			    const int& oldxPosition);
    
    // modify the matrix such that y mouse offset will
    // translate trackball in the Y axis.
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyYtrans(const int& yPosition, 
			    const int& oldyPosition);
    
    // modify the matrix such that y mouse offset will
    // translate trackball in the Z axis.
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyZtrans(const int& yPosition, 
			    const int& oldyPosition);
    
    // modify the matrix such that x,y mouse offset will
    // translate trackball in the XY plane.
    // give - the x,y mouse current and previous positions
    // result - matrix is offset by the deltas
    virtual void applyXYtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition);
    
    // modify the matrix such that x,y mouse offset will
    // translate trackball in the XZ plane.
    // give - the x,y mouse current and previous positions
    // result - matrix is offset by the deltas
    virtual void applyXZtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition);
    
    // modify the matrix such that y mouse offset will
    // rotate trackball about the X axis
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyXrot(const int& xPosition, 
			  const int& oldxPosition);
    
    // modify the matrix such that x mouse offset will
    // rotate trackball about the Y axis
    // give - the x mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyYrot(const int& yPosition, 
			  const int& oldyPosition);
    
    // modify the matrix such that x mouse offset will
    // rotate trackball about the Z axis
    // give - the x mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyZrot(const int& xPosition, 
			  const int& oldxPosition);
    
    // modify the matrix such that x,y mouse offsets will
    // rotate by the derived-class-defined method.
    // give - the x,y mouse current and previous positions
    // result - matrix is offset by the deltas
    virtual void applyGeneralRot(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition);

protected:
    void _recalcMatrix();
    Vec3<float> _eye;
    Vec3<float> _up;
    Vec3<float> _direction;
};

inline void NavBall::_recalcMatrix()
{
    // update the lookat matrix...
    Vec3<float> center = _eye + _direction;
    _lookAtMatrix.makeLookAt( _eye, center, _up );
}

// get the eye position
inline const Vec3<float>& NavBall::getTranslation()
{
    return _eye;
}

// Set the eye position
inline void NavBall::setTranslation( const float& x, const float& y, const float& z )
{
    _eye.set( x, y, z );
    _recalcMatrix(); //TODO: make this more efficient.
}    

//set rotation
inline void NavBall::setRotation(const Matrix4f& rotation)
{
    Vec3<float> forward(0, 0, -1);
    Vec3<float> up(0, 1, 0);
    
    _direction = rotation * forward;
    _up = rotation * up;
    
    _recalcMatrix();
}

// set the translation of the scene
inline void NavBall::setRotation(const float& rad, const float& x, const float& y, 
			    const float& z)
{
    Matrix4f rotation;
    rotation.makeRotation( -rad, x, y, z );
    Vec3<float> forward(0, 0, -1);
    Vec3<float> up(0, 1, 0);
    
    _direction = rotation * forward;
    _up = rotation * up;
    //cout<<_up<<" "<<_direction<<"\n"<<flush;
    
    _recalcMatrix();
}

// Set the eye position
inline void NavBall::setTranslation( const Vec3<float>& eye )
{
    _eye = eye;
    _recalcMatrix(); //TODO: make this more efficient.
}

// get the up direction reletive to your eye
inline const Vec3<float>& NavBall::getUp()
{
    return _direction;
}

inline void NavBall::setUp( const Vec3<float>& up )
{
    _up = up;
    _recalcMatrix();
}

// get the direction you're looking
inline const Vec3<float>& NavBall::getDirection()
{
    return _direction;
}

inline void NavBall::setDirection( const Vec3<float>& dir )
{
    _direction = dir;
    _recalcMatrix();
}

// get the -X cameraspace axis in world space...
inline Vec3<float> NavBall::getX()
{
    return _up.cross(_direction);
}

#endif
