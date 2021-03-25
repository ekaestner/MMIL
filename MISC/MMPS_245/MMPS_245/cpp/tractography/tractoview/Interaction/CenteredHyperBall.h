
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#ifndef _CENTEREDHYPERBALL_H_
#define _CENTEREDHYPERBALL_H_

#include "Vec4.h"
#include "Vec2.h"
#include "Matrix4f.h"

#include "SphereBall.h"

class CenteredHyperBall : public SphereBall
{
public:
    // Constructor
    CenteredHyperBall();
    virtual ~CenteredHyperBall() {;}
    
    // places the center of rotation at the center of the screen
    // -- clears the x and y translational component
    virtual void	center();
    
    // set the global center of rotation in trackball-space coords.
    // set this equal to trackball location (0,0,0)
    //      for an EZ intuitive interface
    // i.e. the trackball is in the center of the screen,
    //      and geometry can only rotate about the trackball's
    //      center.
    // result - geometry will orbit about the point 
    //          regardless of it's translation
    virtual void	setCenterOfRotation(const int& x, const int& y, 
			    const int& z);
    
    // set the center of rotation in world coords.
    virtual void	setCenterOfRotation(const Vec3<float>& cor);
    
    // set translation
    virtual void setTranslation(const Vec3<float>& trans);    
    
    // set translation
    virtual void setTranslation(const float& x, const float& y, const float& z);    
    
    //set rotation
    virtual void setRotation(const Matrix4f& trans);
    
    //set rotation
    virtual void setRotation(const float& rad, const float& x, const float& y, const float& z);
    
    
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
			    
private:
    virtual void _recalcMatrix();
    Matrix4f _rotation;
    Matrix4f _translation;
    Matrix4f _offset;
    Vec3<float>  _centerOfRotation;
    
    // Size as a fraction of window size (default = 0.8).
    float _size;
};

inline void CenteredHyperBall::_recalcMatrix()
{
    // specify the center of rotation
    Matrix4f cor( Matrix4f::identity() );
    cor.translate(_centerOfRotation);
    
    // specify the translation of the scene
    Matrix4f trans( Matrix4f::identity() );
    trans.copyTranslation(-_translation);
    
    // set the camera (lookat) matrix
    _lookAtMatrix = trans * _rotation * cor;
}

//set rotation
inline void CenteredHyperBall::setRotation(const Matrix4f& rotation)
{
    _rotation.setRotation( rotation );
    _recalcMatrix();
}

// set the translation of the scene
inline void CenteredHyperBall::setRotation(const float& rad, const float& x, const float& y, 
			    const float& z)
{
    _rotation.makeRotation( rad, x, y, z );
    _recalcMatrix();
}
			    
// set the translation of the scene
inline void CenteredHyperBall::setTranslation(const float& x, const float& y, 
			const float& z)
{
    _translation.translate(-x, -y, -z);
    _recalcMatrix();
}

// set the translation of the scene
inline void CenteredHyperBall::setTranslation(const Vec3<float>& trans)
{
    _translation.translate(-trans);
    _recalcMatrix();
}

inline void CenteredHyperBall::setCenterOfRotation(const int& x, const int& y, 
			    const int& z)
{
    _centerOfRotation.set((float)x, (float)y, (float)z);
    _recalcMatrix();
}
    
inline void CenteredHyperBall::setCenterOfRotation(const Vec3<float>& cor)
{
    _centerOfRotation = cor;
    _recalcMatrix();
}

#endif
