
/////////////////////////////////////////////////////////////////
///////////////// <auto-copyright END do not edit this line> ///////////////////
#ifndef _HYPERBALL_H_
#define _HYPERBALL_H_

#include <Vec3.h>

#include "SphereBall.h"

class HyperBall : public SphereBall
{
public:
    // Constructor
    HyperBall();
    virtual ~HyperBall() {;}
    
    // places the center of rotation at the center of the screen
    // -- clears the x and y translational component
    virtual void	center();
    
    // set the center of rotation in model-space coords.
    // i.e. set it to a point on your model, and then 
    //      hyperball will rotate about that point
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
    
protected:
    virtual void _recalcMatrix();
    Matrix4f _rotation;
    Vec3<float> _translation;
    Matrix4f _offset;
    Vec3<float>  _centerOfRotation;
    
    // Size as a fraction of window size (default = 0.8).
    float _size;
};

inline void HyperBall::_recalcMatrix()
{
    // specify the center of rotation
    Matrix4f cor;
    cor.makeTranslation( _centerOfRotation );
    
    // specify the translation of the scene
    Matrix4f trans;
    trans.makeTranslation( -_translation );
    
    // set the camera (lookat) matrix
    _lookAtMatrix = trans * _rotation * cor;
}

inline void HyperBall::setCenterOfRotation(const int& x, const int& y, 
			    const int& z)
{
    _centerOfRotation.set(static_cast<float>(x), 
							static_cast<float>(y), 
							static_cast<float>(z));
    _recalcMatrix();
}
    
inline void HyperBall::setCenterOfRotation(const Vec3<float>& cor)
{
    _centerOfRotation = cor;
    _recalcMatrix();
}

//set rotation
inline void HyperBall::setRotation(const Matrix4f& rotation)
{
    _rotation.makeRotation( rotation );
    _recalcMatrix();
}

// set the translation of the scene
inline void HyperBall::setRotation(const float& rad, const float& x, const float& y, 
			    const float& z)
{
    _rotation.makeRotation( rad, x, y, z );
    _recalcMatrix();
}
			    
// set the translation of the scene
inline void HyperBall::setTranslation(const float& x, const float& y, 
			const float& z)
{
    _translation.set(x, y, z);
    _recalcMatrix();
}

// set the translation of the scene
inline void HyperBall::setTranslation(const Vec3<float>& trans)
{
    _translation = trans;
    _recalcMatrix();
}



#endif
