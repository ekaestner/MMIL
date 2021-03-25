
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#ifndef SPHEREBALL_INCLUDED
#define SPHEREBALL_INCLUDED

class SlVec4f;
class SlVec2f;
template <class dataType> class Vec3;
class SlMatrix4f;

#include "TrackBall.h"
class Quatf;
class SphereBall : public TrackBall
{
	//NODE_DECLARE(SphereBall, TrackBall)
public:
    // Constructor
    SphereBall();
    virtual ~SphereBall() {;}
    
    // places the center of rotation at the center of the screen
    // -- clears the x and y translational component
    virtual void	center() = 0;
    
    // set the center of rotation in model-space coords.
    // i.e. set it to a point on your model, and then 
    //      hyperball will rotate about that point
    virtual void	setCenterOfRotation(const int& x, const int& y, 
			    const int& z) = 0;
    
    // set the center of rotation in world coords.
    virtual void	setCenterOfRotation(const Vec3<float>& cor) = 0;
    
protected:

    // Pass the x and y coordinates of the last and current positions of
    // the mouse, scaled so they are from (-1.0 ... 1.0).
    // The resulting rotation is returned as a quaternion rotation in the
    // first paramater.
    
    void _getRotation ( const float& prevX, const float& prevY,
			    const float& currentX, const float& currentY,
			    Quatf& quat );
			    
    // Project an x,y pair onto a sphere of radius r or a hyperbolic sheet
    // if we are away from the center of the sphere.

    static float _projectToSphere ( float r, float x, float y );

			
protected:
    virtual void _recalcMatrix() = 0;
    
    // Size as a fraction of window size (default = 0.8).
    float _size;
};




#endif
