
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#include "Matrix4f.h"

#include "NavBall.h"  // My Header

//NODE_IMPLEMENT(NavBall, TrackBall)


NavBall::NavBall() : TrackBall(), 
					_eye(0, 0, 0), 
					_up(0, 1, 0), 
					_direction(0, 0, -1)
{
}

void NavBall::applyYtrans(const int& yPosition, 
			    const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    
    _eye += _up * (-offset*0.1f);

    _recalcMatrix();
}

void NavBall::applyXtrans(const int& xPosition, 
			    const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    
    //(slide horizontally in user's local coordinate system)
    Vec3<float> localXaxis = _up.cross(_direction);
    _eye += localXaxis * (-offset*0.1f);
    
    _recalcMatrix();
}

void NavBall::applyZtrans(const int& zPosition, 
			    const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    
    //go forward with respect to your view.
    _eye += _direction * (-offset*0.1f);

    _recalcMatrix();
}

void NavBall::applyXYtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    int xOffset = xPosition - oldXposition;
    int yOffset = yPosition - oldYposition;
    
    //strafe
    //(slide up or down in user's local coordinate system)
    Vec3<float> localXaxis = _up.cross(_direction);
    _eye += localXaxis * (-xOffset*0.1f);
    _eye += _up * (-yOffset*0.1f);
    
    _recalcMatrix();
}

//strafe in the XZ plane
void NavBall::applyXZtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    int xOffset = xPosition - oldXposition;
    int yOffset = yPosition - oldYposition;
    
    //(slide horizontally in user's local coordinate system)
    Vec3<float> localXaxis = _up.cross(_direction);
    _eye += localXaxis * (-xOffset*0.1f);
    
    //go forward with respect to your view.
    _eye += _direction * (-yOffset*0.1f);
   
    _recalcMatrix();
}

void NavBall::applyXrot(const int& xPosition, 
			  const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    
    //Pitch
    if( ((int)offset) != 0 )
    {
	float rad = (offset*0.1f) * TO_RAD_F;
	Vec3<float> localXaxis = _up.cross(_direction);
	Matrix4f rotation;
	rotation.makeRotation( rad, localXaxis );
	_up = rotation * _up;
	_direction = rotation * _direction;
    }
    
    _recalcMatrix();
}

void NavBall::applyYrot(const int& yPosition, 
			  const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    
    //offset the view direction
    if( ((int)offset) != 0 )
    {   
	float rad = (-offset*0.1f) * TO_RAD_F;
        Matrix4f rotation;
	rotation.makeRotation( rad, _up );
	_direction = rotation * _direction;
    }
    
    _recalcMatrix();
}

void NavBall::applyZrot(const int& zPosition, 
			  const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    //Roll
    if( ((int)offset) != 0 )
    {
	float rad = (-offset*0.1f) * TO_RAD_F;
	Vec3<float> axis = _direction;
	axis.normalize();
        Matrix4f rotation;
	rotation.makeRotation( rad, axis );
	_up = rotation * _up;
    }
    
    _recalcMatrix();
}

void NavBall::applyGeneralRot(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    float xOffset = (float) (xPosition - oldXposition);
    float yOffset = (float) (yPosition - oldYposition);
    
    //Roll
    Matrix4f rotation;
    if( ((int)xOffset) != 0 )
    {
	float rad = (-xOffset*0.1f) * TO_RAD_F;
	Vec3<float> axis = _direction;
	axis.normalize();
	rotation.makeRotation( rad, axis );
	_up = rotation * _up;
    }
    
    //Pitch
    if( ((int)yOffset) != 0 )
    {
	float rad = (yOffset*0.1f) * TO_RAD_F;
	Vec3<float> localXaxis = _up.cross(_direction);
	rotation.makeRotation( rad, localXaxis );
	_up = rotation * _up;
	_direction = rotation * _direction;
    }
    
    _recalcMatrix();
}

/*
void NavBall::handleMouseEvent( const Mouse& mouse )
{
    SlVec2f	    offset;
    int		    x, y;
 
    mouse.getPosition( x, y );
    
    ///////////////////////////////////////////////////
    // Left button
    
    if( mouse.leftEdgeState() == Mouse::ONDOWN )
    {
    	_lPosition = _lOldPosition = SlVec2f( x, y );
    }
    
    if( mouse.leftBinaryState() == Mouse::ON )
    {
	//get the current mouse position
	_lPosition = SlVec2f( x, y );
	
	//compute the offset that the mouse has moved.
	offset = _lPosition - _lOldPosition;
	
	//save the old position
	_lOldPosition = _lPosition;
	
	//offset the view direction
	Matrix4f rotation;
	if( ((int)offset[0]) != 0 )
	{   
	    float rad = (-offset[0]*0.1) * TO_RAD_F;
	    rotation.makeRotation( rad, _up );
	    _direction = rotation * _direction;
	}
	_eye += _direction * (-offset[1]*0.1);
    }
    
    ///////////////////////////////////////////////////
    // Right button
    
    if( mouse.rightEdgeState() == Mouse::ONDOWN )
    {
    	_rPosition = _rOldPosition = SlVec2f( x, y );
    }
    
    if( mouse.rightBinaryState() == Mouse::ON )
    {
	//get the current mouse position
	_rPosition = SlVec2f( x, y );
	
	//compute the offset that the mouse has moved.
	offset = _rPosition - _rOldPosition;
	
	//save the old position
	_rOldPosition = _rPosition;
	
	//Roll
	Matrix4f rotation;
	if( ((int)offset[0]) != 0 )
	{
	    float rad = (-offset[0]*0.1) * TO_RAD_F;
	    Vec3<float> axis = _direction;
	    axis.normalize();
	    rotation.makeRotation( rad, axis );
	    _up = rotation * _up;
	}
	
	//Pitch
	if( ((int)offset[1]) != 0 )
	{
	    float rad = (offset[1]*0.1) * TO_RAD_F;
	    Vec3<float> localXaxis = _up.cross(_direction);
	    rotation.makeRotation( rad, localXaxis );
	    _up = rotation * _up;
	    _direction = rotation * _direction;
	}
    }
    
    ///////////////////////////////////////////////////
    // Middle button
    
    if( mouse.middleEdgeState() == Mouse::ONDOWN )
    {
    	_mPosition = _mOldPosition = SlVec2f( x, y );
    }
    
    if( mouse.middleBinaryState() == Mouse::ON )
    {
	//get the current mouse position
	_mPosition = SlVec2f( x, y );
	
	//compute the offset that the mouse has moved.
	offset = _mPosition - _mOldPosition;
	
	//save the old position
	_mOldPosition = _mPosition;
	
	//strafe
	//(slide up or down in user's local coordinate system)
	Vec3<float> localXaxis = _up.cross(_direction);
	_eye += localXaxis * (-offset[0]*0.1);
	_eye += _up * (-offset[1]*0.1);
    }  
    
    // update the lookat matrix...
    _recalcMatrix();
}
*/
