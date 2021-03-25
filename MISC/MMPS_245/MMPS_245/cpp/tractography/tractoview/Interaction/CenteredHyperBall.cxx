
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#include "Defines.h"
#include "Vec3.h"
#include "Matrix4f.h"
#include "Quatf.h"

#include "CenteredHyperBall.h"

//NODE_IMPLEMENT(CenteredHyperBall, SphereBall)

CenteredHyperBall::CenteredHyperBall() : SphereBall(), 
			_rotation( Matrix4f::identity() ), 
			_translation( Matrix4f::identity() )
		    
{ 
    _size = 0.65f;
}

void CenteredHyperBall::applyYtrans(const int& yPosition, 
			    const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    
    //define up and xaxis vectors
    Vec3<float> up(0.0f, 1.0f, 0.0f);
    
    //strafe vertically
    //Matrix4f temp( Matrix4f::identity() );
    //temp.translate((up * (-offset*0.1)));
    //_rotation.multLeft(temp);
    _rotation.translate((up * (-offset*0.1f)));
    
    _recalcMatrix();
}

void CenteredHyperBall::applyXtrans(const int& xPosition, 
			    const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    
   
    // define local X axis
    Vec3<float> xaxis(-1, 0, 0);
    
    // straef sideways.
    Matrix4f temp( Matrix4f::identity() );
    temp.translate((xaxis * (-offset*0.1f)));
    _rotation.multLeft(temp);
    
    _recalcMatrix();
}

void CenteredHyperBall::applyZtrans(const int& zPosition, 
			    const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    
    // define forward axis (along Z)
    Vec3<float> forward(0.0f, 0.0f, -1.0f);	
    
    // translate forward.
    _translation.translate((forward * (-offset*0.1f)));

    _recalcMatrix();
}

void CenteredHyperBall::applyXYtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    int xOffset = xPosition - oldXposition;
    int yOffset = yPosition - oldYposition;
    
    //define up and xaxis vectors
    Vec3<float> up(0.0f, 1.0f, 0.0f);
    Vec3<float> xaxis(-1.0f, 0.0f, 0.0f);
    
    //strafe
    Matrix4f temp( Matrix4f::identity() );
    temp.translate((xaxis * (-xOffset*0.1f)));
    temp.translate((up * (-yOffset*0.1f)));
    _rotation.multLeft(temp);
    
    _recalcMatrix();
}

void CenteredHyperBall::applyXZtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    int xOffset = xPosition - oldXposition;
    int yOffset = yPosition - oldYposition;
    
    // define local X axis
    Vec3<float> xaxis(-1.0f, 0.0f, 0.0f);
    
    // straef sideways.
    Matrix4f temp( Matrix4f::identity() );
    temp.translate((xaxis * (-xOffset*0.1f)));
    _rotation.multLeft(temp);
    
    // define forward axis (along Z)
    Vec3<float> forward(0.0f, 0.0f, -1.0f);	
    
    // translate forward.
    _translation.translate((forward * (-yOffset*0.1f)));
    
    _recalcMatrix();
}

void CenteredHyperBall::applyXrot(const int& xPosition, 
			  const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 1, 0, 0);

    _recalcMatrix();
}

void CenteredHyperBall::applyYrot(const int& yPosition, 
			  const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 0, 1, 0);

    _recalcMatrix();
}

void CenteredHyperBall::applyZrot(const int& zPosition, 
			  const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 0, 0, 1);

    _recalcMatrix();
}

void CenteredHyperBall::applyGeneralRot(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    //CenteredHyperBall
    Quatf quat;
    int left, right, top, bottom;
    this->getViewport( TrackBall::top, top );
    this->getViewport( TrackBall::bottom, bottom );
    this->getViewport( TrackBall::right, right );
    this->getViewport( TrackBall::left, left );
    SphereBall::_getRotation(
	    2.0f*( oldXposition/((float)right  ) - 0.5f),
	   -2.0f*( oldYposition/((float)bottom ) - 0.5f), 
	    2.0f*(    xPosition/((float)right  ) - 0.5f),
	   -2.0f*(    yPosition/((float)bottom ) - 0.5f),
	   quat ); 

    //rotate the new lookat by 'quat'.
    _offset = quat;
    _rotation.multLeft( _offset );
    
    _recalcMatrix();
}

/*
void CenteredHyperBall::handleMouseEvent( const Mouse& mouse )
{
    Vec2<float>	    offset;
    int		    x, y;

    mouse.getPosition( x, y );
    
    ///////////////////////////////////////////////////
    // Left button
    
    if( mouse.leftEdgeState() == Mouse::ONDOWN )
    {
    	_lPosition = _lOldPosition = Vec2<float>( x, y );
    }
    
    if( mouse.leftBinaryState() == Mouse::ON )
    {
	//get the current mouse position
	_lPosition = Vec2<float>( x, y );
	
	//compute the offset that the mouse has moved.
	offset = _lPosition - _lOldPosition;
	
	//save the old position
	_lOldPosition = _lPosition;

	// define forward axis (along Z)
	Vec3<float> forward(0, 0, -1);
	
	// define local X axis
	Vec3<float> xaxis(-1, 0, 0);
	
	//Apply the translation
	Matrix4f temp( Matrix4f::identity() );
	temp.translate((xaxis * (-offset[0]*0.1)));
	_rotation.multLeft(temp);
	
	_translation.translate((forward * (-offset[1]*0.1)));
    }
    
    ///////////////////////////////////////////////////
    // Right button
    if( mouse.rightEdgeState() == Mouse::ONDOWN )
    {
    	_rPosition = _rOldPosition = Vec2<float>( x, y );
    }
    
    if( mouse.rightBinaryState() == Mouse::ON )
    {
	//get the current mouse position
	_rPosition = Vec2<float>( x, y );

	//CenteredHyperBall
	Quatf quat;
	int left, right, top, bottom;
	this->getViewport(left, top, right, bottom);
	SphereBall::_getRotation(  
	     2.0f*(_rOldPosition[0] / ((float)right) - 0.5f),
	    -2.0f*(_rOldPosition[1] / ((float)bottom) - 0.5f), 
	     2.0f*(_rPosition[0] / ((float)right) - 0.5f),
	    -2.0f*(_rPosition[1] / ((float)bottom) - 0.5f),
	     quat ); 
	
	//rotate the new lookat by 'quat'.
	_offset = quat;
	_rotation.multLeft( _offset );
	
	//save the old position
	_rOldPosition = _rPosition;
    }
    
    ///////////////////////////////////////////////////
    // Middle button
    
    if( mouse.middleEdgeState() == Mouse::ONDOWN )
    {
    	_mPosition = _mOldPosition = Vec2<float>( x, y );
    }
    
    if( mouse.middleBinaryState() == Mouse::ON )
    {
	//get the current mouse position
	_mPosition = Vec2<float>( x, y );
	
	//compute the offset that the mouse has moved.
	offset = _mPosition - _mOldPosition;
	
	//save the old position
	_mOldPosition = _mPosition;
	
	//define up and xaxis vectors
	Vec3<float> up(0, 1, 0);
	Vec3<float> xaxis(-1, 0, 0);
	
	//strafe
	Matrix4f temp( Matrix4f::identity() );
	temp.translate((xaxis * (-offset[0]*0.1)));
	temp.translate((up * (-offset[1]*0.1)));
	_rotation.multLeft(temp);
    }  
    
    Matrix4f cor( Matrix4f::identity() );
    cor.translate(_centerOfRotation);
    
    // set the camera (lookat) matrix
    _lookAtMatrix =  _translation * cor * _rotation;
}
*/

// places the center of rotation at the center of the screen
void CenteredHyperBall::center()
{
    _centerOfRotation.set(0, 0, 0);
}
