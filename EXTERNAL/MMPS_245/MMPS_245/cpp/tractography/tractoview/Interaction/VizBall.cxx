
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
#include "Matrix4f.h"

#include "VizBall.h"

VizBall::VizBall() : TrackBall(), 
			_rotation( Matrix4f::identity() ), 
			_translation(0, 0, 0), 
			_centerOfRotation(0, 0, 0)
{
	_size = 0.65f;
}


void VizBall::applyXtrans(const int& xPosition, 
			    const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    
    //define up and xaxis vectors
    Vec3<float> xaxis(-1, 0, 0);

    //strafe
    //(slide side to side in user's local coordinate system)
    _translation += ( xaxis * (offset * _sensativity) );

    _recalcMatrix();
}

void VizBall::applyYtrans(const int& yPosition, 
			    const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    
    //define up and xaxis vectors
    Vec3<float> up(0, 1, 0);

    //strafe vertically
    //(slide up in user's local coordinate system)
    _translation += ( up * (offset * _sensativity) );

    _recalcMatrix();
}

void VizBall::applyZtrans(const int& zPosition, 
			    const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    
    // define forward axis (along Z)
    Vec3<float> forward(0, 0, -1);

    //Apply the translation
    _translation += ( forward * (offset * _sensativity) );

    _recalcMatrix();
}

void VizBall::applyXYtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    int xOffset = xPosition - oldXposition;
    int yOffset = yPosition - oldYposition;
    
    //define up and xaxis vectors
    Vec3<float> up(0, 1, 0);
    Vec3<float> xaxis(-1, 0, 0);

    //strafe
    //(slide up or down in user's local coordinate system)
    _translation += ( xaxis * (xOffset * _sensativity) );
    _translation += ( up * (yOffset * _sensativity) );
    
    _recalcMatrix();
}

void VizBall::applyXZtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    int xOffset = xPosition - oldXposition;
    int zOffset = yPosition - oldYposition;
    
    // define forward axis (along Z)
    Vec3<float> forward(0, 0, -1);

    // define local X axis
    Vec3<float> xaxis(-1, 0, 0);

    //Apply the translation
    _translation += ( forward * (zOffset * _sensativity) );
    _translation += ( xaxis * (xOffset * _sensativity) );
    
    _recalcMatrix();
}

void VizBall::applyXrot(const int& xPosition, 
			  const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 1, 0, 0);

    _recalcMatrix();
}

void VizBall::applyYrot(const int& yPosition, 
			  const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 0, 1, 0);

    _recalcMatrix();
}

void VizBall::applyZrot(const int& zPosition, 
			  const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 0, 0, 1);

   cout << "xz trans" << endl; 
    _recalcMatrix();
}

void VizBall::applyGeneralRot(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
	int xDelta = xPosition - oldXposition;
	int yDelta = yPosition - oldYposition;
    Vec3<float> offset( static_cast<float>( xDelta ), 
						static_cast<float>( yDelta ), 
						0.0f );
    
    //if mouse hasn't moved (length == 0), we would get bad data...
    if( offset.length() != 0)
    {
		Vec3<float> zaxis(0,0,1);
		Vec3<float> xyProjected(offset[0], -offset[1], 0);
		float length = xyProjected.length();
		float twist = -length / 600 * 360;
		
		xyProjected.normalize();
		Vec3<float> vector( xyProjected.cross(zaxis) );
		Matrix4f matrix;
		matrix.makeRotation( twist * TO_RAD_F, vector );
		_rotation.multLeft(matrix);

		_recalcMatrix();
    }
}

/*
void VizBall::handleMouseEvent( const Mouse& mouse )
{
	Vec3<float>     offset;
	int         x, y;

	mouse.getPosition( x, y );

	///////////////////////////////////////////////////
	// Left button

	if( mouse.leftEdgeState() == Mouse::ONDOWN )
	{
		_lPosition = _lOldPosition = Vec3<float>( x, y );
	}

	if( mouse.leftBinaryState() == Mouse::ON )
	{
		//get the current mouse position
		_lPosition = Vec3<float>( x, y );

		//compute the offset that the mouse has moved.
		offset = _lPosition - _lOldPosition;

		//save the old position
		_lOldPosition = _lPosition;

		// define forward axis (along Z)
		Vec3<float> forward(0, 0, -1);

		// define local X axis
		Vec3<float> xaxis(-1, 0, 0);

		//Apply the translation
		_translation += (forward * (offset[1]*0.1));
		_translation += (xaxis * (offset[0]*0.1));
	}
	
    
	///////////////////////////////////////////////////
	// Right button

	if( mouse.rightEdgeState() == Mouse::ONDOWN )
	{
		_rPosition = _rOldPosition = Vec3<float>( x, y );
	}

	if( mouse.rightBinaryState() == Mouse::ON )
	{
		//get the current mouse position
		_rPosition = Vec3<float>( x, y );

		//compute the offset that the mouse has moved.
		offset = _rPosition - _rOldPosition;

		//save the old position
		_rOldPosition = _rPosition;

		//if mouse hasn't moved, we would get bad data...
		if(offset.length() != 0)
		{
		    Vec3<float> zaxis(0,0,1);
		    Vec3<float> xyProjected(offset[0], -offset[1], 0);
		    float length = xyProjected.length();
		    float twist = -length / 600 * 360;
		    
		    xyProjected.normalize();
		    Vec3<float> vector( xyProjected.cross(zaxis) );
		    Matrix4f matrix(IDENTITY_MATRIX4F);
		    matrix.makeRotation( twist * TO_RAD_F, vector );
		    _rotation.multLeft(matrix);
		}
   	}

	///////////////////////////////////////////////////
	// Middle button

	if( mouse.middleEdgeState() == Mouse::ONDOWN )
	{
		_mPosition = _mOldPosition = Vec3<float>( x, y );
	}

	if( mouse.middleBinaryState() == Mouse::ON )
	{
		//get the current mouse position
		_mPosition = Vec3<float>( x, y );

		//compute the offset that the mouse has moved.
		offset = _mPosition - _mOldPosition;

		//save the old position
		_mOldPosition = _mPosition;

		//define up and xaxis vectors
		Vec3<float> up(0, 1, 0);
		Vec3<float> xaxis(-1, 0, 0);

		//strafe
		//(slide up or down in user's local coordinate system)
		_translation += (xaxis * (offset[0]*0.1));
		_translation += (up * (offset[1]*0.1));
	}

    // update the lookat matrix...
	_recalcMatrix();
}
*/

// places the center of rotation at the center of the screen
// -- clears the x and y translational component
void VizBall::center()
{
	_translation.set( 0, 0, _translation[2] );
}
