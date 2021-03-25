
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
#include "Defines.h"
#include "Matrix4f.h"
#include "Quatf.h"

#include "HyperBall.h"

HyperBall::HyperBall() : SphereBall(), 
			_rotation(Matrix4f::identity()), 
			_translation(0, 0, 0), 
			_centerOfRotation(0, 0, 0)

{
	_size = 0.65f;
}

void HyperBall::applyYtrans(const int& yPosition, 
			    const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    
    //define up and xaxis vectors
    Vec3<float> up(0, 1, 0);

    //strafe vertically
    //(slide up in user's local coordinate system)
    _translation += ( up * (offset * _sensativity) );

   cout << "y trans" << endl; 

    _recalcMatrix();
}

void HyperBall::applyXtrans(const int& xPosition, 
			    const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    
    //define up and xaxis vectors
    Vec3<float> xaxis(-1, 0, 0);

    //strafe
    //(slide side to side in user's local coordinate system)
    _translation += ( xaxis * (offset * _sensativity) );

   cout << "x trans" << endl; 
    _recalcMatrix();
}

void HyperBall::applyZtrans(const int& zPosition, 
			    const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    
    // define forward axis (along Z)
    Vec3<float> forward(0, 0, -1);

    //Apply the translation
    _translation += ( forward * (offset * _sensativity) );

   
   //Cout << "z trans" << endl; 
   cout << "new pos: " << zPosition << ", old pos:" << oldzPosition<< endl; 

    _recalcMatrix();
}

void HyperBall::applyXYtrans(const int& xPosition, 
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
    
//   cout << "xy trans" << endl; 
    _recalcMatrix();
}

void HyperBall::applyXZtrans(const int& xPosition, 
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
    
//   cout << "xz trans" << endl; 
    _recalcMatrix();
}

void HyperBall::applyXrot(const int& xPosition, 
			  const int& oldxPosition)
{
    int offset = xPosition - oldxPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 1, 0, 0 );

//   cout << "x trans" << endl; 
    _recalcMatrix();
}

void HyperBall::applyYrot(const int& yPosition, 
			  const int& oldyPosition)
{
    int offset = yPosition - oldyPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 0, 1, 0);
//   cout << "y trans" << endl; 

    _recalcMatrix();
}

void HyperBall::applyZrot(const int& zPosition, 
			  const int& oldzPosition)
{
    int offset = zPosition - oldzPosition;
    _rotation.rotate( offset * _sensativity * TO_RAD_F, 0, 0, 1);
//   cout << "z trans" << endl; 

    _recalcMatrix();
}

void HyperBall::applyGeneralRot(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition)
{
    //HyperBall
    Quatf quat;
    int left, right, top, bottom;
    this->getViewport( TrackBall::left, left );
    this->getViewport( TrackBall::top, top );
    this->getViewport( TrackBall::bottom, bottom );
    this->getViewport( TrackBall::right, right );
    
    SphereBall::_getRotation(
            2.0f*( oldXposition/((float)right  ) - 0.5f),
            -2.0f*( oldYposition/((float)bottom ) - 0.5f), 
            2.0f*(    xPosition/((float)right  ) - 0.5f),
            -2.0f*(    yPosition/((float)bottom ) - 0.5f),
            quat );

    //rotate the new lookat by 'quat'.
    _offset = quat;
    _rotation.multLeft( _offset );
    
   //cout << "general rot " << endl; 
    _recalcMatrix();
}



/*
void HyperBall::handleMouseEvent( const Mouse& mouse )
{
	SlVec2f     offset;
	int         x, y;

	mouse.getPosition( x, y );

	///////////////////////////////////////////////////
	// Left button XZ trans

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

		// define forward axis (along Z)
		Vec3<float> forward(0, 0, -1);

		// define local X axis
		Vec3<float> xaxis(-1, 0, 0);

		//Apply the translation
		_translation += (forward * (offset[1]*_sensativity));
		_translation += (xaxis * (offset[0]*_sensativity));
	}

	///////////////////////////////////////////////////
	// Right button General Rotatoin
	if( mouse.rightEdgeState() == Mouse::ONDOWN )
	{
		_rPosition = _rOldPosition = SlVec2f( x, y );
	}

	if( mouse.rightBinaryState() == Mouse::ON )
	{
		//get the current mouse position
		_rPosition = SlVec2f( x, y );

		//HyperBall
		Quatf quat;
		int left, right, top, bottom;
		this->getViewport(left, top, right, bottom);
		SphereBall::_getRotation(
		    2.0f*(_rOldPosition[0]/((float)right) - 0.5f),
		    -2.0f*(_rOldPosition[1]/((float)bottom) - 0.5f), 
		    2.0f*(_rPosition[0]/((float)right) - 0.5f),
		    -2.0f*(_rPosition[1]/((float)bottom) - 0.5f),
		    quat ); 

		//rotate the new lookat by 'quat'.
		_offset = quat;
		_rotation.multLeft( _offset );

		//save the old position
		_rOldPosition = _rPosition;
	}

	///////////////////////////////////////////////////
	// Middle button XY Trans

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

		//define up and xaxis vectors
		Vec3<float> up(0, 1, 0);
		Vec3<float> xaxis(-1, 0, 0);

		//strafe
		//(slide up or down in user's local coordinate system)
		_translation += (xaxis * (offset[0]*_sensativity));
		_translation += (up * (offset[1]*_sensativity));
	}

	// specify the center of rotation
	//Matrix4f cor(IDENTITY_MATRIX4F);
	//cor.translate(_centerOfRotation);

	// specify the translation of the scene
	//Matrix4f trans(IDENTITY_MATRIX4F);
	//trans.setTranslation(-_translation);

	// set the camera (lookat) matrix
	//_lookAtMatrix = trans * _rotation * cor;

    _recalcMatrix();
}
*/


// places the center of rotation at the center of the screen
// -- clears the x and y translational component
void HyperBall::center()
{
	_translation.set(0, 0, _translation[2]);
	_recalcMatrix();
}
