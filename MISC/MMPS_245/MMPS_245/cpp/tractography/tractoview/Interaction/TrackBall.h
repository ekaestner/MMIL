
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
#ifndef TRACKBALL_INCLUDED
#define TRACKBALL_INCLUDED

#include <assert.h>

//TODO: make TrackBall take an external matrix.
//      then TrackBall will act on this matrix.
//      should provide for interesting possibilities, 
//      like switching to another TrackBall in mid-application, 
//      and still stay in the same location, since your matrix hasn't changed yet.
//      otherwise you'd have to sync all the matrices somehow between all the TrackBalls..
class Matrix4f;
class TrackBall //: public Node
{
	//NODE_DECLARE(TrackBall, Node)
public:
    TrackBall();
    virtual ~TrackBall() {}
    
    // modify the matrix such that x mouse offset will
    // translate trackball in the Y axis.
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyXtrans(const int& xPosition, 
			    const int& oldxPosition) = 0;
    
    // modify the matrix such that y mouse offset will
    // translate trackball in the Y axis.
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyYtrans(const int& yPosition, 
			    const int& oldyPosition) = 0;
    
    // modify the matrix such that y mouse offset will
    // translate trackball in the Z axis.
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyZtrans(const int& yPosition, 
			    const int& oldyPosition) = 0;
    
    // modify the matrix such that x,y mouse offset will
    // translate trackball in the XY plane.
    // give - the x,y mouse current and previous positions
    // result - matrix is offset by the deltas
    virtual void applyXYtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition) = 0;
    
    // modify the matrix such that x,y mouse offset will
    // translate trackball in the XZ plane.
    // give - the x,y mouse current and previous positions
    // result - matrix is offset by the deltas
    virtual void applyXZtrans(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition) = 0;
    
    // modify the matrix such that y mouse offset will
    // rotate trackball about the X axis
    // give - the y mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyXrot(const int& xPosition, 
			  const int& oldxPosition) = 0;
    
    // modify the matrix such that x mouse offset will
    // rotate trackball about the Y axis
    // give - the x mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyYrot(const int& yPosition, 
			  const int& oldyPosition) = 0;
    
    // modify the matrix such that x mouse offset will
    // rotate trackball about the Z axis
    // give - the x mouse current and previous positions
    // result - matrix is offset by the delta
    virtual void applyZrot(const int& xPosition, 
			  const int& oldxPosition) = 0;
    
    // modify the matrix such that x,y mouse offsets will
    // rotate by the derived-class-defined method.
    // give - the x,y mouse current and previous positions
    // result - matrix is offset by the deltas
    virtual void applyGeneralRot(const int& xPosition, 
			    const int& yPosition, 
			    const int& oldXposition, 
			    const int& oldYposition) = 0;
        
    // set the window dimensions, 
    // (some trackballs map geometrically to the viewport)
    enum ViewportSide
    {
       top, bottom, left, right
    };
    void		setViewport( const ViewportSide& side, const int& dimension );
    void		getViewport( const ViewportSide& side, int& dimension ) const;
    
    // set translation
    virtual void setTranslation(const Vec3<float>& trans) = 0;    
    
    // set translation
    virtual void setTranslation(const float& x, const float& y, const float& z) = 0;    
    
    //set rotation
    virtual void setRotation(const Matrix4f& trans) = 0;
    
    //set rotation
    virtual void setRotation(const float& rad, const float& x, const float& y, const float& z) = 0;
    
    // get the lookat matrix
    const Matrix4f& getMatrix() const;
    Matrix4f& getMatrix();
    const Matrix4f& getMatrix(int a, int i, int r, int l,int f, int b) const;
    Matrix4f& getMatrix(int a, int i, int r, int l,int f, int b);
    
    // places the center of rotation at the center of the screen
    // -- clears the x and y translational component
    virtual void	center() { _lookAtMatrix.makeIdentity(); }

    void setSens( float sens ) { _sensativity = sens; }

    Matrix4f& setADirection();
    Matrix4f& setIDirection();
    Matrix4f& setRDirection();
    Matrix4f& setLDirection();
    Matrix4f& setFDirection();
    Matrix4f& setBDirection();

    const Matrix4f& setADirection() const;
    const Matrix4f& setIDirection() const;
    const Matrix4f& setRDirection() const;
    const Matrix4f& setLDirection() const;
    const Matrix4f& setFDirection() const;
    const Matrix4f& setBDirection() const;

protected:
    Matrix4f _lookAtMatrix;
   
    Vec3<float>  _lOldPosition; //2D
    Vec3<float>  _lPosition;	//2D
    Vec3<float>  _rOldPosition;	//2D
    Vec3<float>  _rPosition;	//2D
    Vec3<float>  _mOldPosition;	//2D
    Vec3<float>  _mPosition;	//2D
    
    int _windowLeft;
    int _windowRight;
    int _windowTop;
    int _windowBottom;
    float _sensativity;

};


inline TrackBall::TrackBall() : _sensativity(0.1f), _lookAtMatrix( Matrix4f::identity() )
{
}

// set the window dimensions, 
// (some trackballs map geometrically to the viewport)
inline void TrackBall::setViewport( const TrackBall::ViewportSide& side, const int& dimension )
{
   switch (side)
   {
      case left: _windowLeft = dimension; break;
      case top: _windowTop = dimension; break;
      case right: _windowRight = dimension; break;
      case bottom: _windowBottom = dimension; break;
      default:
            assert( false && "invalid param to TrackBall::setViewport" );
   }
}
inline void TrackBall::getViewport( const ViewportSide& side, int& dimension ) const
{
   switch (side)
   {
      case left: dimension = _windowLeft; break;
      case top: dimension = _windowTop; break;
      case right: dimension = _windowRight; break;
      case bottom: dimension = _windowBottom; break;
      default:
            assert( false && "invalid param to TrackBall::getViewport" );
   }
}

// get the lookat matrix
//inline const Matrix4f& TrackBall::getMatrix() const
inline const Matrix4f& TrackBall::getMatrix( ) const
{
    return _lookAtMatrix;
}
inline Matrix4f& TrackBall::getMatrix( ) 
{
    return _lookAtMatrix;
}
inline const Matrix4f& TrackBall::getMatrix(int a, int i, int r, int l, int f, int b) const
{
	if(a) 
		this ->setADirection();
	else if(i) 
		setIDirection();
	else if(r) 
		setRDirection();
	else if(l) 
		setLDirection();
	else if(f) 
		setFDirection();
	else if(b) 
		setBDirection();

    return _lookAtMatrix;
}

// get the lookat matrix
//inline Matrix4f& TrackBall::getMatrix()
inline Matrix4f& TrackBall::getMatrix(int a, int i, int r, int l, int f, int b) 
{
	if(a) 
		setADirection(); // superior
	else if(i) 
		setIDirection(); // inferior
	else if(r) 
		setRDirection(); // right
	else if(l) 
		setLDirection(); // left
	else if(f) 
		setFDirection(); // anterior
	else if(b)
		setBDirection(); // posterior
    return _lookAtMatrix;
}
inline Matrix4f& TrackBall::setADirection()
{
	_lookAtMatrix.set(-1.0, 0.0, 0.0, 0.0, 
			  0.0, 1.0, 0.0, 0.0,
			  0.0, 0.0, -1.0, 0.0,
			  0.0, 0.0, 0.0, 1.0);
    	return _lookAtMatrix;
}
inline Matrix4f&  TrackBall::setIDirection()
{
	_lookAtMatrix.set(1.0, 0.0, 0.0, 0.0, 
			  0.0, 1.0, 0.0, 0.0,
			  0.0, 0.0, 1.0, 0.0,
			  0.0, 0.0, 0.0, 1.0);
    	return _lookAtMatrix;
}
inline Matrix4f&  TrackBall::setRDirection()
{
	_lookAtMatrix.set(0.0, 1.0, 0.0, 0.0, 
			  0.0, 0.0, -1.0, 0.0,
			  -1.0,0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0, 1.0);
    	return _lookAtMatrix;
}
inline Matrix4f&  TrackBall::setLDirection()
{
	_lookAtMatrix.set(0.0, -1.0, 0.0, 0.0, 
			  0.0, 0.0, -1.0, 0.0,
			  1.0, 0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0, 1.0);
    	return _lookAtMatrix;
}
inline Matrix4f&  TrackBall::setFDirection()
{
	_lookAtMatrix.set(1.0, 0.0, 0.0, 0.0, 
			  0.0, 0.0, -1.0, 0.0,
			  0.0, 1.0, 0.0, 0.0,
			  0.0, 0.0, 0.0, 1.0);
    	return _lookAtMatrix;
}
inline Matrix4f&  TrackBall::setBDirection()
{
	_lookAtMatrix.set(-1.0, 0.0, 0.0, 0.0, 
			  0.0, 0.0, -1.0, 0.0,
			  0.0, -1.0, 0.0, 0.0,
			  0.0, 0.0, 0.0, 1.0);
    	return _lookAtMatrix;
}

#endif
