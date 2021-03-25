
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#ifndef MOUSE_INCLUDED
#define MOUSE_INCLUDED

#include <assert.h>

class Mouse
{
// Types
public:
    enum Button
    {
	LEFT, MIDDLE, RIGHT
    };
    
    enum EdgeTriggerState
    {
	ONDOWN, DOWN, ONUP, UP
    };
    
    enum BinaryState
    {
	ON, OFF
    };
    
// Required methods for promised functionality:
public:
    //Constructor.
    Mouse();
    
    //: Call this on every Mouse event
    // Mouse events include position change, and button change
    void		    updateEdgeStates();
    
    // Set the binary state of the button, ON or OFF
    // call this on every Mouse event
    // (a button change or position change)
    // give - b:     Mouse::LEFT or Mouse::MIDDLE or Mouse::RIGHT
    // give - state: Mouse::ON or Mouse::OFF
    // result - updates internal data, 
    //          after this call the edge trigger states are valid.
    void		    setState(const Button& b, 
				const Mouse::BinaryState& state);

    // Set the position of the Mouse object
    // Call this on a Mouse-Move event (when the mouse moves)
    // give - x and y mouse coordinates
    void		    setPosition(const int& x, const int& y);
    
// Preferred methods - total flexibility
public:
    // Get the edge triggered state of the button
    // returns - Mouse::ONDOWN or Mouse::ONUP, etc...
    const EdgeTriggerState& leftEdgeState() const;
    const EdgeTriggerState& middleEdgeState() const;
    const EdgeTriggerState& rightEdgeState() const;
    
    // Get the binary state of the button
    // returns - Mouse::ON or Mouse::OFF
    const BinaryState&	    leftBinaryState() const;
    const BinaryState&	    middleBinaryState() const;
    const BinaryState&	    rightBinaryState() const;
    
    // return the x position of the mouse.
    const int&		    x() const;
    
    // return the y position of the mouse.
    const int&		    y() const;
    
    // return the previous x position of the mouse.
    const int&		    xOld() const;
    
    // return the previous y position of the mouse.
    const int&		    yOld() const;
    
    // get the change between mouse's previous and current x positions
    const int&		    dx() const;
    
    // get the change between mouse's previous and current y positions
    const int&		    dy() const;

// Alternate Methods
public:
    // Get the mouse position
    // result - x and y are in viewport coordinates
    void		    getPosition(int& x, int& y) const;
    
    // Get the previous position of mouse
    void		    getPreviousPosition(int& x, int& y) const;

    // Get the change between mouse's previous and current positions
    void getOffset(int& xOffset, int& yOffset) const;

// Private member data
private:
    int			    _x, _xOld, _xOffset;
    int			    _y, _yOld, _yOffset;
    BinaryState		    _leftBinary, _middleBinary, _rightBinary;
    EdgeTriggerState	    _leftEdge, _middleEdge, _rightEdge;
};



//Constructor.
inline Mouse::Mouse() : _leftBinary(Mouse::OFF), 
		_middleBinary(Mouse::OFF), 
		_rightBinary(Mouse::OFF), 
		 _leftEdge(Mouse::UP), 
		 _middleEdge(Mouse::UP), 
		 _rightEdge(Mouse::UP), 
		 _x(0), _y(0)
{
}

// Get the state of the button
inline const Mouse::EdgeTriggerState& Mouse::leftEdgeState() const
{
    return _leftEdge;
}
inline const Mouse::EdgeTriggerState& Mouse::middleEdgeState() const
{
    return _middleEdge;
}
inline const Mouse::EdgeTriggerState& Mouse::rightEdgeState() const
{
    return _rightEdge;
}
inline const Mouse::BinaryState& Mouse::leftBinaryState() const
{
    return _leftBinary;
}
inline const Mouse::BinaryState& Mouse::middleBinaryState() const
{
    return _middleBinary;
}
inline const Mouse::BinaryState& Mouse::rightBinaryState() const
{
    return _rightBinary;
}

// get the mouse position
inline void Mouse::getPosition(int& x, int& y) const
{
    x = _x;
    y = _y;
}

// return the x position of the mouse.
inline const int& Mouse::x() const
{
    return _x;
}

// return the y position of the mouse.
inline const int& Mouse::y() const
{
    return _y;
}


// return the previous x position of the mouse.
inline const int& Mouse::xOld() const
{
    return _xOld;
}

// return the previous y position of the mouse.
inline const int& Mouse::yOld() const
{
    return _yOld;
}

// get the change between mouse's previous and current x positions
inline const int& Mouse::dx() const
{
    return _xOffset;
}

// get the change between mouse's previous and current y positions
inline const int& Mouse::dy() const
{
    return _yOffset;
}
  
inline void Mouse::setPosition(const int& x, const int& y)
{
    //save the old position
    _xOld = _x;
    _yOld = _y;
    
    //get the current position
    _x = x;
    _y = y;
    
    //get the offset the mouse has moved
    _xOffset = _x - _xOld;
    _yOffset = _y - _yOld;
}

// get the previous position of mouse
inline void Mouse::getPreviousPosition(int& x, int& y) const
{
    x = _xOld;
    y = _yOld;
}

// get the offset between mouse's previous and current positions
inline void Mouse::getOffset(int& xOffset, int& yOffset) const
{
    xOffset = _xOffset;
    yOffset = _yOffset;
}

inline void Mouse::setState(const Button& b, 
			const Mouse::BinaryState& state)
{
    //set the binary state for button 'b'
    switch( b )
    {
    case LEFT:   _leftBinary = state;  break;
    case MIDDLE: _middleBinary = state;break;
    case RIGHT:  _rightBinary = state; break;
    default: assert( false );
    }
}

//: Call this on every Mouse event
// Mouse events include position change, and button change
inline void Mouse::updateEdgeStates()
{
    //find out if...
    switch( _leftBinary )
    {
    case ON:
    //button has been held down or it is just down.
	switch( _leftEdge )
	{	
	case ONUP:
	case UP:    _leftEdge = ONDOWN;  break;
	case ONDOWN:_leftEdge = DOWN; break;
	}
	break;
	
    case OFF:
    //button has been up or it is just up.
	switch( _leftEdge )
	{
	case ONDOWN:
	case DOWN:  _leftEdge = ONUP; break;
	case ONUP:  _leftEdge = UP; break;
	}
	break;
	
    default: assert(false);
    }
    
    //find out if...
    switch( _middleBinary )
    {
    case ON:
    //button has been held down or it is just down.
	switch( _middleEdge )
	{	
	case ONUP:
	case UP:    _middleEdge = ONDOWN;  break;
	case ONDOWN:_middleEdge = DOWN; break;
	}
	break;
	
    case OFF:
    //button has been up or it is just up.
	switch( _middleEdge )
	{
	case ONDOWN:
	case DOWN:  _middleEdge = ONUP; break;
	case ONUP:  _middleEdge = UP; break;
	}
	break;
	
    default: assert(false);
    }
    
    //find out if...
    switch( _rightBinary )
    {
    case ON:
    //button has been held down or it is just down.
	switch( _rightEdge )
	{	
	case ONUP:
	case UP:    _rightEdge = ONDOWN;  break;
	case ONDOWN:_rightEdge = DOWN; break;
	}
	break;
	
    case OFF:
    //button has been up or it is just up.
	switch( _rightEdge )
	{
	case ONDOWN:
	case DOWN:  _rightEdge = ONUP; break;
	case ONUP:  _rightEdge = UP; break;
	}
	break;
	
    default: assert(false);
    }
}

#endif
