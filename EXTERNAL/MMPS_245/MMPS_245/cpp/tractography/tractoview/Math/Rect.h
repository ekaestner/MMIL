
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
#ifndef RECTANGLE_INCLUDED
#define RECTANGLE_INCLUDED

//specifies a rectangle in 2D space.
// top and bottom   are usually 0
// right and left   specify the width and height.

#include "Vec3.h"

class Rect
{
//construction
public:
    inline Rect() : top(0), bottom(0), right(0), left(0) {}
    Rect( const Rect &r);
    inline Rect(const float &l, const float &t, 
	 const float &r, const float &b ) 
	: top(t), bottom(b), right(r), left(l) 
	{}
    inline ~Rect(){};

//interface
public:
    float left;
    float top;
    float right;
    float bottom;
    
    
//operations
public:
	  Rect   operator+  ( const Vec3<float> &p );
 	  Rect   operator-  ( const Vec3<float> &p );
    const Rect & operator=  ( const Rect    &r );
    const Rect & operator+= ( const Vec3<float> &p );
    const Rect & operator-= ( const Vec3<float> &p );
};

inline Rect::Rect( const Rect &r)
{
    this->left	    = r.left;
    this->top	    = r.top;
    this->right	    = r.right;
    this->bottom    = r.bottom;
}

inline Rect Rect::operator+ ( const Vec3<float> &p )
{   
    return Rect( this->left	+ p[0], 
		 this->top	+ p[1], 
		 this->right	+ p[0], 
		 this->bottom	+ p[1] );
}

inline const Rect & Rect::operator= ( const Rect &r )
{
    this->left	    = r.left;
    this->top	    = r.top;
    this->right	    = r.right;
    this->bottom    = r.bottom;
    
    return *this;
}

inline const Rect & Rect::operator+= ( const Vec3<float> &p )
{
    *this = *this + p;
    
    return *this;
}

inline Rect Rect::operator- ( const Vec3<float> &p )
{
    return Rect( this->left	- p[0], 
		 this->top	- p[1], 
		 this->right	- p[0], 
		 this->bottom   - p[1] );
}

inline const Rect & Rect::operator-= ( const Vec3<float> &p )
{
    *this = *this - p;
    
    return *this;
}



#endif
