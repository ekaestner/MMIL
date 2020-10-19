#ifndef LINE_
#define LINE_

#include "Vec3.h"
#include "Matrix4f.h"

class Line
{
public:
	// a line segment is defined by two points
	Vec3<float> p1, p2;

	// multiply matrix * self, returns the resulting line
	void multiply( const Matrix4f& mat, Line& result ) const;
};

// multiply matrix * self, returns the resulting line
inline void Line::multiply( const Matrix4f& mat, Line& result ) 
{
	result.origin = mat * this->origin;
	result.direction = mat * this->direction;
}
