
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "ColorRGBA.h"

//	Default Constructor
ColorRGBA::ColorRGBA()
{
	_channel[0] = 0.0f;
	_channel[1] = 0.0f;
	_channel[2] = 0.0f;
	_channel[3] = 1.0f;
}


//	INTEGER Constructor.
ColorRGBA::ColorRGBA( const unsigned int& red, const unsigned int& green, const unsigned int& blue, const unsigned int& alpha )
{
	float factor = 1.0f / 255.0f;

	_channel[0] = ((float) red) * factor;
	_channel[1] = ((float) green) * factor;
	_channel[2] = ((float) blue) * factor;
	_channel[3] = ((float) alpha) * factor;

	this->validate();
}

//	INTEGER Constructor.
ColorRGBA::ColorRGBA( const int& red, const int& green, const int& blue, const int& alpha )
{
	float factor = 1.0f / 255.0f;

	_channel[0] = ((float) red) * factor;
	_channel[1] = ((float) green) * factor;
	_channel[2] = ((float) blue) * factor;
	_channel[3] = ((float) alpha) * factor;

	this->validate();
}


//	Constructor.
ColorRGBA::ColorRGBA( const float& red, const float& green, const float& blue, const float& alpha )
{
	_channel[0] = red;
	_channel[1] = green;
	_channel[2] = blue;
	_channel[3] = alpha;

	this->validate();
}

//	Copy Constructor.
ColorRGBA::ColorRGBA( const ColorRGBA& color ) 
{
	_channel[0] = color.red();
	_channel[1] = color.green();
	_channel[2] = color.blue();
	_channel[3] = color.alpha();

	this->validate();
}

//: Copy constructor for vector
ColorRGBA::ColorRGBA( const Vec4<float>& color )
{
	_channel = color;

	this->validate();
}


//	Constructor.
ColorRGBA::ColorRGBA( const unsigned int& rgb, const float& alpha )
{
	unsigned int red =   ( ( rgb & 0x000000FF ) >> 0 );
	unsigned int green = ( ( rgb & 0x0000FF00 ) >> 8 );
	unsigned int blue =  ( ( rgb & 0x00FF0000 ) >> 16 );
	
	float factor = 1.0f / 255.0f;

	_channel[0] = ((float) red) * factor;
	_channel[1] = ((float) green) * factor;
	_channel[2] = ((float) blue) * factor;
	_channel[3] = alpha;
}

//  Get the color.
void ColorRGBA::getInteger ( unsigned int& rgb ) const 
{ 
	unsigned int red =   (unsigned int) (_channel[0] * 255.0f);
	unsigned int green = (unsigned int) (_channel[1] * 255.0f);
	unsigned int blue =  (unsigned int) (_channel[2] * 255.0f);

	rgb = ( red << 0 ) | ( green << 8 ) | ( blue << 16 );
}

//  Get the color.
void ColorRGBA::get( float &red, float &green, float &blue ) const 
{ 
	red =   _channel[0];
	green = _channel[1];
	blue =  _channel[2];
}

//  Get the color.
void ColorRGBA::get( float &red, float &green, float &blue, float &alpha ) const 
{ 
	red =   _channel[0];
	green = _channel[1];
	blue =  _channel[2];
	alpha = _channel[3];
}
	
//  Sets the ColorRGBA on the left (*this) of the "=" equal to the 
//  one on the right (color).
const ColorRGBA& ColorRGBA::operator=( const ColorRGBA& color )
{
	_channel[0] = color[0];
	_channel[1] = color[1];
	_channel[2] = color[2];
	_channel[3] = color[3];

	this->validate();

	return *this; 
}

//  Set the rgba value.
void ColorRGBA::set( const float& red, const float& green, const float& blue, const float& alpha ) 
{
	_channel[0] = red;
	_channel[1] = green;
	_channel[2] = blue;
	_channel[3] = alpha;

	this->validate();
}

//  Set the rgba value.
void ColorRGBA::setInteger( const int& red, const int& green, const int& blue, const int& alpha ) 
{
	float factor = 1.0f / 255.0f;

	_channel[0] = ((float) red) * factor;
	_channel[1] = ((float) green) * factor;
	_channel[2] = ((float) blue) * factor;
	_channel[3] = ((float) alpha) * factor;

	this->validate();
}

//  Set the rgba value.
void ColorRGBA::setInteger ( const unsigned int& red, const unsigned int& green, 
						   const unsigned int& blue, const unsigned int& alpha ) 
{
	float factor = 1.0f / 255.0f;

	_channel[0] = ((float) red) * factor;
	_channel[1] = ((float) green) * factor;
	_channel[2] = ((float) blue) * factor;
	_channel[3] = ((float) alpha) * factor;

	this->validate();
}

//  Set the rgba value.
void ColorRGBA::set( const ColorRGBA &color ) 
{
	_channel[0] = color[0];
	_channel[1] = color[1];
	_channel[2] = color[2];
	_channel[3] = color[3];

	this->validate();
}

//  Set the color.
void ColorRGBA::setInteger ( const unsigned int& rgb, const float& alpha ) 
{
	unsigned int red =   ( ( rgb & 0x000000FF ) >> 0 );
	unsigned int green = ( ( rgb & 0x0000FF00 ) >> 8 );
	unsigned int blue =  ( ( rgb & 0x00FF0000 ) >> 16 );
	
	float factor = 1.0f / 255.0f;

	_channel[0] = ((float) red) * factor;
	_channel[1] = ((float) green) * factor;
	_channel[2] = ((float) blue) * factor;
	_channel[3] = alpha;
}

	
//  Make sure the values are in range.
void ColorRGBA::validate() 
{
	if		( _channel[0] > 1.0f ) _channel[0] = 1.0f;
	else if ( _channel[0] < 0.0f ) _channel[0] = 0.0f;

	if		( _channel[1] > 1.0f ) _channel[1] = 1.0f;
	else if ( _channel[1] < 0.0f ) _channel[1] = 0.0f;

	if		( _channel[2] > 1.0f ) _channel[2] = 1.0f;
	else if ( _channel[2] < 0.0f ) _channel[2] = 0.0f;

	if		( _channel[3] > 1.0f ) _channel[3] = 1.0f;
	else if ( _channel[3] < 0.0f ) _channel[3] = 0.0f;
}

//  Output operator.
ostream& operator<<( ostream &out, const ColorRGBA &color )
{
	out << color[0] << " " 
		<< color[1] << " " 
		<< color[2] << " " 
		<< color[3];

	return out;
}

//  Input operator.
istream& operator>>( istream &in, ColorRGBA &color )
{
	in >> color[0] 
	   >> color[1] 
	   >> color[2] 
	   >> color[3];

	return in;
}

//  Get the hsv.
void ColorRGBA::getHSV ( float &hue, float &saturation, float &value ) const
{
	ColorRGBA::RGB2HSV ( _channel[0], _channel[1], _channel[2], hue, saturation, value );
}

//  Set the hsv.
void ColorRGBA::setHSV ( const float& hue, const float& saturation, const float& value, const float& alpha)
{
	ColorRGBA::HSV2RGB( hue, saturation, value, _channel[0], _channel[1], _channel[2] );
	_channel[3] = alpha;
}
