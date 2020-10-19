
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#ifndef _COLORf_H_
#define _COLORf_H_

#include <iostream.h>
#include "Vec4.h"
#include "Defines.h" // kev::max kev::min

   class ColorRGBA
   {
   protected:
	   // ColorRGBA "has" 4 channels.
	   Vec4<float> _channel;
	   void validate();

   public:
	   //: Channel names.  
	   // ColorRGBA is in RGBA format always - use getHSV or etc...  to convert formats.
	   enum Channel
	   {
		   HUE,SAT,VAL,
		   RED,GREEN,BLUE
	   };

	   //: Default Constructor
	   ColorRGBA();

	   //: Initialize the color with floating point numbers between [0.0f, 1.0f]
	   ColorRGBA( const float& red, const float& green, const float& blue, const float& alpha = 1.0f );
	   
	   //: Initialize the color with integer numbers between [0,255]
	   ColorRGBA( const unsigned int& red, const unsigned int& green, const unsigned int& blue, const unsigned int& alpha = 255);
	   
	   //: Initialize the color with integer numbers between [0,255]
	   ColorRGBA( const int& red, const int& green, const int& blue, const int& alpha = 255);
	   
	   //: Initialize the color with one integer number 0x00rrggbb where:
	   // r - unsigned character representing the red channel
	   // g - unsigned character representing the green channel
	   // b - unsigned character representing the blue channel
	   // alpha - floating point alpha value between [0,1]
	   ColorRGBA( const unsigned int& rgb, const float& alpha = 1.0f );
	   
	   //: Copy constructor
	   ColorRGBA( const ColorRGBA& color );

	   //: Copy constructor for vector
	   ColorRGBA( const Vec4<float>& color );
	   
	   
	   

	   //: Get the color values between 0-1.
	   void					get( float &red, float &green, float &blue, float &alpha ) const;
	   void					get( float &red, float &green, float &blue ) const;

           // get a copy of the ColorRGBA
	   void					get( ColorRGBA &color ) const;

	   //: Get the color values between 0-255.
	   void					getInteger ( unsigned int &red, unsigned int &green, unsigned int &blue, unsigned int &alpha ) const;
	   void					getInteger ( unsigned int &red, unsigned int &green, unsigned int &blue ) const;
	   void					getInteger ( int &red, int &green, int &blue, int &alpha ) const;
	   void					getInteger ( int &red, int &green, int &blue ) const;
	   void					getInteger ( unsigned int &rgb ) const;
	   

	   //: Access to the "data" in this class.  (read only)
	   //  this function can be used as if it "is" color channel data. (that's why there is no "get")
	   inline const float *	data() const { return _channel.data(); }

	   //: Access to the "Vec4<float>" vector in this class, 
	   //  useful for arithmetic operations, this vector holds the 4 channel color data.
	   inline float *			data() { return _channel.data(); }

	   //: Access to the "Vec4<float>" vector in this class, (read only)
	   //  useful for arithmetic operations, this vector holds the 4 channel color data.
	   inline const Vec4<float>&	vector() const { return _channel; }

	   //: Access to the "Vec4<float>" vector in this class, 
	   //  useful for arithmetic operations, this vector holds the 4 channel color data.
	   inline Vec4<float>&			vector() { return _channel; }

	   //: Get the hsv. Hue is defined as [0,360) or SL_UNDEFINED_HUE_COLOR.
	   //  All others are defined as [0,1];
	   void					getHSV( float &hue, float &saturation, float &value ) const;

	   //: Assign the color.
	   const ColorRGBA &		operator=( const ColorRGBA& color );
	   const ColorRGBA &		operator=( const unsigned int& color );
	   const ColorRGBA &		operator=( const int& color );

	   //: Equality operator.
	   bool					operator==( const ColorRGBA &color ) { return this->_channel == color._channel; }

      // Linear Interpolation between two colors, includes alpha.
	   void			      Lerp(const ColorRGBA& from, const ColorRGBA& to, 
					               const float& lerp )
      {
         Vec4<float> offset = to._channel - from._channel;
         this->_channel = from._channel + offset * lerp;
      }


   // Methods to access channel data: in floating point between [0,1]
   public:
	   // Conversion operator (read only) - from SlColor to float*
	   inline					operator const float* () const { return _channel.data(); }

	   // Conversion operator - from SlColor to float*
	   inline					operator float* () { return _channel.data(); }
	   
	   // Index an individual color channel [0,3] --> red/green/blue/alpha
	   inline float &			operator[] ( const int& i ) { return _channel[i]; }
	   // Index an individual color channel [0,3] (read only) --> red/green/blue/alpha
	   inline const float &	operator[] ( const int& i ) const { return _channel[i]; }
	   
	   //: Access the red channel
	   inline float&			red() { return _channel[0]; }
	   //: Access the red channel (read only)
	   inline const float&		red() const { return _channel[0]; }
	   
	   //: Access the green channel
	   inline float&			green() { return _channel[1]; }
	   //: Access the green channel (read only)
	   inline const float&		green() const { return _channel[1]; }
	   
	   //: Access the blue channel
	   inline float &			blue() { return _channel[2]; }
	   //: Access the blue channel (read only)
	   inline const float&		blue() const { return _channel[2]; }
	   
	   //: Access the alpha channel
	   inline float &			alpha() { return _channel[3]; }
	   //: Access the alpha channel (read only)
	   inline const float&		alpha() const { return _channel[3]; }

   public:
	   //: Input/output.
	   friend  ostream &		operator<< ( ostream &out, const ColorRGBA &color );
	   friend  istream &		operator>> ( istream &in, ColorRGBA &color );

	   //: Set the color values between 0.0f - 1.0f.
	   void					set( const float& red, const float&  green, const float&  blue, const float& alpha = 1.0f);
	   void					set( const ColorRGBA& color );

	   //: Set the color values between 0 - 255.
	   void					setInteger ( const unsigned int& red, const unsigned int& green, const unsigned int& blue, const unsigned int& alpha = 255 );
	   void					setInteger ( const int& red, const int& green, const int& blue, const int& alpha = 255 );
	   void					setInteger ( const unsigned int& rgb, const float& alpha = 1.0f );
	   
	   //: Set the color in Hue Saturation Brightness format.
	   //  All values are defined in floating point between [0,1];
	   void 					setHSV ( const float &hue, const float &saturation, const float &value, const float& alpha = 1.0f );

	   ///////////////////////////////////////////////////////////////////////////////////
	   //Convert Hue-Saturation-Value format to Red-Green-Blue format
	   //
	   //r, g, b, s, v are from 0 to 1
	   //h is from 0 to 360
	   static void				HSV2RGB(const float &h, const float &s, const float &v, float &r, float &g, float &b);
	   static void				HSV2RGB(const ColorRGBA &hsv, ColorRGBA &rgb);

	   ///////////////////////////////////////////////////////////////////////////////////
	   //Convert Red-Green-Blue format to Hue-Saturation-Value format
	   //
	   //r, g, b, s, v are from 0 to 1
	   //h is from 0 to 360
	   static void				RGB2HSV(const float &r, const float &g, const float &b, float &h, float &s, float &v);
	   static void				RGB2HSV(const ColorRGBA &rgb, ColorRGBA &hsv);

	   //function returns "color" in the hex format specified with bpp (see below COLORXBIT functions)
	   // bpp can be 15, 16, 24, or 32.
	   // red, green, blue, can be from [0, 1]
	   static void				RGB2HEX(const float &red, const float &green, const float &blue, const int &bpp, int &color);
	   static void				RGB2HEX(const ColorRGBA &rgb, const int &bpp, int &color);

	   //function returns "color" in the hex format specified with bpp (see below COLORXBIT functions)
	   // bpp can be 15, 16, 24, or 32.
	   // red, green, blue, can be from [0, 1]
	   static void				HSV2HEX(const float &hue, const float &saturation, const float &value, const int &bpp, int &color);
	   static void				HSV2HEX(const ColorRGBA &hsv, const int &bpp, int &color);

	   ///////////////////////////////////////////////////////////////////////////////////
	   //Conversion of RGB to 15bpp, (32768 colors possible: "High Color")
	   //
	   //Info:
	   // - 5 bits stored for each color
	   // - 32 reds, 32 greens, 32 blues.
	   // - The blue mask is 0x001F, the green mask is 0x03E0, and the red mask is 0x7C00. 
	   //
	   //Usage:
	   // - Pass in 0x1f, 0x1f, 0x1f as the highest values.
	   //
	   inline static short int COLOR15BIT(const int &r, const int &g, const int &b)
	   {
		   return (r<<10 & 0x7C00) | (g<<5 & 0x03E0) | (b & 0x001F);
	   }

	   ///////////////////////////////////////////////////////////////////////////////////
	   //Conversion of RGB to 16bpp, (65536 colors possible: "High Color")
	   //
	   //Info:
	   // - 6 bits for the green, 5 red, 5 blue
	   // - 32 reds, 64 greens, 32 blues.
	   // - The blue mask is 0x001F, the green mask is 0x07E0, and the red mask is 0xF800. 
	   //
	   //Usage:
	   // - Pass in 0x1f, 0x3f, 0x1f as the highest values.
	   //
	   inline static int		COLOR16BIT(const int &r, const int &g, const int &b)
	   {
		   return (r<<11 & 0xF800) | (g<<5 & 0x07E0) | (b & 0x001F);
	   }

	   ///////////////////////////////////////////////////////////////////////////////////
	   //Conversion of RGB to 24bpp, (16 million colors possible: "True Color")
	   //
	   //Info:
	   // - 8 bits for each color RGB
	   // - 256 reds, 256 greens, 256 blues.
	   // - The blue mask is 0x0000FF, the green mask is 0x00FF00, the red mask is 0xFF0000. 
	   //
	   //Usage:
	   // - Pass in 0xFF, 0xFF, 0xFF as the highest values.
	   //
	   inline static int		COLOR24BIT(const int &r, const int &g, const int &b)
	   {
		   return (r<<16 & 0xFF0000) | (g<<8 & 0x00FF00) | (b & 0x0000FF);
	   }

	   ///////////////////////////////////////////////////////////////////////////////////
	   //Conversion of RGB to 32bpp, (16 million colors possible: "True Color")
	   //
	   //Info:
	   // - 8 bits for each color RGB
	   // - 256 reds, 256 greens, 256 blues.
	   // - The blue mask is 0x000000FF, the green mask is 0x0000FF00, the red mask is 0x00FF0000. 
	   //
	   //Usage:
	   // - Pass in 0xFF, 0xFF, 0xFF as the highest values.
	   //
	   inline static int		COLOR32BIT(const int &r, const int &g, const int &b)
	   {
		   return (r<<16 & 0x00FF0000) | (g<<8 & 0x0000FF00) | (b & 0x000000FF);
	   }
   };


   //////////////////////////////////////////////////////////////////////////
   //
   //	Algorithm from: "Computer Graphics, Principles and Practice", 
   //	2nd edition in "C", Foley, van Dam, page 593.
   //
   //  Convert Hue/Saturation/Value to Red/Green/Blue
   //  r, g, b, h, s, v are from 0 to 1
   //
   //////////////////////////////////////////////////////////////////////////
   inline void ColorRGBA::HSV2RGB(const float &h, const float &s, const float &v, float &r, float &g, float &b)
   {
	   float _h = (h * 360.0f);
	   float _s = s, _v = v, f, p, q, t;
	   int i;

	   //the color is on the black and white center line
	   if(_s == 0.0f)
	   {
		   //achromatic color: there is no hue
		   if(_h == 0.0f) // == undefined) 
		   {
			   r = g = b = _v;
		   } //Error!!, this is a good approx though..
			   else
			   r = g = b = _v;
	   } else {
		   //chromatic color: s != 0, so there is a hue.

		   if(_h == 360.0f)
			   _h = 0.0f;
		   _h /= 60.0f;  //h is now in [0,6)  ??
		   i = (int) floor(_h);
		   f = _h - i;	//fractional part of h
		   p = _v * (1.0f - _s);
		   q = _v * (1.0f - (_s * f));
		   t = _v * (1.0f - (_s * (1.0f - f)));
		   switch(i)
		   {
		   case 0: r = _v; g =  t; b =  p; break;
		   case 1: r =  q; g = _v; b =  p; break;
		   case 2: r =  p; g = _v; b =  t; break;
		   case 3: r =  p; g =  q; b = _v; break;
		   case 4: r =  t; g =  p; b = _v; break;
		   case 5: r = _v; g =  p; b =  q; break;
		   }
	   }
   }

   /////////////////////////////////////////////////////////////////////////////////
   //
   //	Algorithm from: "Computer Graphics, Principles and Practice", 
   //	2nd edition in "C", Foley, van Dam, page 593.
   //
   //  Convert Red/Green/Blue to Hue/Saturation/Value
   //  r, g, b, h, s, v are from 0 to 1
   //
   /////////////////////////////////////////////////////////////////////////////////
   inline void ColorRGBA::RGB2HSV(const float &r, const float &g, const float &b, float &h, float &s, float &v)
   {
	   float delta = 0.0f;
	   float max = kev::max(r, g, b);
	   float min = kev::min(r, g, b);

	   //get the value...
	   v = max;

	   // get the saturation...
	   s = (max != 0.0f) ? ((max - min) / max) : 0.0f;
	   //if( max != 0.0f )
	   //	s = (max - min) / max;


	   //chromatic case: saturation is not 0, so determine hue.
	   if(s != 0.0f)
	   {
		   delta = max - min;

		   if(r == max)		
			   h = (g - b)/delta;
		   else if(g == max)	
			   h = 2.0f + (b - r)/delta;
		   else if(b == max)	
			   h = 4.0f + (r - g)/delta;

		   h *= 60.0f;
		   if(h < 0.0f)
			   h += 360.0f;
	   } else 
		   //== undefined
		   h = 0.0f; 

	   //scale h to be between 0 and 1
	   h /= 360.0f;
   }

   //function returns "color" in the hex format specified with bpp
   // bpp can be 15, 16, 24, or 32.
   // red, green, blue, can be from [0, 1]

   inline void	ColorRGBA::RGB2HEX(const float &red, const float &green, const float &blue, const int &bpp, int &color)
   {
	   int r, g, b;

	   color = r = g = b = 0;

	   switch( bpp )
	   {
	   case 15:
		   r = static_cast<int>(red * 31.0f);
		   g = static_cast<int>(green * 31.0f);
		   b = static_cast<int>(blue * 31.0f);
		   color = (short int)ColorRGBA::COLOR15BIT(r, g, b); 
		   break;
	   case 16: 
		   r = static_cast<int>(red * 31.0f);
		   g = static_cast<int>(green * 63.0f);
		   b = static_cast<int>(blue * 31.0f);
		   color = (short int)ColorRGBA::COLOR16BIT(r, g, b); 
		   break;
	   case 24:
	   case 32: 
		   r = static_cast<int>(red * 255);
		   g = static_cast<int>(green * 255);
		   b = static_cast<int>(blue * 255);
		   color = ColorRGBA::COLOR32BIT(r, g, b); 
		   break;
	   }
   }

   //function returns "color" in the hex format specified with bpp
   // bpp can be 15, 16, 24, or 32.
   // red, green, blue, can be from [0, 1]

   inline void	ColorRGBA::HSV2HEX(const float &hue, const float &saturation, const float &value, const int &bpp, int &color)
   {
	   float red, green, blue;

	   //convert hsv value to rgb.
	   ColorRGBA::HSV2RGB( hue, saturation, value, red, green, blue );

	   //convert rgb value to hex
	   ColorRGBA::RGB2HEX( red, green, blue, bpp, color );
   }

   //////////////////////////////////////////////////////////////////
   // Cloned functions:
   inline void ColorRGBA::HSV2RGB(const ColorRGBA &hsv, ColorRGBA &rgb)
   {
	   ColorRGBA::HSV2RGB( hsv[0], hsv[1], hsv[2], rgb[0], rgb[1], rgb[2] );
   }

   inline void ColorRGBA::RGB2HSV(const ColorRGBA &rgb, ColorRGBA &hsv)
   {
	   ColorRGBA::RGB2HSV( rgb[0], rgb[1], rgb[2], hsv[0], hsv[1], hsv[2] );
   }

   inline void ColorRGBA::RGB2HEX(const ColorRGBA &rgb, const int &bpp, int &color)
   {
	   ColorRGBA::RGB2HEX( rgb[0], rgb[1], rgb[2], bpp, color );
   }

   inline void ColorRGBA::HSV2HEX(const ColorRGBA &hsv, const int &bpp, int &color)
   {
	   ColorRGBA::HSV2HEX( hsv[0], hsv[1], hsv[2], bpp, color );
   }
   //////////////////////////////////////////////////////////////////


#endif
