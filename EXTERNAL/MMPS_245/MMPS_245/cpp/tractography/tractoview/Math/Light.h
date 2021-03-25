
#ifndef LIGHT_OBJECT
#define LIGHT_OBJECT

#include "Vec4.h"
#include "ColorRGBA.h"

class Light
{
public:
     enum LightType
     {
       ambient, diffuse, specular
     };
     Light();

     // Light::setColor accepts a token for the color attribute to set (-
     // Light::ambient, Light::diffuse, or Light::specular) and three floating point
     // values (r, g, and b) in the range [0.0 .. 1.0] defining values for the
     // red, green, and blue components of the indicated attribute of the light
     // source.  By default, the r, g, and b values are all 1.0.
     void              setColor( const LightType& which, const float& r = 1.0f, const float& g = 1.0f, const float& b = 1.0f );
     void              getColor( const LightType& which, float& r, float& g, float& b ) const;

     // Light::setAtten sets the attenuation parameters of the pfLight. The
     // light intensity is scaled at each vertex by:
     //
     //      1.0 / (constant + linear * dist + quadratic * dist^2)
     //
     // where 'dist' is the distance from the light position to the lit vertex.
     // Note that 'dist' is 1.0 for infinite light sources. The default
     // attenuation values are constant = 1.0, linear = 0.0, quadratic = 0.0,
     // i.e., light attenuation is disabled.
     void              setAtten( const float& constant = 1.0f, const float& linear = 0.0f, const float& quadratic = 0.0f );
     
     // Light::getAtten returns the
     //     attenuation parameters of the pfLight in constant, linear, and quadratic.
     void              getAtten( float& constant, float& linear, float& quadratic ) const;

     // Light::setPos receives four floating point values to set the x, y, z,
     // and w, coordinates for the position of the light source.  Typically, the
     // homogeneous coordinate w is 0.0 to indicate that the light position is
     // infinitely far from the origin in the direction (x, y, z).  Local light
     // sources are specified by a non-zero value for w and usually incur a
     // performance penalty.
     void              setPos( const float& x, const float& y, const float& z, const float& w );

     // Light::getPos copies the x, y, z and w
     // coordinates of the light source into the parameters x, y, z and w,
     // respectively.
     void              getPos( float& x, float& y, float& z, float& w ) const;


     // Light::setSpotDir specifies the direction in which a spot light source
     // emits its light.  It receives three floating point values, x, y, and z,
     // specifying the x, y, and z direction vectors.  
     // It is significant only when is spread is not 180, which it is by default. 
     // The default direction is (0,0,-1).
     void              setSpotDir( const float& x = 0.0f, const float& y = 0.0f, const float& z = -1.0f );

     // Light::getSpotDir copies the x, y, and z direction vectors into 
     // the parameters x, y, and z.
     void              getSpotDir( float& x, float& y, float& z ) const;

     // Light::setSpotCone specifies the exponent and spread of the spot light
     // cone, and receives two floating point values, f1 and f2, to set the
     // exponent for the intensity, and the spread of the cone, respectively.
//
// Exponent is a single integer or floating-point value that specifies the intensity distribution of the
// light.  Integer and floating-point values are mapped directly.  Only values in the range [0,128] are
// accepted.
//
// Effective light intensity is attenuated by the cosine of the angle between the direction of the light and
// the direction from the light to the vertex being lighted, raised to the power of the spot exponent.
// Thus, higher spot exponents result in a more focused light source, regardless of the spot cutoff angle
// (see next paragraph).  The default spot exponent is 0, resulting in uniform light distribution.
//
// Spread (cutoff) is a single integer or floating-point value that specifies the maximum spread angle of a light
// source.  Integer and floating-point values are mapped directly.  Only values in the range [0,90], and the
// special value 180, are accepted.  If the angle between the direction of the light and the direction
// from the light to the vertex being lighted is greater than the spot cutoff angle, the light is completely
// masked.  Otherwise, its intensity is controlled by the spot exponent and the attenuation factors.  The
// default spot cutoff is 180, resulting in uniform light distribution.
     void              setSpotCone( const float& exponent = 0.0f, const float& spread = 180.0f );

     // Light::getSpotCone copies the current exponent and spread of the cone
     // into the parameters f1 and f2.
     void              getSpotCone( float& exponent, float& spread ) const;

     void              on();
     void              off();
     bool              isOn() const;

     //int               pfGetCurLights( pfLight *lights[PF_MAX_LIGHTS] );

     void setNumber( int num ) { mLightnumber = num; }
     int number() const { return mLightnumber; }
private:
    int                 mLightnumber;
    bool                mIsOn;

    // color properties
    ColorRGBA           mAmbient;
    ColorRGBA           mDiffuse;
    ColorRGBA           mSpecular;

    // set position[3] == 1 for positional light
    // set position[3] == 0 for directional light
    Vec4<float>         mPosition;

    // attenuation
    float               mConstantAttenuation;
    float               mLinearAttenuation;
    float               mQuadraticAttenuation;

//spotlight
private:
    //Angle of the light cone
    float               mCutoff;

    //direction, default is 0,0,-1
    Vec3<float>         mDirection;

    float               mExponent;
};


#endif
