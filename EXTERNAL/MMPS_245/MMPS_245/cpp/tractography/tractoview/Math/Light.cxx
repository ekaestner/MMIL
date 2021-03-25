#include <assert.h>
#include "Light.h"

Light::Light() : mLightnumber( 0 ),
                mAmbient( 0.1f, 0.1f, 0.1f, 0.1f ),
                mDiffuse( 0.1f, 0.1f, 0.1f, 0.1f ),
                mSpecular( 0.1f, 0.1f, 0.1f, 0.1f ),
                mConstantAttenuation( 1.0f ),
                mLinearAttenuation( 0.0f ),
                mQuadraticAttenuation( 0.0f ),
                mIsOn( true ),
                mCutoff( 180.0f ),
                mDirection( 0.0f, 0.0f, -1.0f ),
                mExponent( 0.0f )
{

}
     // Light::setColor accepts a token for the color attribute to set (-
     // Light::ambient, Light::diffuse, or Light::specular) and three floating point
     // values (r, g, and b) in the range [0.0 .. 1.0] defining values for the
     // red, green, and blue components of the indicated attribute of the light
     // source.  By default, the r, g, and b values are all 1.0.
void Light::setColor( const Light::LightType& which, const float& r, const float& g, const float& b )
{
   switch( which )
   {
   case ambient: mAmbient.set( r, g, b ); break;
   case diffuse: mDiffuse.set( r, g, b ); break;
   case specular: mSpecular.set( r, g, b ); break;
   default:
   assert( which != ambient && which != diffuse && which != specular && "invalid LightType" );
   }
}

void Light::getColor( const LightType& which, float& r, float& g, float& b ) const
{
   switch( which )
   {
   case ambient: mAmbient.get( r, g, b ); break;
   case diffuse: mDiffuse.get( r, g, b ); break;
   case specular: mSpecular.get( r, g, b ); break;
   default:
   assert( which != ambient && which != diffuse && which != specular && "invalid LightType" );
   }
}
// Light::setAtten sets the attenuation parameters of the pfLight. The
// light intensity is scaled at each vertex by:
//
//      1.0 / (constant + linear * dist + quadratic * dist^2)
//
// where 'dist' is the distance from the light position to the lit vertex.
// Note that 'dist' is 1.0 for infinite light sources. The default
// attenuation values are constant = 1.0, linear = 0.0, quadratic = 0.0,
// i.e., light attenuation is disabled.
void Light::setAtten( const float& constant, const float& linear, const float& quadratic )
{
   mConstantAttenuation = constant;
   mLinearAttenuation = linear;
   mQuadraticAttenuation = quadratic;
}

// Light::getAtten returns the
//     attenuation parameters of the pfLight in constant, linear, and quadratic.
void Light::getAtten( float& constant, float& linear, float& quadratic ) const
{
   constant = mConstantAttenuation;
   linear = mLinearAttenuation;
   quadratic = mQuadraticAttenuation;
}

// Light::setPos receives four floating point values to set the x, y, z,
// and w, coordinates for the position of the light source.  Typically, the
// homogeneous coordinate w is 0.0 to indicate that the light position is
// infinitely far from the origin in the direction (x, y, z).  Local light
// sources are specified by a non-zero value for w and usually incur a
// performance penalty.
void Light::setPos( const float& x, const float& y, const float& z, const float& w )
{
   if (w == 0.0f)
   {
      Vec3<float> temp( x, y, z );
      temp.normalize();
      mPosition.set( temp[0], temp[1], temp[2], 0.0f );
   }
   else
   {
      mPosition.set( x, y, z, w );
   }
}

// Light::getPos copies the x, y, z and w
// coordinates of the light source into the parameters x, y, z and w,
// respectively.
void Light::getPos( float& x, float& y, float& z, float& w ) const
{
   mPosition.get( x, y, z, w );
}


// Light::setSpotDir specifies the direction in which a spot light source
// emits its light.  It receives three floating point values, x, y, and z,
// specifying the x, y, and z direction vectors.
void Light::setSpotDir( const float& x, const float& y, const float& z )
{
   mDirection.set( x, y, z );
}

     // Light::getSpotDir copies the x, y, and z direction vectors into
     // the parameters x, y, and z.
     void Light::getSpotDir( float& x, float& y, float& z ) const
{
   mDirection.get( x, y, z );
}

     // Light::setSpotCone specifies the exponent and spread of the spot light
     // cone, and receives two floating point values, f1 and f2, to set the
     // exponent for the intensity, and the spread of the cone, respectively.
     void Light::setSpotCone( const float& exponent, const float& spread )
{
if ((spread >= 0.0f && spread <= 90.0f) || spread == 180.0f)
   mCutoff = spread;

if (spread >= 0.0f && spread <= 128.0f)
   mExponent = exponent;
}

     // Light::getSpotCone copies the current exponent and spread of the cone
     // into the parameters f1 and f2.
     void Light::getSpotCone( float& exponent, float& spread ) const
{
spread = mCutoff;
exponent = mExponent;
}

     void Light::on()
{
     mIsOn = true;
}
     void Light::off()
{
mIsOn = false;
}
     bool Light::isOn() const
{
return mIsOn;
}
