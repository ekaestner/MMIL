#ifndef SPRITE_XFORM
#define SPRITE_XFORM

#include "Matrix4f.h"
#include "Camera.h"

// transform for sprites (billboards).
// draw your sprites in the XY plane, and the SpriteTransform will do the rest.
class SpriteTransform
{
public:
   SpriteTransform() : mSpriteMatrix( Matrix4f::identity() ), mPosition( 0.0f, 0.0f, 0.0f ), mAxis( 0.0f, 0.0f, 0.0f )
   {
   }

   // must call this after every camera change, but before rendering
   inline void rotate( const Camera& camera )
   {
      if (mAxis == Vec3<float>( 0.0f, 0.0f, 0.0f ))
      {
         mSpriteMatrix.makeTranslation( mPosition[0], mPosition[1], mPosition[2] );
         mSpriteMatrix.copyRotation( camera.position() );
      }
      else
      {
         assert( false && "axis sprites not implemented" );
      }
   }

   // set the position of your sprite.
   void setPosition( const Vec3<float>& position )
   {
      mPosition = position;
   }

   // set the axis for the sprite to rotate about (good for trees for example)
   // set to 0,0,0 for a point sprite (good for smoke for example)
   void setAxis( const Vec3<float>& axis )
   {
      mAxis = axis;
   }

   const Matrix4f& matrix() const { return mSpriteMatrix; }

private:
   Vec3<float> mPosition, mAxis;
   Matrix4f mSpriteMatrix;
};



#endif

