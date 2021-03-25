#ifndef VIEW_VOLUME_GEOMETRY
#define VIEW_VOLUME_GEOMETRY

#include "Matrix4f.h"

class ViewVolumeGeometry
{
public:
   ViewVolumeGeometry()
   {
      this->init();
   }
   
   void init()
   {
      mUnitCube.resize( 8 );
      mUnitCube[0].set( -1, -1, -1 );
      mUnitCube[1].set( -1, 1, -1 );
      mUnitCube[2].set( 1, 1, -1 );
      mUnitCube[3].set( 1, -1, -1 );
      
      mUnitCube[4].set( -1, -1, 1 );
      mUnitCube[5].set( -1, 1, 1 );
      mUnitCube[6].set( 1, 1, 1 );
      mUnitCube[7].set( 1, -1, 1 );
      
      mVerts = mUnitCube;
   }
   
   // set to the inverse of your view matrix (the one set by Matrix4f::makePerspective())
	void make( const Matrix4f& invViewMat )
	{
      for (int x = 0; x < 8; ++x)
      {
         mVerts[x] = invViewMat * mUnitCube[x];
      }
	}

   // this is the inverse of your camera matrix (the one set by Trackball, or Matrix4f::makeLookAt())
   // pass in a camera -> world transform here.
   void xform( const Matrix4f& invCameraMatrix )
   {
      for (int x = 0; x < 8; ++x)
      {
         mVerts[x] = invCameraMatrix * mVerts[x];
      }
   }

   Vec3<float> getEye() const
   {
      Vec3<float> viewerPosition( (mVerts[0] + mVerts[1] + mVerts[2] + mVerts[3]) * 0.25 );
      return viewerPosition;
   }   
   
   const std::vector< Vec3<float> >& verts() const { return mVerts; }
  
private:
   std::vector< Vec3<float> > mVerts;
   std::vector< Vec3<float> > mUnitCube;
};


#endif

