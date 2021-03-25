#ifndef CAMERA_INCLUDED_H
#define CAMERA_INCLUDED_H

#include "ViewVolumeGeometry.h"
#include "Matrix4f.h"

// Camera (as defined here) is a location and a frustum (two space transforms)
// it takes objects from world space -> camera space -> view space
// (view space is what is then transformed to your screen)
// NOTE: don't confuse cameraspace and the class Camera - they are very different.
//       The class Camera just does a transform similar to a RL camera that we are all used to.
class Camera
{
public:
   Camera() : mFrustumGeometryInCameraspace(), mCameraToProjection( Matrix4f::identity() ), 
                     mProjectionToCamera( Matrix4f::identity() ), 
                     mWorldToCamera( Matrix4f::identity() ), 
                     mCameraToWorld( Matrix4f::identity() ),
                     mWorldToProjection(Matrix4f::identity() ),
                     mFrustumGeomNeedsUpdate( true )
   {
   }
public:
   // set the cameraspace->projectionspace transform
   // this is your projection matrix
   inline void setProjectionXform( const Matrix4f& cameraToProjectionXform ) 
   { 
      mCameraToProjection = cameraToProjectionXform;
      mProjectionToCamera = mCameraToProjection;
      mProjectionToCamera.invert();
      
      // FG is already unit cube (projectionspace), put it into cameraspace
      mFrustumGeometryInCameraspace.make( mProjectionToCamera ); 
      this->update(); 
   }
	
   // set the worldspace->cameraspace transform
   // this is your camera/navigation matrix
   inline void setCameraXform( const Matrix4f& worldToCameraXform ) 
   { 
      mWorldToCamera = worldToCameraXform;
      mCameraToWorld = mWorldToCamera;
      mCameraToWorld.invert();
      this->update();
   }

   // set the camera to some position/orientation in world space
   inline void setPosition( const Matrix4f& position ) 
   { 
      mCameraToWorld = position;
      mWorldToCamera = mCameraToWorld;
      mWorldToCamera.invert();
      this->update();
   }

   // get the geometry for the view volume in world coordinate space.
   inline ViewVolumeGeometry& frustumGeometry() 
   { 
      if (mFrustumGeomNeedsUpdate == true)
      {
         mFrustumGeometryInWorldspace = mFrustumGeometryInCameraspace;
         mFrustumGeometryInWorldspace.xform( mCameraToWorld );
         mFrustumGeomNeedsUpdate = false;
      }

      return mFrustumGeometryInWorldspace;
   }

public:
   // the Camera xform (takes worldspace to view space)
   inline const Matrix4f& matrix() const { return mWorldToProjection; }

   // just the location (takes worldspace to cameraspace)
   inline const Matrix4f& cameraXform() const { return mWorldToCamera; }

   // just the frustum (takes cameraspace to view space)
   inline const Matrix4f& projectionXform() const { return mCameraToProjection; }

   // camera's position in world coordinates. (takes cameraspace to worldspace)
   // for example, this could be used to get the camera's position for rendering camera model in other viewports
   // or comparing the eye point to objects in the scene for back to front ordering.
   inline const Matrix4f& position() const { return mCameraToWorld; }

private:
   void update() { Matrix4f::multiply( mWorldToProjection, mCameraToProjection, mWorldToCamera ); mFrustumGeomNeedsUpdate = true; }
   Matrix4f mWorldToCamera;
   Matrix4f mCameraToWorld; // inv

   Matrix4f mCameraToProjection;
   Matrix4f mProjectionToCamera; // inv
   
   Matrix4f mWorldToProjection;

   ViewVolumeGeometry mFrustumGeometryInCameraspace, mFrustumGeometryInWorldspace;
   bool mFrustumGeomNeedsUpdate;
};


#endif

