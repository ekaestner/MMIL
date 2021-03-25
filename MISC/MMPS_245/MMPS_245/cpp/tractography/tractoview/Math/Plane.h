#ifndef PLANE_INCLUDED
#define PLANE_INCLUDED

#include <Utils/Defines.h> // defines isEqual( float float ... );
#include "Vec3.h"
#include "Plane.h"

namespace kev
{

   template <class dataType>
   class Plane
   {
   public:
      Plane();
      Plane( const Vec3<dataType>& normal, const dataType& constant );
      Plane( const Vec3<dataType>& normal, const Vec3<dataType>& pointOnPlane );
      Plane( const Plane<dataType>& plane );

      // d = Ax + By + Cz + D
      dataType getDistance( const Vec3<dataType>& p ) const;

      void getNormal( Vec3<dataType>& normal ) const;
      void getConstant( dataType& constant ) const;

      const Vec3<dataType>& normal() const;
      const dataType& constant() const;

      // define a plane with a point and a normal
      void set( const Vec3<dataType>& normal, const Vec3<dataType>& pointOnPlane );

      // define a plane with a normal and the offset from origin to plane (shortest distance)
      void set( const Vec3<dataType>& normal, const dataType& constant );
      void set( const Plane<dataType>& plane );
      void flip();

   private:
      // plane equation:
      // D = Ax + By + Cz
      Vec3<dataType> mNormal;   // A, B, C
      dataType       mConstant; // D
   };

   template <class dataType>
   inline Plane<dataType>::Plane() : mNormal(), mConstant(0) {}

   template <class dataType>
   inline Plane<dataType>::Plane( const Vec3<dataType>& normal, const dataType& constant ) : mNormal( normal ), mConstant( constant ) 
   {
      mNormal.normalize();
   }
      
   
   template <class dataType>
   inline Plane<dataType>::Plane( const Vec3<dataType>& normal, const Vec3<dataType>& pointOnPlane ) : mNormal(), mConstant(0) 
   {
      this->set( normal, pointOnPlane );
   }

   template <class dataType>
   inline Plane<dataType>::Plane( const Plane<dataType>& plane ) : mNormal(plane.normal()), mConstant(plane.constant())
   {
   }

   // d = Ax + By + Cz - D
   template <class dataType>
   inline dataType Plane<dataType>::getDistance( const Vec3<dataType>& p ) const
   {
      return p.dot( mNormal ) - mConstant;
   }

   template <class dataType>
   inline void Plane<dataType>::getNormal( Vec3<dataType>& normal ) const 
   { 
      normal = mNormal; 
   }

   template <class dataType>
   inline void Plane<dataType>::getConstant( dataType& constant ) const 
   { 
      constant = mConstant; 
   }

   template <class dataType>
   inline const Vec3<dataType>& Plane<dataType>::normal() const 
   { 
      return mNormal; 
   }

   template <class dataType>
   inline const dataType& Plane<dataType>::constant() const 
   { 
      return mConstant;
   }

   // define a plane with a point and a normal
   template <class dataType>
   inline void Plane<dataType>::set( const Vec3<dataType>& normal, const Vec3<dataType>& pointOnPlane )
   {
      mNormal = normal;
      mNormal.normalize();

      // plane equation: distance from origin = Normal dot any point on the plane
      //   D = (A*x+B*y+C*z)    which also == |N||p|cos(t)
      //   D =  N.p
      mConstant = mNormal.dot( pointOnPlane );
   }

   // define a plane with a point and a normal
   template <class dataType>
   inline void Plane<dataType>::set( const Vec3<dataType>& normal, const dataType& constant )
   {
      mNormal = normal;
      mNormal.normalize();
      mConstant = constant;
   }

   template <class dataType>
   inline void Plane<dataType>::set( const Plane<dataType>& plane )
   {
      mNormal = plane.normal();
      mConstant = plane.constant();
   }

   template <class dataType>
   inline void Plane<dataType>::flip()
   {
      mNormal = -mNormal;
      mConstant = -mConstant;
   }

}; //namespace kev

#endif
