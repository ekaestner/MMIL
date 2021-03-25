#ifndef TRIANGLE_INCLUDED
#define TRIANGLE_INCLUDED

#include "Vec3.h"
#include "Tri3.h"
#include "Plane.h"

namespace kev
{

   template <class dataType>
   class Tri3
   {
   public:
      const Vec3<dataType>& p0() const { return a; }
      const Vec3<dataType>& p1() const { return b; }
      const Vec3<dataType>& p2() const { return c; }

      void p0( dataType x, dataType y, dataType z ) { a.set( x, y, z ); }
      void p1( dataType x, dataType y, dataType z ) { b.set( x, y, z ); }
      void p2( dataType x, dataType y, dataType z ) { c.set( x, y, z ); }

      // edge 0
      Vec3<dataType> e0() const { return a - b; }
      Vec3<dataType> e1() const { return b - c; }
      Vec3<dataType> e2() const { return c - a; }

      inline void getNormal( Vec3<dataType>& normal ) const
      {
         normal = (a - b).cross(a - c);
         normal.normalize();
      }
      inline void getPlane( kev::Plane<dataType>& p ) const
      {
         Vec3<dataType> triNorm, triPoint;
         this->getNormal( triNorm );
         p.set( triNorm, a );
      }
   private:
      Vec3<dataType> a, b, c;
   };

};

#endif


