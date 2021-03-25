#ifndef COPLANAR_FUNCS_INCLUDED
#define COPLANAR_FUNCS_INCLUDED

#include "Plane.h"
#include "Polygon.h"

namespace kev
{
   template <class dataType>
   inline bool isCoplanar( const Plane<dataType>& plane1, const Plane<dataType>& plane2, dataType tolerance = 0.0f )
   {
      bool result = plane1.normal().isEqual( plane2.normal(), tolerance ) &&
         kev::isEqual( plane1.constant(), plane2.constant(), tolerance );
      return result;
   }

   template <class dataType>
   inline bool isCoplanarOpposite( const Plane<dataType>& plane1, const Plane<dataType>& plane2, dataType tolerance = 0.0f )
   {
      Plane<dataType> flipped = plane1;
      flipped.flip();
      return isCoplanar( flipped, plane2, tolerance );
   }

   template <class dataType>
   inline bool isCoplanar( const Polygon<dataType> polygon, const Plane<dataType>& plane, const dataType& tolerance = 0.0f )
   {
      Plane<dataType> polygonPlane;
      polygon.getPlane( polygonPlane );

      return isCoplanar( polygonPlane, plane, tolerance );
   }

   template <class dataType>
   inline bool isCoplanarOpposite( const Polygon<dataType> polygon, const Plane<dataType>& plane, const dataType& tolerance = 0.0f )
   {
      Plane<dataType> polygonPlane;
      polygon.getPlane( polygonPlane );

      return isCoplanarOpposite( polygonPlane, plane, tolerance );
   }

}; //namespace kev



#endif


