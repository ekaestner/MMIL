#include "Defines.h"
#include "Vec3.h"
#include "Tri3.h"
#include "Ray3.h"
#include "Seg3.h"
#include "Plane.h"

#include "convert.h"

namespace kev
{

   void intersectRayPlane( const Ray3<float>& ray, 
					   const Plane<float>& plane, 
					   Vec3<float>& intersectedPoint, bool& result )
   {
      const Vec3<float> &N = plane.normal(), &O = ray.origin(), &D = ray.direction();
      const float &d = plane.constant();

      // if plane and ray are parallel (N.D = 0) the intersection is rejected
      float nDotD = N.dot(D);
      if (nDotD == 0)
      {
         result = false;
         return;
      }

      // 1.) use the plane equation:
      //     N.P = d
      //
      // 2.) and the ray equation:
      //     r(t) = O + D*t
      //
      // 3.) evaluation of t corresponding to the intersection point can be 
      //     obtained using eq 1 and 2
      //             d + N.O
      //     t = - -----------
      //               N.D
      float t = -(N.dot(O) - d) / nDotD;

      // if the intersection is behind the origin of the ray (t < 0) the intersection is rejected
      // TODO: is there any reason i should instead do a "behind or equal"? (t <= 0)
      if (t < 0.0f)
      {
         result = false;
         return;
      }
   
      // 4.) using the parametric ray equation r(t) = O + D*t
      //     find the intersection point
      intersectedPoint = O + D * t;

      // intersection successful
      result = true;
   };



   void intersectSegPlane( const Seg3<float>& seg, 
					   const Plane<float>& plane, 
					   Vec3<float>& intersectedPoint, bool& result )
   {
      float d0 = plane.getDistance( seg.p0() );
      float d1 = plane.getDistance( seg.p1() );

      bool negPos = (d0 <= 0 && d1 >= 0);
      bool posNeg = (d0 >= 0 && d1 <= 0);
      if (!negPos && !posNeg)
      {
         result = false;
         //cout<<"both points are on the same side of the plane\n"<<flush;
         return;
      } // after this the seg could still be parallel d0==0 && d1==0

      // reuse the ray/plane intersect algorithm... (this will fail only if the seg is parallel)
      Ray3<float> ray;
      seg2ray( seg, ray );
      intersectRayPlane( ray, plane, intersectedPoint, result );
   };

   void intersectSegPlane( const Vertex<float>& v0, const Vertex<float>& v1,
                           const Plane<float>& plane,
                           Vertex<float>& intersectedPoint, bool& result )
   {
      Seg3<float> seg( v0.coord(), v1.coord() );
      kev::intersectSegPlane( seg, plane, intersectedPoint.coord(), result );
      
      if (result == true)
      {
         float length = (v1.coord() - v0.coord()).length();
         float lengthToIsect = (intersectedPoint.coord() - v0.coord()).length();
         //float lerp = (length - lengthToIsect) / length;
         float lerp = (lengthToIsect / length);

         // get the intersected texcoord
         kev::Lerp( v0.texcoord(), v1.texcoord(), lerp, intersectedPoint.texcoord() );

         // get the intersected normal
         kev::Lerp( v0.normal(), v1.normal(), lerp, intersectedPoint.normal() );

         // get the intersected color
         kev::Lerp( v0.color(), v1.color(), lerp, intersectedPoint.color() );
      }
   }


   void intersectSegTri( const Seg3<float>& seg, 
					   const Tri3<float>& tri, 
					   Vec3<float>& intersectedPoint, bool& result )
   {
      Plane<float> plane;
      tri.getPlane( plane );
      intersectSegPlane( seg, plane, intersectedPoint, result );

      // angleAtPoint(x) - the Tri3's angle at point(x)
       // iVector - the vector from point(x) to intersectedPoint
       // iAngle - the angle between edge(x) and iVector 
       // if iAngle is less (wider) than the angleAtPoint(x) 
       //    then point is not in the Tri3. 
    
    
       //: try at the (point0) location
    
      Vec3<float> edge0 = tri.p1() - tri.p0();
      Vec3<float> edge2 = tri.p2() - tri.p0();
      edge0.normalize();
      edge2.normalize();

       float angleAtPoint0( edge0.dot(edge2) );
    
       Vec3<float> iVector( intersectedPoint - tri.p0() );
       iVector.normalize();

       float iAngle( iVector.dot(edge0) );
    
       if (iAngle <= angleAtPoint0)
       {
         result = false;
         return;
       }	
    
    
       //: try at the (point1) location
   
       // we didn't need edge1 until now...
       Vec3<float> edge1 = tri.p2() - tri.p1();
       edge1.normalize();
    
       float angleAtPoint1( edge1.dot(-edge0) );
    
       iVector.set( intersectedPoint - tri.p1() );
       iVector.normalize();
    
       iAngle = iVector.dot( edge1 );
    
       if (iAngle <= angleAtPoint1)
       {
         result = false;
         return;
       }	
    
    
       //: try at the (point2) location
    
       float angleAtPoint2( edge1.dot(edge2) );
    
       iVector.set( intersectedPoint - tri.p2() );
       iVector.normalize();
    
       iAngle = iVector.dot(-edge2);
    
       if (iAngle <= angleAtPoint2)
       {
         result = false;
         return;
       }	
    
       //: if we made it here, then the point is inside the Tri3.
       result = true;
   }

}; // namespace kev
