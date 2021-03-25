#ifndef SPLITTING_FUNCTIONS
#define SPLITTING_FUNCTIONS

#include <vector>
#include <assert.h>

#include "Vec3.h"
#include "Tri3.h"
#include "Ray3.h"
#include "Seg3.h"
#include "Plane.h"
#include "Polygon.h"
#include "coplanar.h"
#include "outputDebug.h"

#include "intersect.h"

namespace kev
{
   const float DEFAULT_SPLIT_TOLERANCE = 0.0001f;
   // we have 9 cases for d:
   //
   //    d0 d1
   //   -------
   // 1| 0  0
   // 2| -  -
   // 3| +  +
   // 4| -  +
   // 5| +  -
   // 6| 0  -
   // 7| 0  +
   // 8| +  0
   // 9| -  0
   //
   // args: given a segment, get the distance each point is from the plane, 
   //       these dists are d0 and d1
   int splitclassify( float d0, float d1 )
   {
      // find which case we have
      if (d0 == 0 && d1 == 0) return 1;
      else if (d0 < 0 && d1 < 0) return 2;
      else if (d0 > 0 && d1 > 0) return 3;
      else if (d0 < 0 && d1 > 0) return 4;
      else if (d0 > 0 && d1 < 0) return 5;
      else if (d0 == 0 && d1 < 0) return 6;
      else if (d0 == 0 && d1 > 0) return 7;
      else if (d0 > 0 && d1 == 0) return 8;
      else if (d0 < 0 && d1 == 0) return 9;
      else 
      {
         assert( false && "this point should be unreachable, how did you get here?" );
         return -9999;
      }
   }

   // we have 9 cases for d:
   //
   //    d0 d1 d2
   //   ----------
   // 11| 0  0  0
   // 12| 0  0  +
   // 13| 0  0  -
   // 81| +  0  0
   // 82| +  0  +
   // 83| +  0  -
   // 91| -  0  0
   // 92| -  0  -
   // 93| -  0  +
   //
   // this function is only for sub classification of cases 6,7,8,9 
   // where a segment begins or ends on the plane (d1=0)
   // for a segment in case 6 or 7, d0 is the previous vertex distance
   // for a segment in case 8 or 9, d2 is the following vertex distance
   int splitsubclassify( float d0, float d1, float d2 )
   {
      assert( d1 == 0 && "this function is only for sub classification of cases 6,7,8,9 where a segment begins or ends on the plane (d1=0)" );
   
      // find which case we have
      if (d0 > 0 && d1 == 0 && d2 == 0) return 81;
      else if (d0 > 0 && d1 == 0 && d2 > 0) return 82;
      else if (d0 > 0 && d1 == 0 && d2 < 0) return 83;
      else if (d0 < 0 && d1 == 0 && d2 == 0) return 91;
      else if (d0 < 0 && d1 == 0 && d2 < 0) return 92;
      else if (d0 < 0 && d1 == 0 && d2 > 0) return 93;
      else if (d0 == 0 && d1 == 0 && d2 == 0) return 11; // colinear verts!!!
      else if (d0 == 0 && d1 == 0 && d2 > 0) return 12; 
      else if (d0 == 0 && d1 == 0 && d2 < 0) return 13; 
      else 
      { 
         assert( false && "this point should be unreachable, how did you get here?" ); 
         return -9999;
      }
   }


   // returns indecies for x-1, x+0, x+1, and x+2
   void getVectorNeighbors( const int& size, const int& x, int& x_minus_one, int& x_0, int& x_plus_one, int& x_plus_two )
   {
      x_minus_one = (x - 1 + size) % size;
      x_0         =  x             % size;
      x_plus_one  = (x + 1)        % size;
      x_plus_two  = (x + 2)        % size;
   }

   void getSubDistNext( const Plane<float>& plane, const Vec3<float>& vnext, float& d0, float& d1, float& e0, float& e1, float& e2, float tolerance = DEFAULT_SPLIT_TOLERANCE)
   {
      e0 = d0;
      e1 = d1;
      e2 = plane.getDistance( vnext );
      if (kev::isEqual( e2, 0.0f, tolerance ))
            e2 = 0.0f;
   }
   void getSubDistPrev( const Plane<float>& plane, const Vec3<float>& vprev, float& d0, float& d1, float& e0, float& e1, float& e2, float tolerance = DEFAULT_SPLIT_TOLERANCE)
   {
      e0 = plane.getDistance( vprev );
      if (kev::isEqual( e0, 0.0f, tolerance ))
            e0 = 0.0f;
      e1 = d0;
      e2 = d1;
   }
   #include <iostream.h>
   // last vertex must be a duplicate of the first
   // coincident verts are not allowed.
   // collinear are allowed.
   // from: Graphics Gems III p 219 "partitioning a 3d convex polygon with an arbitrary plane"
   // polygon's GState may be referenced after this function is done...
   void planeSplitsPolygon( const Plane<float>& plane, 
                           Polygon<float>& polygon, 
                           Polygon<float>& positive, 
                           Polygon<float>& negative,
                           float tolerance = DEFAULT_SPLIT_TOLERANCE )
   {
      positive.clear();
      negative.clear();
      int polysize = polygon.size();
      assert( polygon.size() >= 3 && "not a triangle" );
      assert( polygon.detectIdenticalVerts() == false && "the polygon has identical verticies, cannot continue" );
      assert( !kev::isCoplanar( polygon, plane, 0.0f ) && "test for coplanar before trying to split");
   

      for ( int x = 0; x < polygon.size(); ++x)
      {
         // get indicies for x-1, x+0, x+1, and x+2
         int x_minus_one, x_current, x_plus_one, x_plus_two;
         getVectorNeighbors( polygon.size(), x, x_minus_one, x_current, x_plus_one, x_plus_two );

         const Vertex<float>& v0 = polygon[ x_current ];
         const Vertex<float>& v1 = polygon[ x_plus_one ];
      
         // get distances
         float d0 = plane.getDistance( v0.coord() );
         float d1 = plane.getDistance( v1.coord() );

         if (kev::isEqual( d0, 0.0f, tolerance ))
            d0 = 0.0f;

         if (kev::isEqual( d1, 0.0f, tolerance ))
            d1 = 0.0f;

         // may or may not be used later...
         const Vertex<float>& vprev = polygon[ x_minus_one ];
         const Vertex<float>& vnext = polygon[ x_plus_two ];
         float       dSub0, dSub1, dSub2;
      
         int classification = splitclassify( d0, d1 );
         switch (classification)
         {
         // The edge is embedded in the plane, since +00- or -00+ cannot 
         // happen (convex polygons only), we know that the edge is either on 
         // the positive side or the negative side.  Use subclassify to find out which:
                     // v0 v1
                     //------
         case 1:     // 0  0
            getSubDistNext( plane, vnext.coord(), d0, d1, dSub0, dSub1, dSub2 );
            switch (splitsubclassify( dSub0, dSub1, dSub2 ))
            {
                     // v0 v1 vn          v0 is...
                     //------ --          ------------
            case 12: // 0  0   +          -- positive
               positive.push_back( v0 );
               break;
            case 13: // 0  0   -          -- negative
               negative.push_back( v0 );
               break;
            case 11: // 0  0   0          -- crap...   argghhhh!  colinear verts!!!  what to do?!?!?
               positive.push_back( v0 );     // just do this, (the dumb thing), and later,
               negative.push_back( v0 );     // when pos or neg polygon is [size < 3] then clear() it.
                                             //: TODO: the only time you have an edge on the plane is when the rest of 
                                             //  the points are to one side so we can just put v0 on the same side that 
                                             //  most of the verts lay...
                                             //  TODO: here, determine that side, then assign v0 to it...
                                             // really just do a preprocess for collinear, and strip them out
                                             // since they aren't useful to the polygon's shape.
               cout<<"WARNING: detected colinear verts\n"<<flush;
               break;
            }
            break;

         // The edge is totally on the negative, or positive side of the plane.  
         // so there is no intersection
                     // v0 v1             v0 is...
                     //------             ------------
         case 2:     // -  -              -- negative
            negative.push_back( v0 ); 
            break;
         case 3:     // +  +              -- positive
            positive.push_back( v0 ); 
            break;


         // The edge crosses the plane, so there is an intersection, compute it, 
         // splitting the seg into each half...
         // v0 inserted once, vI inserted twice.  v1 will be added next time around...
                     // v0 v1             v0 is...
                     //------             ------------
         case 4:     // -  +              -- negative plus we add an intersected point vI
         case 5:     // +  -              -- positive plus we add an intersected point vI
            {
            Vertex<float> vI;
            bool        result( false );
            intersectSegPlane( v0, v1, plane, vI, result );
            if (result == false)
            {
               cout<<"v0:    "<<v0.coord()[0]<<" "<<v0.coord()[1]<<" "<<v0.coord()[2]<<"\n"<<flush;
               cout<<"v1:    "<<v1.coord()[0]<<" "<<v1.coord()[1]<<" "<<v1.coord()[2]<<"\n"<<flush;
               cout<<"plane: "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
               cout<<"vI:    "<<vI.coord()[0]<<" "<<vI.coord()[1]<<" "<<vI.coord()[2]<<"\n"<<flush;
            }
            assert( result && "intersection should be true by definition of case 4 or 5" );
      
            if (classification == 4) 
               negative.push_back( v0 );
            else if (classification == 5)
               positive.push_back( v0 );
            else 
            {
               assert( false && "this point should be unreachable, how did you get here?" );
            }

            negative.push_back( vI );
            positive.push_back( vI );
            }
            break;
         
         // the edge (v0,v1) is touching the plane at its start vertex, v0
         // to determine we need to know the classification of the [prev] vertex before the edge: v(-1)
         // case 6 is handled by case 83 (intersection) or 92 (no intersection)
         // case 7 is handled by case 93 (intersection) or 82 (no intersection) 
            //: we need to classify [d(-1) d0 d1]
            //  get the previous (v-1) distance == d0 (swapping all values to fit in d0 d1 d2...)
                     //    v0 v1
                     //    -----
         case 6:     //    0  -
            getSubDistPrev( plane, vprev.coord(), d0, d1, dSub0, dSub1, dSub2 );
            switch (splitsubclassify( dSub0, dSub1, dSub2 ))
            {
                     // vp v0 v1          v0 is...
                     // -- -----          -----------
            case 83: // +  0  -           -- both (intersected)
               negative.push_back( v0 );
               positive.push_back( v0 );
               break;
            case 13: // 0  0  -           -- negative
            case 92: // -  0  -           -- negative
               negative.push_back( v0 );
               break;
            default: assert( false && "this point should be unreachable, how did you get here?" );
            }
            break;
         case 7:     //    0  +
            getSubDistPrev( plane, vprev.coord(), d0, d1, dSub0, dSub1, dSub2 );
            switch (splitsubclassify( dSub0, dSub1, dSub2 ))
            {
                     // vp v0 v1          v0 is...
                     // -- -----          -----------
            case 12: // 0  0  +           -- positive
            case 82: // +  0  +           -- positive
               positive.push_back( v0 );
               break;
            case 93: // -  0  +           -- both (intersected)
               positive.push_back( v0 );
               negative.push_back( v0 );
               break;
            default: assert( false && "this point should be unreachable, how did you get here?" );
            }
            break;

         // the edge (v0,v1) is touching the plane at its ending vertex, v1
         // to determine we need to know the classification of the [next] vertex after the edge: v2
         // case 8 is handled by 81 82 83
         // case 9 is handled by 91 92 93
         // only case 83 and 93 indicate intersection
         case 8:    // +  0
            getSubDistNext( plane, vnext.coord(), d0, d1, dSub0, dSub1, dSub2 );
            switch (splitsubclassify( dSub0, dSub1, dSub2 ))
            {
                    // v0 v1 vn                v0 is...
                    // ----- --             -----------
            case 81:// +  0   0             -- positive
            case 82:// +  0   +             -- positive
            case 83:// +  0   -             -- positive (next will be 6/83, an intersection)
               positive.push_back( v0 ); 
               break;
            default: assert( false && "this point should be unreachable, how did you get here?" );
            }
            break;

         case 9:    // -  0
            //: we need to sub classify the distances using the next distance "d2"
            //  get the next (v2) distance == dSub2
            getSubDistNext( plane, vnext.coord(), d0, d1, dSub0, dSub1, dSub2 );
            switch (splitsubclassify( dSub0, dSub1, dSub2 ))
            {
                    // v0 v1 vn                v0 is...
                    // ----- --             -----------
            case 91:// -  0   0             -- negative
            case 92:// -  0   -             -- negative
            case 93:// -  0   +             -- negative (next will be 7/93, an intersection)
               negative.push_back( v0 );
               break;
            default: assert( false && "this point should be unreachable, how did you get here?" );
            }
            break;

         // unreachable (if all was coded well)
         default: 
            assert( false && "this point should be unreachable, how did you get here?" );
            break;
         }
      }

      int wantedNumVerts = 3;
      int polygonSize = polygon.size();
      int negSize = negative.size();
      int posSize = positive.size();
      if (
          !((positive.size() >= wantedNumVerts || positive.size() == 0) && (negative.size() >= wantedNumVerts || negative.size() == 0))
          || 
          (positive.size() == 0 && negative.size() == 0)
         )
      {
         cout<<"plane:\n"<<flush;
         outputDebug( plane ); cout<<"\n"<<flush;
      
         cout<<"poly:\n"<<flush;
         polygon.outputDebug(); cout<<"\n"<<flush;

         cout<<"negative:\n"<<flush;
         negative.outputDebug();cout<<"\n"<<flush;
      
         cout<<"positive:\n"<<flush;
         positive.outputDebug();cout<<"\n"<<flush;
      }

      
      if (positive.size() < wantedNumVerts)
         positive.clear();
      if (negative.size() < wantedNumVerts)
         negative.clear();
      
      assert( ((positive.size() >= wantedNumVerts || positive.size() == 0) && (negative.size() >= wantedNumVerts || negative.size() == 0)) && "both polys are not even triangles, this is wrong, because the input polygon was at least a triangle (pre)" );

      int negSize2 = negative.size();
      int posSize2 = positive.size();

      assert( positive.detectIdenticalVerts() == false && "positive polygon has identical verticies, something is wrong with the splitting algorithm - try changing the tolerance a little higher" );
      assert( negative.detectIdenticalVerts() == false && "negative polygon has identical verticies, something is wrong with the splitting algorithm - try changing the tolerance a little higher" );
      assert( positive.size() >= wantedNumVerts || positive.size() == 0 && "positive tri is not a triangle" );
      assert( negative.size() >= wantedNumVerts || negative.size() == 0 && "negative tri is not a triangle" );
      assert( ((positive.size() >= wantedNumVerts || positive.size() == 0) && (negative.size() >= wantedNumVerts || negative.size() == 0)) && "both polys are not even triangles, this is wrong, because the input polygon was at least a triangle (post)" );
      assert( !(positive.size() == 0 && negative.size() == 0) && "splitting a polygon should give at least 1 polygon" );
   

      // make both positive and negative inherit polygons's GState
      positive.refGState( polygon.gstate() );
      negative.refGState( polygon.gstate() );
   }


}; //namespace kev

#endif
