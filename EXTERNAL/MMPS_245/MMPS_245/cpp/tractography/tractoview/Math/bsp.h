#ifndef BSP_TREE_INCLUDED
#define BSP_TREE_INCLUDED

#include <list>
#include <Utils/random.h>
#include <Math/Polygon.h>
#include <Math/split.h>


namespace kev
{
   class BspNode
   {
   public:
      BspNode() : mNegSide( NULL ), mPosSide( NULL ), mPlane(), mSameBucket(), mOppBucket()
      {
         assert( mNegSide == NULL && mPosSide == NULL );
      }
      ~BspNode() 
      {
         if (mNegSide != NULL) 
         {
            delete mNegSide; 
            mNegSide = NULL;
         }
         if (mPosSide != NULL) 
         {
            delete mPosSide; 
            mPosSide = NULL;
         }
      }

      inline float whichSideIsViewer( const Vec3<float>& viewPosition ) const
      {
          return mPlane.getDistance( viewPosition );
      }
      
      inline float whichSideIsViewFrustum( const std::vector<Vec3<float> >& vf ) const
      {
         int neg( 0 ), pos( 0 );

         for (int x = 0; x < vf.size(); ++x)
         {
            float d = this->whichSideIsViewer( vf[x] );
            if (d < 0.0f)         // negative
            {
               ++neg;
            }
            else                  // positive
            {
               ++pos;
            }
         }

         if (neg >= 1 && pos == 0)
            return -1;
         else if (neg == 0 && pos >= 1)
            return 1;
         else 
            return 0;
      }   
      
      void constructTree( std::list<Polygon<float> >& polygons )
      {
         if (polygons.size() == 0)
            return;

         //: a bsp node is a plane subdividing a scene of polygons...
         // it holds its coplanar polygons, and pointers to pos and neg sides of the BSP tree

         // set my plane
         this->choosePlane( polygons, mPlane );

         // partition a list of polygons by the plane into 4 buckets:
         // - buckets of polygons in negative/positive sides of the plane
         // - buckets of coplanar polygons facing in same/opposite directions as plane
         std::list<Polygon<float> > neg_bucket, pos_bucket;
         this->partitionPolygonsWithPlane( mPlane, polygons, neg_bucket, pos_bucket, mSameBucket, mOppBucket );
         
         // recursively process remaining Polygons, if any, on either side
         if (neg_bucket.size() != 0)
         {
            mNegSide = new BspNode();
            mNegSide->constructTree( neg_bucket );
         }
         if (pos_bucket.size() != 0)
         {
            mPosSide = new BspNode();
            mPosSide->constructTree( pos_bucket );
         }
      }
   
      // using the given plane, organize the given list of polygons into 4 buckets: neg, pos, same, opp
      void partitionPolygonsWithPlane( const Plane<float>& plane, 
         std::list<Polygon<float> >& polygons, 
         std::list<Polygon<float> >& negativeBucket, 
         std::list<Polygon<float> >& positiveBucket, 
         std::list<Polygon<float> >& sameBucket, 
         std::list<Polygon<float> >& oppositeBucket, const float& tolerance = 0 )
      {
         std::list<Polygon<float> >::iterator it;
         for (it = polygons.begin(); it != polygons.end(); ++it)
         {
            Polygon<float>& polygon = (*it);
            Plane<float> polygonPlane;
            polygon.getPlane( polygonPlane );

            if (kev::isCoplanar( polygonPlane, plane, tolerance ))
            {
               sameBucket.push_back( polygon );
            }
            else if (kev::isCoplanarOpposite( polygonPlane, plane, tolerance ))
            {
               oppositeBucket.push_back( polygon );
            }
            else
            {
               Polygon<float> positivePolygon, negativePolygon;
               planeSplitsPolygon( plane, polygon, positivePolygon, negativePolygon );

               if (positivePolygon.size() != 0)
                  positiveBucket.push_back( positivePolygon );
            
               if (negativePolygon.size() != 0)
                  negativeBucket.push_back( negativePolygon );
            }
         }
      }

      // function to best choose a plane in a group of Polygon<float>s
      void choosePlane( const std::list<Polygon<float> >& polys, Plane<float>& plane )
      {
         int size = polys.size();
         
         // this number could be changed to some heuristic based on examination of the polygon data
         int which = kev::getRandom( 0, size - 1 ); 
         
         std::list<Polygon<float> >::const_iterator it = polys.begin();
         int x = 0;
         for (; it != polys.end() && x < which; ++it, ++x)
         {
            // advancing the iterator to the random position...
         }
         (*it).getPlane( plane );
      }

      Plane<float>& plane() { return mPlane; }
      const BspNode* negative() const { return mNegSide; }
      const BspNode* positive() const { return mPosSide; }
      std::list<Polygon<float> >&  same() { return mSameBucket; }
      std::list<Polygon<float> >&  opposite() { return mOppBucket; }

      const std::list<Polygon<float> >&  same() const { return mSameBucket; }
      const std::list<Polygon<float> >&  opposite() const { return mOppBucket; }

      
   private:
      Plane<float> mPlane;
      std::list<Polygon<float> > mSameBucket, mOppBucket;
      BspNode* mNegSide;
      BspNode* mPosSide;
   };

}; //namespace kev

#endif
