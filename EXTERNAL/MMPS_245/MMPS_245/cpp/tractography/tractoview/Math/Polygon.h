#ifndef POLYGON_INCLUDED
#define POLYGON_INCLUDED

#include "Vertex.h"
#include "Vec3.h"
#include "Vec2.h"
#include "Vec4.h"
#include "Plane.h"
#include <vector>
#include <assert.h>

#include "GState.h"

namespace kev
{
   template <class dataType>
   class Polygon
   {
   public:
      Polygon() : mGState( NULL ) {}
      inline void set( const Polygon& polygon ) 
      {
         this->refGState( polygon.mGState );
         mVerts = polygon.mVerts;
      }
      inline Polygon( const Polygon& polygon ) : mGState( NULL )
      {
         this->set( polygon );
      }
      inline Polygon& operator=( const Polygon& polygon ) 
      {
         this->set( polygon );
         return *this;
      }
      ~Polygon() 
      { 
         if (mGState != NULL) 
            mGState->deref();
         mGState = NULL;
      }
      inline void refGState( GState* gstate )
      { 
         // if was already referenceing a GState, then unreference it.
         if (mGState != NULL)
         {
            mGState->deref(); 
            mGState = NULL;
         }
         
         // set the gstate
         mGState = gstate; 
         
         // if the newly set GState is non-NULL, then reference it
         if (mGState != NULL) 
         {
            mGState->ref(); 
         }
      }
      inline GState* gstate() { return mGState; }
      inline const GState* gstate() const { return mGState; }
      inline Vertex<dataType>&       operator[]( int x )       { assert( x < mVerts.size() ); return mVerts[x]; }
      inline const Vertex<dataType>& operator[]( int x ) const { assert( x < mVerts.size() ); return mVerts[x]; }
      inline int                     size() const              { return mVerts.size(); }
      inline void                    resize( int x )           { mVerts.resize(); }
      inline void                    clear() { mVerts.clear(); }
      inline void                    push_back( const Vertex<dataType>& vert ) { mVerts.push_back( vert ); }
      
      inline Vec3<dataType>&         coord( int x )       { assert( x < mVerts.size() ); return mVerts[x].coord(); }
      inline const Vec3<dataType>&   coord( int x ) const { assert( x < mVerts.size() ); return mVerts[x].coord(); }
      inline Vec3<dataType>&         normal( int x )       { assert( x < mVerts.size() ); return mVerts[x].normal; }
      inline const Vec3<dataType>&   normal( int x ) const { assert( x < mVerts.size() ); return mVerts[x].normal; }
      inline Vec4<dataType>&         color( int x )       { assert( x < mVerts.size() ); return mVerts[x].coord(); }
      inline const Vec4<dataType>&   color( int x ) const { assert( x < mVerts.size() ); return mVerts[x].coord(); }
      inline Vec2<dataType>&         texcoord( int x )       { assert( x < mVerts.size() ); return mVerts[x].texcoord; }
      inline const Vec2<dataType>&   texcoord( int x ) const { assert( x < mVerts.size() ); return mVerts[x].texcoord; }
      
      inline void getNormal( Vec3<dataType>& normal ) const
      {
         // needs to at least be a triangle
         assert( this->size() >= 3 );
      
         normal = (this->coord(2) - this->coord(1)).cross(this->coord(0) - this->coord(1));
         normal.normalize();
      }
      inline void getPlane( Plane<dataType>& plane ) const
      {
         // needs to at least be a triangle
         int size = this->size();
         assert( size >= 3 );
      
         Vec3<dataType> norm;
         this->getNormal( norm );
         plane.set( norm, this->coord(0) );
      }

      bool detectIdenticalVerts() const
      {
         int x;
         
         for (x = 0; x < mVerts.size(); ++x)
         {
            Vec3<float> x1 = mVerts[x].coord(), 
                        x2 = mVerts[(x + 1) % mVerts.size()].coord();
            float l = (x1 - x2).length();
            if (x1 == x2)
               return true;
         }

         return false;
      }

      void fixIdenticalVerts()
      {
         int x;
         std::vector<Vertex<dataType> > temp;

         for (x = 0; x < mVerts.size(); ++x)
         {
            Vertex<float> &x1 = mVerts[x], 
                        &x2 = mVerts[(x + 1) % mVerts.size()];
            if (x1.coord() != x2.coord())
               temp.push_back( x1 );
         }

         mVerts.clear();
         mVerts = temp;
      }

      bool detectCollinear() const
      {
         int x;
         
         if (this->detectIdenticalVerts() == true)
            return true;

         for (x = 0; x < mVerts.size(); ++x)
         {
            Vec3<float> x1 = mVerts[((x + mVerts.size()) - 1) % mVerts.size()].coord(), 
                        x2 = mVerts[x].coord(), 
                        x3 = mVerts[(x + 1) % mVerts.size()].coord();
            Vec3<float> v1(x1 - x2),
                        v2(x2 - x3);
            v1.normalize();
            v2.normalize();
            float costheta = v1.dot(v2);
            if (costheta == 1.0f || costheta == -1.0f)
            {
               return true;
            }
         }

         return false;
      }
      
      void fixCollinear()
      {
         std::vector<Vertex<dataType> > temp2;

         this->fixIdenticalVerts();
         int x;

         for (x = 0; x < mVerts.size(); ++x)
         {
            Vertex<float> x1 = mVerts[((x + mVerts.size()) - 1) % mVerts.size()], 
                        x2 = mVerts[x], 
                        x3 = mVerts[(x + 1) % mVerts.size()];
            Vec3<float> v1(x1.coord() - x2.coord()),
                        v2(x2.coord() - x3.coord());
            v1.normalize();
            v2.normalize();
            float costheta = v1.dot(v2);
            if (costheta == 1.0f || costheta == -1.0f)
            {
            }
            else
            {
               temp2.push_back( x2 );
            }
         }

         mVerts.clear();
         mVerts = temp2;
      }

      void outputDebug() const 
      {
         std::vector<Vertex<dataType> >::const_iterator it;
         for (it = mVerts.begin(); it != mVerts.end(); ++it)
         {
            (*it).outputDebug();
         }
      }
      
      void normalize()
      {
         Vec3<float> normal;
         this->getNormal( normal );
         int x;
         for (x = 0; x < mVerts.size(); ++x)
         {
            mVerts[x].setNormal( normal );
         }         
      }

   public:
      // CCW ordering...
      std::vector<Vertex<dataType> > mVerts;
      GState* mGState;
   };
}; //namespace kev

#endif
