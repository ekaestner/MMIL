#ifndef VERTEX_H
#define VERTEX_H

#include "Vec3.h"
#include "Vec2.h"
#include "Vec4.h"

namespace kev
{
      template <class dataType>
      class Vertex
      {
      public:
         Vertex() : mCoord(), mTexCoord(), mColor(), mNormal()
         {
         }
         
         Vertex( dataType x, dataType y, dataType z ) : mCoord(x, y, z), mTexCoord(), mColor(), mNormal()
         {
         }
         
         void outputDebug() const
         {
            cout<<"["<<mCoord[0]<<" "<<mCoord[1]<<" "<<mCoord[2]<<"], "<<flush;
         }

         inline void setCoord( const Vec3<dataType>& coord ) { mCoord = coord; }
         inline Vec3<dataType>& coord() { return mCoord; }
         inline const Vec3<dataType>& coord() const { return mCoord; }
         
         inline void setTexCoord( const Vec2<dataType>& texcoord ) { mTexCoord = texcoord; }
         inline Vec2<dataType>& texcoord() { return mTexCoord; }
         inline const Vec2<dataType>& texcoord() const { return mTexCoord; }
         
         inline void setColor( const Vec4<dataType>& color ) { mColor = color; }
         inline Vec4<dataType>& color() { return mColor; }
         inline const Vec4<dataType>& color() const { return mColor; }
         
         inline void setNormal( const Vec3<dataType>& normal ) { mNormal = normal; }
         inline Vec3<dataType>& normal() { return mNormal; }
         inline const Vec3<dataType>& normal() const { return mNormal; }
         

         Vertex( const Vertex& vertex ) : mCoord(vertex.coord()), mTexCoord(vertex.texcoord()), mColor(vertex.color()), mNormal(vertex.normal())
         {
         }
         
         Vertex& operator=( const Vertex& vertex )
         {
            this->setCoord( vertex.coord() );
            this->setColor( vertex.color() );
            this->setNormal( vertex.normal() );
            this->setTexCoord( vertex.texcoord() );
            return *this;
         }

                  
      private:
         Vec3<dataType> mCoord;
         Vec2<dataType> mTexCoord;
         Vec4<dataType> mColor;
         Vec3<dataType> mNormal;
      };
};


#endif
