#ifndef SEGMENT_INCLUDED
#define SEGMENT_INCLUDED

template <class dataType>
class Seg3
{
public:
   Seg3() : mPoint0(), mPoint1() {}
   Seg3( const Vec3<dataType>& p0, const Vec3<dataType>& p1 ) : mPoint0( p0 ), mPoint1( p1 ) {}
   void set( const Vec3<dataType>& p0, const Vec3<dataType>& p1 )
   {
      mPoint0 = p0;
      mPoint1 = p1;
   }
   inline const Vec3<dataType>& p0() const { return mPoint0; }
   inline const Vec3<dataType>& p1() const { return mPoint1; }

   dataType length() const { return (mPoint0 - mPoint1).length(); }
private:
   Vec3<dataType> mPoint0, mPoint1;
};


#endif


