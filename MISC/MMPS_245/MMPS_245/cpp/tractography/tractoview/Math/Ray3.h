#ifndef RAY_INCLUDED
#define RAY_INCLUDED

template <class dataType>
class Ray3
{
public:
   Ray3() : mOrigin(), mDirection() {}
   Ray3( const Vec3<dataType>& origin, const Vec3<dataType>& direction ) : mOrigin( origin ), mDirection( direction )
   {
      mDirection.normalize();
   }
   void set( const Vec3<dataType>& o, const Vec3<dataType>& d )
   {
      mOrigin = o;
      mDirection = d;
      mDirection.normalize();
   }
   inline const Vec3<dataType>& origin() const { return mOrigin; }
   inline const Vec3<dataType>& direction() const { return mDirection; }
private:
   Vec3<dataType> mOrigin, mDirection; // O, D
};



#endif
