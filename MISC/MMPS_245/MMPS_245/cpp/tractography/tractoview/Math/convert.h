#ifndef CONVERSION_ROUTINES
#define CONVERSION_ROUTINES

template<class dataType>
inline void seg2ray( const Seg3<dataType>& seg, Ray3<dataType>& ray )
{
   Vec3<dataType> origin( seg.p0() );
   Vec3<dataType> dir( seg.p1() - seg.p0() );
   dir.normalize();

   ray.set( origin, dir );
}

#endif
