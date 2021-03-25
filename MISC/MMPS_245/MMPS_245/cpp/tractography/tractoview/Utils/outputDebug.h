#ifndef OUTPUTDEBUG
#define OUTPUTDEBUG

namespace kev
{   
   void outputDebug( const Plane<float>& plane, const char* const name = "" )
   {
      cout<<name<<":  "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   }

   void outputDebug( const Vec3<float>& vec, const char* const name = "" )
   {
      cout<<name<<":  "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<"\n"<<flush;
   }
};

#endif
