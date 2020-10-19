#ifndef KEV_RANDOM_H
#define KEV_RANDOM_H

#include <vector>
#include <stdlib.h> //rand
#include <math.h>


namespace kev
{
   // return a random number between 0.0f and 1.0f
   inline float getRandom()
   {
      // rand returns int from  0 - 2^15-1
      const float two15minone( pow( 2.0f, 15 ) - 1 );
      float r = static_cast<float>( rand() );
      r /= two15minone;
      return r;
   }

   // return a random number between x1 and x2
   inline float getRandom( float x1, float x2 )
   {
      float r = kev::getRandom();
      float size = x2 - x1;
      return r * size + x1;
   }

   template <class dataType>
   inline void randomize( std::vector<dataType>& items )
   {
      for (int x = items.size() - 1; x > 0; --x)
      {
         // get a random one to swap items[x] with
         int random = getRandom( 0, x );

         // do the swap
         dataType temp = items[x];
         items[x] = items[random];
         items[random] = temp;
      }      
   }
};

#endif
