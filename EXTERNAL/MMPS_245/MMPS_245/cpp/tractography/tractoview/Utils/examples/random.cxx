#include <vector>
#include "../random.h"

void main()
{
   cout<<"0-1: "<<kev::getRandom( 0, 1 )<<"\n"
         <<"0-5: "<<kev::getRandom( 0, 5 )<<"\n"
         <<"0-100: "<<kev::getRandom( 0, 100 )<<"\n"
         <<"50-1000: "<<kev::getRandom( 50, 1000 )<<"\n"
         <<"50-51: "<<kev::getRandom( 50, 51 )<<"\n"
         <<"\n"<<flush;
   
   std::vector<float> vec;
   vec.resize( 10 );
   for (int x=0; x < vec.size(); ++x)
      vec[x] = x;
   
   kev::randomize( vec );
   
   for (x=0; x < vec.size(); ++x)
      cout<<x<<": "<<vec[x]<<"\n"<<flush;
}
