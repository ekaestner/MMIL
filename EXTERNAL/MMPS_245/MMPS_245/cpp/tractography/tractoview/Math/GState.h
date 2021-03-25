#ifndef GSTATE_INCLUDED
#define GSTATE_INCLUDED

#include "Material.h"
#include "Texture.h"
#include <string>
#include "refobj.h"

namespace kev
{
   class GState : public refobj
   {
   public:
      GState()
      {
      }
      std::string name;
      Material mat;

      std::string mapName;
      
      mutable Texture texture;
   };
};

#endif
