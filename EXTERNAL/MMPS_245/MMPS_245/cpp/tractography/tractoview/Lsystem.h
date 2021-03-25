
////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////
#ifndef LSYSTEM_INCLUDED
#define LSYSTEM_INCLUDED

#include <string>
#include <kstack.h>
#include <Matrix4f.h>

class LsystemFunctor;
class Lsystem
{
public:
    Lsystem();
    virtual ~Lsystem();
    
    //set the current string that describes the look of the plant.
    //the characters are set up with "addFunctor"
    void setString(const char* string);
    
    //in the current string:
    // replaces every instance of "source" with "dest"
    void searchAndReplace(const char* source, const char* dest);

    //render the current lsystem.
    virtual void render();
    
protected:
    //add a function to the system.
    // specify functor object that derives from LsystemFunctor.
    // if one was previously set for that character, then 
    // addRule will return a pointer to this functor.  
    // You should delete it if you created it with "new"
    void addFunctor(const char& character, 
		    LsystemFunctor*& newFunctor, 
		    LsystemFunctor*& oldFunctor);
    
    //get a character-associated functor.
    LsystemFunctor* getFunctor(const char& character);
    
    
protected:
    //this works for 8bit ASCII text only, wont work for 16bit UNICODE
    LsystemFunctor* _functorTable[256];
    
	std::string _currentString;
    
public:
    friend class LsystemFunctor;
	static kstack<Matrix4f> matrixStack;
};

#endif
