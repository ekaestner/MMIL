
///////////////////////////////////////////////////////////////
// 
#include <Defines.h>
#include <Matrix4f.h>
#include <kstack.h>
#include "LsystemFunctor.h"

#include "Lsystem.h"

kstack<Matrix4f> Lsystem::matrixStack;

Lsystem::Lsystem()
{
    memset( _functorTable, NULL, sizeof(_functorTable) );
}
    
Lsystem::~Lsystem()
{
    
}

//set the current string that describes the look of the plant.
void Lsystem::setString(const char* string)
{
    _currentString = string;
    //cout<<"L-System: "<<_currentString<<"\n"<<flush;
}

//in the current string:
// replaces every instance of "source" with "dest"
void Lsystem::searchAndReplace(const char* source, const char* dest)
{
	std::string newString;
    int count = 0;
    
    std::string sourceString = source;
    std::string replacementString = dest;
    
    //cout<<"Source: "<<sourceString<<"\n"<<flush;
    //cout<<"Replacement: "<<replacementString<<"\n"<<flush;
    
    // The lsystem wont be degenerating today ...
    assert( sourceString.size() <= replacementString.size() );
    
    //serious performance hits will occur if you don't 
    // pre-allocate the memory for the new string.
    newString.reserve(_currentString.size()*2);
    
    // search and replace:
    // build a new string, don't modify the old one.  
    // then copy the new one into the old one.
    for(int x = 0; x < _currentString.size(); ++x)
    { 
	//copy the chars into the new string.
	newString += _currentString[x];
	
	// if the two chars don't match, then reset
	// we need to begin matching again.
	if( _currentString[x] != sourceString[count] )
	    count = 0;
	
	// if the two chars do match, increment the count
	// so we can check the next two chars.
	if( _currentString[x] == sourceString[count] )
	    ++count;
	
	// if "count" is the size of the source string
	// then we've matched a whole substring in the 
	// currentString with the "sourceString".
	// Modify and add to newString, the "replaceString".
	if( count == sourceString.size() )
	{
	    //reset "count"
	    count = 0;
	    
	    //Find the beginning of the recognized source text.
	    // reverse "sourceString.size()" number of chars 
	    // in newString
	    int sourceBegin = newString.size() - sourceString.size();
	    
	    // start replacing there, with "replaceString"
	    for(int y = 0; y < replacementString.size(); ++y)
	    {
		if(y < sourceString.size() )
		    newString[sourceBegin + y] = replacementString[y];
		else
		    newString += replacementString[y];
	    }
	}
    }
    
    // replace the current string with the new string.
    _currentString = newString;
}

//render the current lsystem.
void Lsystem::render()
{
    Lsystem::matrixStack.clear();
    Matrix4f currentMatrix( Matrix4f::identity() );
    Lsystem::matrixStack.push( currentMatrix );
    
    // fill the vector with each char of the string
    for(int x = 0; x < _currentString.size(); ++x)
    { 
	LsystemFunctor* functor = this->getFunctor( _currentString[x] );
	if( functor != NULL )
	    functor->render();
    }
}

//add a function to the system.
// specify functor object that derives from LsystemFunctor.
// if one was previously set for that character, then 
// addRule will return a pointer to this functor.  
// You should delete it if you created it with "new"
void Lsystem::addFunctor(const char& character, 
			LsystemFunctor*& newFunctor, 
			LsystemFunctor*& oldFunctor)
{
    oldFunctor = NULL;
    
    // we're replacing a functor, return the old one to the client
    if( _functorTable[ static_cast<int>(character) ] )
	oldFunctor = _functorTable[ static_cast<int>(character) ];
    
    // set the new functor
    _functorTable[ static_cast<int>(character) ] = newFunctor;
}

//get a character-associated functor.
LsystemFunctor* Lsystem::getFunctor(const char& character)
{
    return _functorTable[ static_cast<int>(character) ];
}



