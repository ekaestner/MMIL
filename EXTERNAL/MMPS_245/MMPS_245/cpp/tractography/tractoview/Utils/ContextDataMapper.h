
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////

#include "ContextData.h"

//: a class to map data to contexts. (an example explains this best... read on)
//
//  usage:
//     // allocate a new mapper for your data type (int)
//     ContextDataMapper<int> mapper;
//
//     // when in the nth context, write to the ith data:
//     mapper(n, i) = (int)3;
//
//     // get a new ID to reference data with:
//     int id = genLookupId(1);
// 
//     // get your data from slot 'id' for the nth context:
//     int data = mapper(n, id);
//
template< class Data = int >
class ContextDataMapper
{
public:
    ContextDataMapper() : _nextSlotId(0) {}
      
    //: generate a slot number to put your data into
    // NOTE: you can call this sizeof(unsigned int) times for now.
    unsigned int getNextOpenSlot(const unsigned int& numberSlotsRequested);
    
    //: look up a display list id associated with a "Renderer" display list ID.
    //  give - the current context you're in.
    //  give - the slot number you wish to access data from.
    //  returns - a graphics API (OpenGL) display list ID number.
    Data& operator() ( const int& contextId, const unsigned int& lookupId );
    
private:
    // I use ContextData class for auto-growing vectors.
    ContextData< ContextData<Data> > _contextDataMapping;
    
    // the next slot available
    unsigned int _nextSlotId;
};

//: generate a slot number to put your data into
// NOTE: you can call this sizeof(unsigned int) times for now.
template< class Data >
inline unsigned int ContextDataMapper<Data>::getNextOpenSlot(const unsigned int& numberSlotsRequested)
{
    int nextAvailableSlot = _nextSlotId;
    
    // calc the next available slot.
    _nextSlotId += numberSlotsRequested;
    
    return nextAvailableSlot;
}

//: look up a display list id associated with a "Renderer" display list ID.
//  give - the current context you're in.
//  give - the slot number you wish to access data from.
//  returns - a graphics API (OpenGL) display list ID number.
template< class Data >
inline Data& ContextDataMapper<Data>::operator() ( const int& contextId, const unsigned int& lookupId )
{
    // make sure we know what the "next" ID will always be.
    if (_nextSlotId <= lookupId)
	_nextSlotId = lookupId + 1;
	
    return (_contextDataMapping( contextId ))( lookupId );
}

