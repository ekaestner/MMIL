
//////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#ifndef STACK_INCLUDED
#define STACK_INCLUDED

#include <list>

template<class Data> 
class kstack
{
public:   
    //default constructor. 
    kstack();
    
    //destructor
    // NOTE: does not delete data
    ~kstack();
    
    // get a reference to the newest object on the stack.
    Data& top();
    
    // get a reference to the newest object on the stack.
    const Data& top() const;
    
    // remove the newest item from the stack. 
    // NOTE: does not 'delete' data
    void pop();
    
    // check to see if an item is present in the stack.
    bool isPresent(const Data& item);// const;
    
    // insert a new item into the stack.
    void push(const Data& item);
    
    //clear the stack.
    // NOTE: does not delete data
    void clear();
    
    //return the size.
    inline int size() const { return _stack.size(); }
    
protected:
    std::list<Data> _stack;
};


template<class Data>
inline kstack<Data>::kstack()
{
}

template<class Data>
inline kstack<Data>::~kstack()
{
}

// get a reference to the current object on the stack.
template<class Data>
inline Data& kstack<Data>::top()
{
    //assert( _stack.size() > 0 );
    return _stack.back();
}

// get a const reference to the current object on the stack.
template<class Data>
inline const Data& kstack<Data>::top() const
{
    //assert( _stack.size() > 0 );
    return _stack.back();
}

// pop an item off the stack.
template<class Data>
inline void kstack<Data>::pop()
{
    assert( _stack.size() > 0 );
    _stack.pop_back();
}

// push an item onto the stack.
template<class Data>
inline void kstack<Data>::push(const Data& item)
{
    _stack.push_back(item);
}

// check to see if an item is present in the stack.
template<class Data>
inline bool kstack<Data>::isPresent(const Data& item) //const
{
    //list<Data>::const_iterator itemIT;
	
    const Data itemIT;
    for(itemIT = _stack.begin(); itemIT != _stack.end(); ++itemIT)
    {
	//found a match
	if( item == *itemIT )
	    return true;
    }
    
    //couldn't find a match
    return false;
}

//clear the stack.
// NOTE: does not delete data
template<class Data>
inline void kstack<Data>::clear()
{
    _stack.empty();
}

#endif
