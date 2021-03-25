
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
#ifndef CONTEXT_DATA_H
#define CONTEXT_DATA_H


//-----------------------------------------------------------------------
//: map: int --> <Type>
//  Lets you associate data of any type to an integer "context"
//
// Ex: Map a float
//   ContextData<float> myData;
//   myData( 0 ) = 3.14596f;
//   myData( 1 ) = 6.29082f;
//   for (int context = 0; context < 2; ++context)
//   {
//      cout << "Data in context " << context << " is: " << myData( context ) << "\n";
//   }
//
// Output:
//   Data in context 0 is: 3.14596f
//	 Data in context 1 is: 6.29082f
//
//! NOTE: Requires that the type of the context data provide a default
//+  constructor used to initialize all of the copies of the data.
//-----------------------------------------------------------------------
#include <vector>
template< class ContextDataType >
class ContextData
{
// client interface
public:
	//: Returns reference to user data for the contextId
    ContextDataType& operator() ( const int& contextId );
    
protected:
    //: Return a ptr to the correct data element in the current context
    //! PRE: We are in the draw function
    ContextDataType&  getItem( const int& contextId );

private:
	//: Vector of user data
    std::vector< ContextDataType >    _contextDataVector;
};

//: Return a ptr to the correct data element in the current context
//! PRE: We are in the draw function
template<class ContextDataType>
inline ContextDataType& ContextData<ContextDataType>::getItem( const int& contextId )
{   
    // Make sure that we will reference a valid element
    while(_contextDataVector.size() <= contextId)
    {
		_contextDataVector.push_back(ContextDataType());
    }
    
    return _contextDataVector[contextId];
}

//: Returns reference to user data for the contextId
template<class ContextDataType>
inline ContextDataType& ContextData<ContextDataType>::operator() ( const int& contextId )
{ 
    return this->getItem( contextId ); 
}

#endif
