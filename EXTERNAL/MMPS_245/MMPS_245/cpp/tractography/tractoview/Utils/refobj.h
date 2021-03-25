// this refobjer is better(general) than refptr.h because it allows polymorphism.
// i may choose this method, but Ref.h is better(safer), since it is smart pointer'ish
//

// this class is useful for inheritance, but is less automatic, 
// and will be harder to debug your problems resulting from misuse... since you are 
// trusted to call ref and deref correctly


#ifndef MANAGED_BY_REF_COUNT__INHERITANCE
#define MANAGED_BY_REF_COUNT__INHERITANCE
#include <assert.h>

// you MUST call deref to delete, do not call delete!!!!!
// TODO: hide new and delete somehow, since delete should never be used, and new makes it tempting to call delete.
class refobj
{
public:
    //: default constructor
    // count() == 1 after creation.
    refobj();

    refobj( const refobj& r ) : ___mNumTimesReferenced( 1 ), ___mBadcount( -69 )
    {
       //dont copy ref number.
    }

public:
    //: destructor
    // you MUST call deref to delete, do not call delete!!!!!
    virtual ~refobj()
    {
       assert( ___mNumTimesReferenced != ___mBadcount && "this data has already been deleted, it cannot be deleted twice (this is a weird bug)." );
       ___mNumTimesReferenced = ___mBadcount;
    }

public:
    refobj& operator=( const refobj& r )
    {
       //dont copy ref number.
       return *this;
    }

   //: increase the reference count
   // call this if you need
   void ref() const;

   //: decrease the reference count.
   void deref() const;
 
   //: get the reference count
   const int& refCount() const;

   bool refIsInvalid() { return ___mNumTimesReferenced == ___mBadcount; }

private:
    mutable int ___mNumTimesReferenced;
//   static const int ___mBadpointer;
   const int ___mBadcount;
};

//const int refobj::___mBadpointer( 0xdeadface );

/*
// use this class to create an object of type T, that is also referenceable.
template <class T>
class Countable : public refobj, public T
{
};
*/
          
//: default constructor
inline refobj::refobj() : ___mBadcount( -69 ), ___mNumTimesReferenced( 1 )
{
    
}

//: increase the reference count
inline void refobj::ref() const
{
   assert( ___mNumTimesReferenced != ___mBadcount && "this data has been deleted, someone is probably holding on to the data after they called deref()" );
   assert( ___mNumTimesReferenced != 0 && "this data should have been deleted, someone is probably holding on to the data after they called deref()\n" );
    ++___mNumTimesReferenced;
}

//: decrease the reference count.
inline void refobj::deref() const
{
   assert( ___mNumTimesReferenced != ___mBadcount && "this data has been dereferenced more times than referenced, someone is probably holding on to the data after they called deref(), or someone called deref too many times." );
   --___mNumTimesReferenced;
   if (___mNumTimesReferenced <= 0)
   {
      delete this;
   }
}

//: get the reference count
inline const int& refobj::refCount() const
{
   assert( ___mNumTimesReferenced != ___mBadcount && "this data has been deleted, someone is probably holding on to the data after they called deref()" );
   return ___mNumTimesReferenced;
}

#endif
