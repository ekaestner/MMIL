#ifndef DATA_REFERENCE_MANAGEMENT
#define DATA_REFERENCE_MANAGEMENT

#include <assert.h>

// NOTE This is autopointer-like... inheritance will not be good here...

template< class _dataType >
class Ref
{
public:
   Ref() : mDataPointer( NULL )
   {
      this->newref();
   }
   
   virtual ~Ref() 
   {
      assert( mDataPointer != NULL && mDataPointer->mRefCount > 0 && "in deleting reference to data, found that somehow that reference was already destroyed" );
      this->deref();
   }
   // copy constructor.
   // adds (or becomes an additional) refence to data in 'r'
   Ref( const Ref& r ) : mDataPointer( NULL )
   {
      this->ref( r );
   }
   _dataType* operator->()
   {
      assert( mDataPointer != NULL && mDataPointer->mRefCount > 0 && "operator->(): the data you are accessing is unrefferenced, or this ref object has gone out of scope." );
      return &(mDataPointer->mData);
   }

   const _dataType* const operator->() const
   {
      return &(mDataPointer->mData);
   }
   
   _dataType* operator*()
   {
      return &(mDataPointer->mData);
   }

   const _dataType* const operator*() const
   {
      return &(mDataPointer->mData);
   }
   
   Ref& operator=( const Ref& r )
   {
      this->ref( r );
      return *this;
   }
   
   int count()
   {
      if (mDataPointer == NULL)
      {
         return 0;
      }
      else
      {
         return mDataPointer->mRefCount;
      }
   }
private:
   struct DataPointer
   {
      DataPointer() : mRefCount( 0 ), mData()
      {
      }
      _dataType mData;
      int mRefCount;
   };
   DataPointer* mDataPointer;
   void deref()
   {
      if (mDataPointer != NULL)
      {
         --mDataPointer->mRefCount;
         // if we're the last to deref this memory, 
         // then delete it.
         if (mDataPointer->mRefCount <= 0)
         {
            delete mDataPointer;
         }

         // since I am derefed, I reference nothing.
         // this keeps me from derefing more than once on the same data
         mDataPointer = NULL;
      }
   }
   // reference new memory of type _dataType
   void newref()
   {
      // im not sure it would hurt if this was called many times...
      // it may mess up the caller though if they for some reason
      // expected the refcount to go up, but why would they?
      //assert( mDataPointer == NULL && "newref() can only be called once" );
      if (mDataPointer == NULL)
      {
         mDataPointer = new DataPointer;
         ++mDataPointer->mRefCount;
      }
   }
   // make a reference to the same stuff ref is reffering to.
   void ref( const Ref& r )
   {
      // dereference any data that I was already referenceing
      this->deref();
      
      assert( mDataPointer == NULL && r.mDataPointer != NULL && r.mDataPointer->mRefCount > 0 && "you can only reference data that has reference count > 0" );
      mDataPointer = r.mDataPointer;
      ++( mDataPointer->mRefCount );
   }
};

#endif
