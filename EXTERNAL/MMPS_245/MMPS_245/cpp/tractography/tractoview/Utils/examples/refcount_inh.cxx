
#include "refobj.h"
#include <iostream.h>
#include <list>

class bogusData : public refobj
{
public:
   bogusData() : refobj()
   {  
   }
   
protected:
   virtual ~bogusData() {}

private:
   char data[1024];
};

class mg : public refobj
{
public:
   virtual ~mg() 
   { 
      int x = 0;
      int  notinv = 0;
      std::list<bogusData*>::iterator it;
      for (it = mBd.begin(); it != mBd.end(); )
      {
         (*it)->deref();
         if ((*it)->refIsInvalid() == false)
         {
            notinv++;;
         }         
         std::list<bogusData*>::iterator dit = it;
         ++it;
         mBd.erase( dit );
         ++x;
      }
      cout<<"erased and derefed "<<x<<" elmnts from the list\n"<<flush;
      if (notinv == true)
         cout<<notinv<<" elts are still valid, someone else must deref them\n"<<flush;
   }
   void add( bogusData* bd )
   {
      mBd.push_back( bd );
      bd->ref();
   }
   
   std::list<bogusData*> mBd;
};

void main()
{
   bogusData* me;
   
   me = new bogusData;
   cout<<me->refCount()<<" == 1\n"<<flush;
   me->ref();
   cout<<me->refCount()<<" == 2\n"<<flush;
   me->deref();
   cout<<me->refCount()<<" == 1\n"<<flush;
   me->deref();
      
   int x;
   for ( x = 0; x < 128; ++x)
   {
      me = new bogusData;
      me->deref();
   }
   
   mg* mglist = new mg;
   for ( x = 0; x < 128; ++x)
   {
      me = new bogusData;
      
      // add the data to the list
      mglist->add( me );
      
      // dereference my copy, because now, I am not responsible for it.
      me->deref();
   }
  
   cout<<"\n\nlist will deallocate, and will have 0 elts still valid\n"<<flush;
   mglist->deref();
   
   // this should not compile...
   //delete me;
   
   cout<<"\n\nlist will deallocate, but will have 1 elt still valid\n"<<flush;
   mglist = new mg;
      me = new bogusData;
      me->ref();
      mglist->add( me );
   mglist->deref();
      me->deref();
      
    cout<<"\n\nThis should crash or assert:\n"<<flush;
      me = new bogusData;
      
      // dereference my copy, because now, I am not responsible for it.
      me->deref();
      
      // add the data to the list
      mglist->add( me );
}
