
#include <stdio.h>
#include <iostream.h>
#include "autoref.h"

class bogusData
{
public:
   bogusData()
   {
      sprintf( array, "I am the Walrus! aaaaaaaaaaaaaaaabbbbbbbbbbbbbbbccccccccccccccccccc\n" );
      array[1000] = 't';
      array[1001] = 'e';
      array[1002] = 's';
      array[1003] = 't';
   }
   char array[1024];
};

void functionTest1( autoref<bogusData> dataref )
{
   if (dataref.count() == 6)
      cout<<"Passed function pass-by-value, function has claimed a reference to the data\n"<<flush;
   else
      cout<<"Failed pass-by¯value\n"<<flush;

   // increase the count by a couple, just so we'll have a bunch die,
   // and thus deref back to 5.
   autoref<bogusData> thisWillDie0 = dataref;
   autoref<bogusData> thisWillDie1 = dataref;
   autoref<bogusData> thisWillDie2 = dataref;
   autoref<bogusData> thisWillDie3 = dataref;
   autoref<bogusData> thisWillDie4 = dataref;
   autoref<bogusData> thisWillDie5 = dataref;
   if (thisWillDie5.count() == 12 && thisWillDie4->array[1000] == 'b')
      cout<<"Passed making many data refs (this will never fail if previous tests passed)\n"<<flush;
}

void main()
{
   autoref<bogusData> bogus_data;
   if (bogus_data.count() == 1)
      cout<<"Passed add ref on creation of reference\n"<<flush;
   else
      cout<<"!! Failed\n";
   if (bogus_data->array[1000] == 't')
      cout<<"Passed data access of newly created reference to data\n"<<flush;
   else
      cout<<"!! Failed\n";
   
   autoref< bogusData > a_reference( bogus_data );
   if (bogus_data.count() == 2 && a_reference.count() == 2)
      cout<<"Passed copy constructor increaded refcount\n"<<flush;
   else
      cout<<"!! Failed\n";
   
   bogus_data->array[1000] = 'b';
   if (a_reference->array[1000] == 'b')
      cout<<"Passed data change in one reference, leads to data change in another reference (of same data)\n"<<flush;
   else
      cout<<"!! Failed\n";
      
   autoref< bogusData > b_reference = bogus_data;
   if (b_reference->array[1000] == 'b' && b_reference.count() == 3)
      cout<<"Passed make reference using (=)copy constructor, this check should be redundant of the copyconstructor test since they use the same funcs\n"<<flush;
   else
      cout<<"!! Failed\n";
   
   autoref< bogusData > c_reference;
   c_reference = bogus_data;
   if (c_reference->array[1000] == 'b' && c_reference.count() == 4)
      cout<<"Passed make reference using operator==(), data and count both checks out ok.\n"<<flush;
   else
      cout<<"!! Failed\n";

   autoref< bogusData > d_reference;
   d_reference = c_reference;
   if (d_reference->array[1000] == 'b' && d_reference.count() == 5)
      cout<<"Passed make reference using a 2nd generation reference.\n"<<flush;
   else
      cout<<"!! Failed\n";

   functionTest1( c_reference );

   if (c_reference->array[1000] == 'b' && d_reference.count() == 5)
      cout<<"Passed all references owned by function dereferenced ok.\n"<<flush;
   else
      cout<<"!! Failed\n";
}
