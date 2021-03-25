#include <math.h>
#include <iostream.h>

#include "Endian.h" // for endian funcs.

void main()
{
	if (kev::littleEndian() == true)
		cout<<"Little Endian\n"<<flush;
	
	if (kev::littleEndian() == false)
		cout<<"Big Endian\n"<<flush;
}
