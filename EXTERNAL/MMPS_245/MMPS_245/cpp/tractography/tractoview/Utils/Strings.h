#ifndef KEVIN_STRING_UTILITIES
#define KEVIN_STRING_UTILITIES

#include <string.h>

namespace kevnStringUtils
{
	// returns the last token in the a string
	inline char* strtokLast( char* string, char* after )
	{
		// see what the extension is.
		char* ext = ::strtok( string, after );
		char* temp = ext;
		
		while (temp != NULL)
		{
			temp = ::strtok( NULL, after );
			if (temp != NULL)
				ext = temp;
		}

		return ext;
	}

	// returns the last token in the string, or returns string if there is only one token.
	// doesn't modify the source string
	inline const char* const getLastStringToken( const char* const string, char delimiter )
	{
		for (int x = ::strlen(string) - 1; x > 0; --x)
		{
			if (string[x] == delimiter)
			{
				return &string[x];
			}
		}
		
		return string;
	}

	// returns just the .xxx part of a filename.
	// if there are no .'s in the name, then it returns a pointer to '\0' (the end)
	inline char* fileExtension( char* filename )
	{
		// see if there are no .'s in the name
		if ( ::strchr( filename, '.' ) == NULL )
			return &filename[ ::strlen(filename) - 1 ];;

		// return the extension after the last .
		return strtokLast( filename, "." );
	}
};
using namespace kevnStringUtils;

#endif
