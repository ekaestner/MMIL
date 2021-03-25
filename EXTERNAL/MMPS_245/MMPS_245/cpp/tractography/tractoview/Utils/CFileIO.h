
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#ifndef CFILEIO_INCLUDED
#define CFILEIO_INCLUDED

#include "Endian.h" //needed for kev::bigEndian/kev::littleEndian funcs
#include <fstream.h>
#include <string>

namespace CFileIO
{
	//: true - 
	inline bool fileExists( const char* const name )
	{
		FILE* file = ::fopen( name, "r" );
		if (file == NULL)
		{
			return false;
		}
		
		else
		{
			::fclose( file );
			return true;
		}
	}
	
   enum Endianness
   {
      LITTLE, BIG
   };

	//: Read
	// this function is designed to work with little endian files, such as TGA and BMP.
   // this will return data that is correctly translated for your system architecture.
   // i.e. reading a little endian BMP on a big endian mips system, 
   //      returns (translated) bigendian data
	// i.e. reading a little endian TGA on a little endian intel system, 
   //      returns (untranslated) littleendian data
	template<class typeT>
	inline static int ReadData( Endianness fileByteOrdering, FILE* fp, typeT& data )
	{
		int size = ::fread( &data, sizeof(typeT), 1, fp );

		// if we're not on a little endian machine (intel is little endian) then reverse the bytes.
		if (fileByteOrdering == LITTLE && kev::bigEndian() ||
         fileByteOrdering == BIG && kev::littleEndian())
      {
         kev::byteReverse( data );
		}

		return size;
	}

	//: Write
	// this function is designed to work with little endian files, such as TGA and BMP.
	template<class typeT>
	inline static int WriteData( Endianness fileByteOrdering, FILE* fp, const typeT& data )
	{
		register typeT tempData = data;

		// if we're not on a little endian machine (i.e. intel is little endian, mips is big) 
      // then reverse the bytes.
		if (fileByteOrdering == LITTLE && kev::bigEndian() ||
         fileByteOrdering == BIG && kev::littleEndian())
      {
			kev::byteReverse( tempData );
		}
			
		int size = ::fwrite( &tempData, 1, sizeof(typeT), fp );

		return size;
	}



	// Read one byte of data from the file stream
	inline static int ReadByte( FILE *fp, unsigned char& value )
	{
		return ::fread( &value, 1, sizeof(unsigned char), fp );
	}

	// Read one short of data from the file stream
	// NOTE: these enhance readability and type checking, try to use them over the generic ReadData function
	inline static int ReadShort( Endianness fileByteOrdering, FILE* fp, unsigned short& value )
	{
		return ReadData( fileByteOrdering, fp, value );
	}

	// Read one long of data from the file stream
	// NOTE: these enhance readability and type checking, try to use them over the generic ReadData function
	inline static int ReadLong( Endianness fileByteOrdering, FILE *fp, unsigned long& value )
	{
		return ReadData( fileByteOrdering, fp, value );
	}

	// Write one byte of data to the file stream
	inline static int WriteByte( FILE *fp, const unsigned char& value )
	{
		return ::fwrite( &value, 1, sizeof(unsigned char), fp );
	}

	// Write one short of data to the file stream
	// NOTE: these enhance readability and type checking, try to use them over the generic WriteData function
	inline static int WriteShort( Endianness fileByteOrdering, FILE* fp, const unsigned short& value )
	{
		return WriteData( fileByteOrdering, fp, value );
	}

	// Write one long of data to the file stream
	// NOTE: these enhance readability and type checking, try to use them over the generic WriteData function
	inline static int WriteLong( Endianness fileByteOrdering, FILE *fp, const unsigned long& value )
	{
		return WriteData( fileByteOrdering, fp, value );
	}

   inline static void getLine( ifstream& f, std::string& text  )
   {
      char buffer[2049];
      f.getline( buffer, 2048, '\n' );
      buffer[2048] = '\0';
      if (f.gcount() < 2048)
      {
         buffer[f.gcount()] = '\0';
      }
      text = buffer;
   }

   inline static void getAll( ifstream& f, std::string& buffer )
   {
	   //cout<<"Reading:["<<flush;
      while ((f.eof() == false) && (f.fail() == 0) )
      //while (f.gcount() >= 2048)
      {
	      //cout<<"."<<flush;
         char buf[2049];
         f.read( buf, 2048 );
         buf[2048] = '\0';
         if (f.gcount() < 2048)
         {
            buf[f.gcount()] = '\0';
         }
         buffer += buf;
      }
	   //cout<<"]\n"<<flush;
      //cout << "Gcount == " << f.gcount() << "\n"<<flush;
   }
};

#endif
