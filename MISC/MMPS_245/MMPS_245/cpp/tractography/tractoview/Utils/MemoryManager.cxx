#define MEMORY_MANAGER
#include <malloc.h>  // malloc
#include <crtdbg.h>



#include <string.h>  // strtok
#include <process.h> // exit
#include <stdio.h>   // sprintf
#include <assert.h>  // assert

//: Assert
//  if expression "verify" is false, then "text" is logged, and code asserts only in DEBUG builds
#if defined(_DEBUG) || defined(DEBUG)
	void Assert( const bool& verify, const char* const text )
	{
		if (verify == false)						
		{											
			#ifdef WIN32
				_RPT1( _CRT_WARN, "ASSERT: %s\n", text );
			#else
				//cout<<"ASSERT: "<<text<<"\n"<<flush;
			#endif
		}											
		assert( verify );
	}
#else
	#define Assert( verify, text ) ((void)0)
#endif

#define Fatal Assert

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

// Set to 1 to test unloading all objects.  System will run very slow, but it should work.
#define UNLOAD_ALL 0

#if defined(_DEBUG) || defined(DEBUG)
	// Memory header used for debug builds.
	struct MemoryHeader 
	{
		char			file[48];				// filename
		int				line;					// line number
		long			magic;					// magic number
		MemoryHeader*	next;					// pointer to next in list
		MemoryHeader*	previous;				// pointer to previous in list
		int				size;					// # of unsigned chars allocated, excluding debug info
	};

	//
	// Memory trailer used for debug builds.
	//
	struct MemoryTrailer {
		long magic;
	};
	const long MEMORY_MAGIC_NUMBER = 0x12345678;	// magic number for memory trailer

	//
	// Linked list for reporting memory leaks.
	//
	struct MemoryHeader *memoryLinkedList = NULL;		// list of allocated objects

#else
	//
	// Memory header used for release builds.
	//
	struct MemoryHeader {
		unsigned int  size;							// # of unsigned chars allocated, excluding debug info
	};
#endif


static unsigned long volatile memoryInUse = 0;		// total memory currently allocated

#if defined(_DEBUG) || defined(DEBUG)
	//----------------------------------------------------------------------------
	// new... Global operator new.  (DEBUG version)
	//
	void *operator new(unsigned int sz, const char *file, int line)
	{
		// Allocate memory.
		MemoryHeader *header = (MemoryHeader *)malloc(sz + sizeof(MemoryHeader) + sizeof(MemoryTrailer));
		if (header == NULL) 
		{
			Assert(0, "");
			exit(1);
		}


		const char* filename = NULL;
		if (::strlen(file) > sizeof(header->file))
		{
			filename = getLastStringToken( file, '\\' );
		}
		else filename = file;

		// Fill in header & trailer.
		if (filename) 
			::strncpy(header->file, filename, sizeof(header->file));
		else     
			header->file[0] = '\0';
		header->magic = MEMORY_MAGIC_NUMBER;
		header->line = line;
		header->size = sz;
		MemoryTrailer *trailer = (MemoryTrailer *)((unsigned char *)header + sizeof(MemoryHeader) + sz);
		trailer->magic = MEMORY_MAGIC_NUMBER;

		// Link to head of list.
		header->next = memoryLinkedList;
		header->previous = NULL;
		memoryLinkedList     = header;
		if (header->next) 
			header->next->previous = header;

		// Update global size.
		memoryInUse += sz;

		// Return pointer to data.
		return ((void *)((unsigned char *)header + sizeof(MemoryHeader)));
	}

	// standard operator new (overloaded for our wonderful memory manager)
	void*	operator new( unsigned int sz ) 
	{
		return ::operator new(sz, "(unknown -- didn't use the MemoryManager)", 0);
	}

	// delete... Global operator delete.  (DEBUG version)
	void operator delete(void *m)
	{
		if (m) 
		{
			MemoryHeader *header  = (MemoryHeader *)((unsigned char *)m - sizeof(MemoryHeader));

			// Ensure magic numbers are correct.
			MemoryTrailer *trailer = (MemoryTrailer *)((unsigned char *)header + sizeof(MemoryHeader) + header->size);
			if ((header->magic != MEMORY_MAGIC_NUMBER) || (trailer->magic != MEMORY_MAGIC_NUMBER)) 
			{
				Assert( false, "Attempted to delete memory not alocated by the MemoryManager" );
			}

			else
			{
				// Unlink from list.
				if (header->previous) 
					header->previous->next = header->next;
				else
					memoryLinkedList = header->next;
				
				if (header->next) 
					header->next->previous = header->previous;

				// Free memory.
				memoryInUse -= header->size;

				::memset( header, 0x69696969, header->size );
				free( header );
			}
		}
	}

	

	//----------------------------------------------------------------------------
	// ReportMemoryLeaks... Report any memory leaks.
	//
	void ReportMemoryLeaks()
	{
		if (memoryLinkedList == NULL) {
			_RPT0(_CRT_WARN, "MemoryManager: No memory leaks detected\n");
			return;
		}
		long sz = 0;
		_RPT0(_CRT_WARN, "\n\nMemoryManager: Memory leaks detected (some of these could be globals which will die after the MM.):\n");
		_RPT0(_CRT_WARN, "               Call ReportMemoryLeaks from your global destructors for more up to date info.\n");
		for (MemoryHeader *h = memoryLinkedList; h; h = h->next) 
		{
			_RPT3(_CRT_WARN, "    size=%d, file=%s, line=%d\n", h->size, h->file, h->line);
			sz += h->size;
		}
		_RPT0(_CRT_WARN, "    -----------------------------------------------------------------------\n");
		_RPT1(_CRT_WARN, "    %d bytes total\n\n", sz);
	}

class AutoReport
{
public:
	~AutoReport()
	{
		ReportMemoryLeaks();
	}
};
AutoReport ar;
	//----------------------------------------------------------------------------
	// CheckMemoryList... Check the list of allocated objects and ensure it's
	// valid.
	//
	void CheckMemoryList()
	{
		//
		// Check integrity of c-runtime library's heap.
		//
		//long flag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
		//flag |= _CRTDBG_CHECK_ALWAYS_DF;
		//_CrtSetDbgFlag(flag);
		//void CheckMemoryList();
		//char *test = (char *)malloc(1);
		//free(test);
		//flag &= ~_CRTDBG_CHECK_ALWAYS_DF;		// turn off
		//_CrtSetDbgFlag(flag);

		// Now check our list.
		for (MemoryHeader *m = memoryLinkedList; m; m = m->next) 
		{
			// Ensure magic numbers are correct.
			MemoryTrailer *trailer = (MemoryTrailer *)((unsigned char *)m + sizeof(MemoryHeader) + m->size);
			if ((m->magic != MEMORY_MAGIC_NUMBER) || (trailer->magic != MEMORY_MAGIC_NUMBER)) 
			{
				char buf[256];
				::sprintf(buf, "Memory overwrite detected.  File='%s', line=%d", m->file, m->line );
				Assert(0, buf);
			}
		}
	}
#else


	//----------------------------------------------------------------------------
	// new... Global operator new.  (RELEASE version)
	//
	void *operator new(unsigned int sz)
	{
		//
		// Allocate memory.
		//
		MemoryHeader *header = (MemoryHeader *)malloc(sz + sizeof(MemoryHeader));
		if (header == NULL) 
		{
			Fatal(false, "Out of memory");
			exit(1);
		}
		header->size = sz;

		//
		// Update global size.
		//
		memoryInUse += sz;

		//
		// Return pointer to data.
		//
		return ((void *)((unsigned char *)header + sizeof(MemoryHeader)));
	}


	//----------------------------------------------------------------------------
	// delete... Global operator delete.  (RELEASE version)
	//
	void operator delete(void *m)
	{
		if (m) 
		{
			MemoryHeader *header = (MemoryHeader *)((unsigned char *)m - sizeof(MemoryHeader));
			memoryInUse -= header->size;
			free(header);
		}
	}

#endif



//----------------------------------------------------------------------------
// GetMemoryUsage... Get memory usage.
//
long GetMemoryUsage()
{
	return memoryInUse;
}
