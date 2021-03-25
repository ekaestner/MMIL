#ifndef MEMORY_MANAGER
#define MEMORY_MANAGER


	#if defined(_DEBUG) || (DEBUG)
		// operator new
		// specifies byte size, filename, and linenumber.
		void*			operator new(unsigned int sz, const char *file, int line);

		// standard operator new (overloaded for our wonderful memory manager)
		void*			operator new( unsigned int sz );

		#define new		new (__FILE__, __LINE__)

		// standard operator delete (overloaded for our wonderful memory manager)
		void			operator delete(void *m);

		// operator delete, needed to compliment the operator new.
		//inline void		operator delete(void* memoryPtr, const char *file, unsigned int line)
		//{
		//	operator delete(memoryPtr);
		//}


		void			ReportMemoryLeaks();
	#else
		void*			operator new(unsigned int sz);

		// standard operator delete (overloaded for our wonderful memory manager)
		void			operator delete( void* m );

		inline void*	operator new(unsigned int s, int, const char *, int)
        { 
			return ::operator new(s);
		}

		#define New		new
	#endif



#endif
