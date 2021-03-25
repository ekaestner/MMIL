#ifndef SINGLETON_INCLUDED
#define SINGLETON_INCLUDED

// NOTES:
// + _instance needs to be a pointer
//   We need to be able to force construction with the classname::instance() method.
//
// + _instance as auto_ptr wont work in every case: "static auto_ptr<classname> instance;"
//   It can self construct after classname::instance() was called sometimes
//
// + Question:
//   Can we really guarentee that the _instance pointer will always be NULL before instance can be called?
//
// + Who deletes _instance when application is done?
//   auto_ptr looked good for awhile (see above)
//   a "Deleter" class is much better since it does nothing to affect the construction of _instance.

#define SINGLETON_DECLARE( classname )										\
	public:																	\
		static classname& instance();										\
		class classname##Deleter											\
		{																	\
		public:																\
			~classname##Deleter()											\
			{																\
				/*Log(3)<<"Killing "<<#classname<<"\n"<<flush;*/				\
				delete classname::_instance;								\
			}																\
		};																	\
	protected:																\
		classname();														\
		virtual ~classname();												\
	private:																\
		friend class classname##Deleter;									\
		static classname* _instance;										
																			
#define SINGLETON_IMPLEMENT( classname )									\
	classname* classname::_instance = NULL;									\
																			\
	classname& classname::instance()										\
	{																		\
		if (_instance == NULL)												\
		{																	\
			_instance = new classname;										\
		}																	\
		return *_instance;													\
	}																		\
																			\
	static classname::classname##Deleter classname##DeleterInstance;


#endif
