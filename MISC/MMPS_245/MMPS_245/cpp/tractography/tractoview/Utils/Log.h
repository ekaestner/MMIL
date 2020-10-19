// Message logger.  
// 
#ifndef LOG_INCLUDED
#define LOG_INCLUDED



#define Log(level) \
	if(level <= Logger::instance().getLevel()) \
		Logger::instance().getStream(level)
#define LOG_BEGIN(val) \
	if(val <= Logger::instance().getLevel()) \
	    Logger::instance().getStream(val, 1)
#define LOG_END(val) \
	if(val <= Logger::instance().getLevel()) \
	    Logger::instance().getStream(val, -1)

//------------------------------------------
//: Class to support debug output
//
// Suggested use of val/debugLevel
//
// 1 - critical messages / Config data
// 2 -
// 3 - Object construction
// 4 -
// 5 - Highly verbose debug output
// 6 - Function entry and exit
// 7 - In house only type debug output
//-----------------------------------------
#include "SINGLETON.h"
class Logger
{
SINGLETON_DECLARE( Logger )

private:
    ostream& outputStream;

public:
   ostream& getStream(int level, int indentChange = 0);
   int getLevel();

private:
   int _debugLevel;      // Debug level to use
   int _indentLevel;     // Amount to indent
};



#endif
