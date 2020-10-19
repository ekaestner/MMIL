
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Log.h"

SINGLETON_IMPLEMENT( Logger )

Logger::Logger() : outputStream( cout )//"/home/users/kevn/debug.txt", ios::out, filebuf::openprot )
{
	_indentLevel = 0;     // Initialy don't indent
	_debugLevel = 0;      // Should actually try to read env variable

	char* debug_lev = getenv("DEBUG_LEVEL");
	if(debug_lev != NULL)
	{
		 _debugLevel = atoi(debug_lev);
		 //outputStream << "DEBUG_LEVEL: Set to " << _debugLevel << endl << flush;
	} else {
		 //outputStream << "DEBUG_LEVEL: Not found. " << endl << flush;
		 //outputStream << "DEBUG_LEVEL: Defaults to " << _debugLevel << endl << flush;
	}
}

Logger::~Logger()
{
}

ostream& Logger::getStream( int level, int indentChange )
{
  if(indentChange < 0)
     _indentLevel += indentChange;

  // Insert the correct number of tabs into the stream for indenting
  for(int i=0;i<_indentLevel;i++)
     outputStream << "\t";

  if(indentChange > 0)
     _indentLevel += indentChange;

  return outputStream;
}

int Logger::getLevel()
{
  return _debugLevel;
}
