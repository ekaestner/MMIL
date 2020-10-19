default:
	@echo "specify platform: win32 or unix"

win32:
	cl.exe refcount_tpl.cxx /I../. /Ferefcount_tpl.exe
	cl.exe memman.cxx /I../. /Fememman.exe
	cl.exe refcount_inh.cxx /I../. /Ferefcount_inh.exe
	cl.exe endian.cxx /I../. /Feendian.exe
   
unix:
	CC -g -woff 3303 random.cxx -I../. -orandom -lm
	CC -g -woff 3303 refcount_inh.cxx -I../. -orefcount_inh
	CC -g -woff 3303 refcount_tpl.cxx -I../. -orefcount_tpl
	CC -g -woff 3303 endian.cxx -I../. -oendian
