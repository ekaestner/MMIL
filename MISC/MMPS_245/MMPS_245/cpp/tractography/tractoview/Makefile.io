# objects to be included in the project....
IO_SRC = FileIO.c
#GLUTKERNEL = glutMain.cxx GlutKernelApplication.cxx
#DICOMLIB = DicomReader.cxx
SRCS =$(IO_SRC) $(INTERACTIONSRCS) $(GLUTKERNEL) $(DICOMLIB)  

# project name...
PROJECT_NAME = FileIO

# dicom toolkit directory
OFFIS_ROOT = /home/suabe/tool/dcmtk-3.5.4

OBJDIR = .
DEPDIR = 

PATH_TO_SRCS = .:fileIO/old
#PATH_TO_SRCS = .:Math:Interaction:OpenGL:glutkernel:Utils
INCLUDE_PATHS =  -I. -IMath -IInteraction -IOpenGL -Iglutkernel -Iexamples -IUtils -I/home/suabe/tool/opengl/includenew  \
		 -I$(OFFIS_ROOT)/dcmdata/include/dcmtk/dcmdata -I$(OFFIS_ROOT)/config/include/dcmtk/config \
		 -I$(OFFIS_ROOT)/ofstd/include/dcmtk/ofstd \
        	 -I$(OFFIS_ROOT)/dcmdata/include -I$(OFFIS_ROOT)/config/include -I$(OFFIS_ROOT)/ofstd/include

MK_DIR = mk

LDFLAGS = -L$(OFFIS_ROOT)/dcmdata/libsrc -L$(OFFIS_ROOT)/ofstd/libsrc -L/home/suabe/zlib-1.2.3

########################################################################
# system dependent stuff, that users will want to change...
# stuff that cant be abstracted in the build system...
##########################################################################
include ${MK_DIR}/mk.platform
include ${MK_DIR}/mk.hosttype

# IRIX
ifeq ($(HOSTTYPE), IRIX)
   LIBS = -lglut -lGLU -lGL -lX11 -lXmu
   include ${MK_DIR}/mk.irix
endif

#LiNUX

ifeq ($(HOSTTYPE), linux)
   LIBS = ./glutobj/glut_menu2.o ./glutobj/glut_hel12.o \
	  ./glutobj/glut_9x15.o ./glutobj/glut_bwidth.o \
	  -L/home/suabe/tool/opengl/libnew -lglut -lGLU -lGL  \
	  -L/home/suabe/tool/opengl/libnew -lmui \
	  -L/usr/X11R6/lib64 -lX11 -lXmu -lXt  -lSM -lICE -lXext -lXi $(LDFLAGS) -ldcmdata -lofstd  -lz -lm -lpthread \
	  -L/home/suabe/tool/opengl/libnew -lglut -lGLU -lGL  
   include ${MK_DIR}/mk.gnu
endif

# Win32
ifeq ($(HOSTTYPE), Win32)
#   LIBS = 
   DEPENDFLAGS = -D__cplusplus -D_WIN32 -DWIN32 $(INCLUDE_PATHS)
   include ${MK_DIR}/mk.win32
endif

########################################################################

include ${MK_DIR}/mk.objs

# Compile then Link
build: $(OBJS)
	-@echo ""
	-@echo "-------------------------------------------"
	-@echo "Linking... (.$(OBJ_FILE_SUFFIX)'s --> .so)"
	-@echo "-------------------------------------------"
	$(EXE_LINKER) $(EXE_FLAGS) $(OBJS) ./DicomReader.o $(OUTPUT_EXE_FLAG) $(EXE_FILENAME) $(LIBS)
	-@echo ""
	-@echo "-------------------------------------------"
	-@echo "$(EXE_FILENAME) done"
	-@echo "-------------------------------------------"

tarball: clean clobber
	make-tarball.sh Fractal src ${TARBALLPATH}

# Remove the compiled stuff from the system
clean:
	-@echo ""
	-@echo "-------------------------------------------"
	-@echo "Removing compiled stuff from the system	 "
	-@echo "-------------------------------------------"
	-rm -r $(OBJS) $(EXE_FILENAME) *.ncb *.opt *.plg *.ilk *.idb *.pdb *.pch Debug/ Release/ ii_files/



clobber: clean
	-rm *.d

include ${MK_DIR}/mk.depend
