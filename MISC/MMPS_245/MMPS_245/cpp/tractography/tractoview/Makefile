# objects to be included in the project....
MATHSRCS = Matrix4f.cxx Light.cxx Material.cxx ColorRGBA.cxx Quatf.cxx
INTERACTIONSRCS = TrackBall.cxx HyperBall.cxx SphereBall.cxx
SRCS =$(MATHSRCS) $(INTERACTIONSRCS) Lsystem.cxx tractomain.cxx

ROOTDIR = /home/mmildev/tracto/tracto_dh

# project name...
PROJECT_NAME = tractoview

OBJDIR = .
DEPDIR = 

PATH_TO_SRCS = .:tractomain:Utils:Math:Interaction:OpenGL

TOOLS_DIR = /home/mmildev/tracts/tools

INCLUDE_PATHS =  \
  -I. -IMath -IInteraction -IOpenGL -Itractomain -IUtils \
  -I$(TOOLS_DIR)/opengl/includenew  \
  -I$(TOOLS_DIR)/qt-x11-free-3.3.3/mkspecs/linux-g++ \
  -I$(TOOLS_DIR)/qt-x11-free-3.3.3/include \
  -I/usr/include/freetype2 \
  -I/usr/X11R6/include \
  -Itiff-v3.6.1/include \

MK_DIR = mk

LDFLAGS =

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
  LIBS = \
    ./glutobj/glut_menu2.o ./glutobj/glut_hel12.o \
    ./glutobj/glut_9x15.o ./glutobj/glut_bwidth.o \
    $(TOOLS_DIR)/opengl/libnew/libglut.a  -lGLU -lGL \
    $(TOOLS_DIR)/opengl/libnew/libmui.a \
    -L$(TOOLS_DIR)/qt-x11-free-3.3.3/lib -lqt -lXrandr \
    -lXcursor -lXinerama \
    -lXft -lfreetype -lfontconfig \
    -L/usr/X11R6/lib64 -lX11 -lXmu -lXt  -lSM -lICE \
    -lXext -lXi $(LDFLAGS) -lz -lm -lpthread \
    -Ltiff-v3.6.1/lib  -ltiff
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
	$(EXE_LINKER) $(EXE_FLAGS) $(OBJS) $(OUTPUT_EXE_FLAG) $(EXE_FILENAME) $(LIBS)
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
