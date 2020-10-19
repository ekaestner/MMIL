# need next one to get -cckr at -c time for analyse/fourier/deconv
#CFLAGS =        -O2 -mips2 -I$(LOCAL_COMPILE_DIR)/include -cckr
CFLAGS =	-g -mips2 -I$(LOCAL_COMPILE_DIR)/include

LDFLAGS =	-g -L$(LOCAL_COMPILE_DIR)/lib

LLIBS =		-ltk4.0 -ltcl7.4
LIBS =		-lX11 -lgl -lm

CC =		cc

OBJ1 = \
  rcxcopy.o \
  surflib.o

OBJ2 = \
  rcxcombine.o \
  surflib.o \
  ttestlib.o

OBJ3 = \
  bfloat2txt.o \
  surflib.o

OBJ4 = \
  bresize.o \
  surflib.o

OBJ5 = \
  smoothsurfvals.o \
  surflib.o

OBJ6 = \
  w2bfloat.o \
  surflib.o

OBJ7 = \
  readsurfnverts.o \
  surflib.o

OBJ8 = \
  tal2w.o \
  surflib.o

OBJ9 = \
  wmask.o \
  surflib.o

default:  rcxcopy rcxcombine bfloat2txt bresize \
          smoothsurfvals w2bfloat readsurfnverts \
          tal2w wmask

rcxcopy:	$(OBJ1)
		$(CC) $(CFLAGS) $(OBJ1) -o $@ $(LIBS)

rcxcombine:	$(OBJ2)
		$(CC) $(CFLAGS) $(OBJ2) -o $@ $(LIBS)

bfloat2txt:	$(OBJ3)
		$(CC) $(CFLAGS) $(OBJ3) -o $@ $(LIBS)

bresize:	$(OBJ4)
		$(CC) $(CFLAGS) $(OBJ4) -o $@ $(LIBS)

smoothsurfvals:	$(OBJ5)
		$(CC) $(CFLAGS) $(OBJ5) -o $@ $(LIBS)

w2bfloat:	$(OBJ6)
		$(CC) $(CFLAGS) $(OBJ6) -o $@ $(LIBS)

readsurfnverts:	$(OBJ7)
		$(CC) $(CFLAGS) $(OBJ7) -o $@ $(LIBS)

tal2w:	$(OBJ8)
		$(CC) $(CFLAGS) $(OBJ8) -o $@ $(LIBS)

wmask:	$(OBJ9)
		$(CC) $(CFLAGS) $(OBJ9) -o $@ $(LIBS)

clean:
		rm -f *.o core rcxcopy rcxcombine bfloat2txt \
                   smoothsurfvals w2bfloat readsurfnverts \
                   tal2w wmask

