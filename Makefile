SHELL = /bin/sh

# libraries
Linux_LIB = 
Linux_CLIB = 
Linux_MYLIB = 
LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_MYLIB)

# flags
Linux_OPTFLAGS = -fast -xSSE3
Linux_CFLAGS = $(${OS}_OPTFLAGS)

# executable
MAINNAM = nanocrystal

# compilers
Linux_CC = icc
Linux_LD = icc

# sources
OBJECTS = \
	main.o build.o cut.o passivate.o chiral.o bonds.o size.o atom.o vector.o read.o write.o errorHandling.o  

# make and make clean
.c.o:
	$(${OS}_CC) $(${OS}_CFLAGS) -c  $*.c 	

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
