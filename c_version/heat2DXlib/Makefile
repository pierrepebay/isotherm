# Makefile pour les applications X

CC = gcc 
CFLAGS = -I/usr/X11/include
LibFLAGS = -L/usr/X11/lib
LDLIBS = -lX11 -lm

SRC = dataMain.c processingModule.c
OBJ = $(SRC:%.c=%.o)

dataMain : $(OBJ) processingModule.h
	$(CC)  $(LibFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o : %.c dataProcessing.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean

clean:	
	\rm  dataMain *.o saveData

