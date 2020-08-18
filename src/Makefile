IDIR =../libStatGen/include

CC=g++ -std=c++0x -O4
CFLAGS=-I$(IDIR)   -L../libStatGen

ODIR=obj
LDIR =../libStatGen/

LIBS=-lStatGen -lz -lpthread

#_DEPS = hellomake.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fastqproc: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
