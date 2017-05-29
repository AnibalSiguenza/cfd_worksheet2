# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c LBDefinitions.c 

# Compiler
# --------
CC=gcc

CFLAGS=-Werror -pedantic -Wall

#if mode variable is empty, setting release build mode                                                                                                
ifeq ($(mode),debug)                                                                                                                                  
CFLAGS += -O0 -g
endif    


# Linker flags
# ------------
LDFLAGS=-lm

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
