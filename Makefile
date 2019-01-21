TARGET := ausgleichsproblem
SOURCES := $(TARGET).c
OBJS := $(SOURCES:.c=.o)
CFLAGS := -g -Wall -O1 
LIBS := -lm -lgsl
INCLUDES := -I/usr/local/include
CC := gcc
#$(shell which gcc)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean :
	rm $(TARGET) $(TARGET).o