
include ../../Makeheader

INCLUDE = -I../include
CPPFLAGS += $(INCLUDE)

OBJDIR = ../../obj/fft
SRCS = $(wildcard *.c)
OBJECTS = $(SRCS:%.c=$(OBJDIR)/%.o)

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

all: $(OBJECTS)