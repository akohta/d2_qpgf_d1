CC      =gcc
CFLAGS  =-O3 -Wall
LDFLAGS =-lm
SRCDIR  =src
OBJDIR  =src
TARGET1 =example.out
TRGSRC1 =example.c

SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS)) ))
HEAD=$(wildcard $(SRCDIR)/*.h)

TRGOBJ1=$(OBJS) $(patsubst %.c,%.o,$(TRGSRC1))

all : $(TARGET1)

$(TARGET1) : $(TRGOBJ1) 
	$(CC) $(CFLAGS) $^  -o $@ $(LDFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I $(SRCDIR) -c $< -o $@

.c.o :
	$(CC) $(CFLAGS) -I$(SRCDIR) -c $<


clean:
	@rm -rf $(TARGET1) *.o $(OBJDIR)/*.o
        
$(OBJS) : $(HEAD)

