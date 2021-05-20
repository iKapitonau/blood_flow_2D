CC ?= gcc
LDFLAGS := $(LDFLAGS) -lconfig -lm
CFLAGS := $(CFLAGS) -Wall -Wextra -pedantic -Iinclude/ -g

srcs = $(wildcard src/*.c)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.dep)

.PHONY: all clean

all: $(objs)
	$(CC) $(objs) $(LDFLAGS) -o all

ifneq ($(MAKECMDGOALS), clean)
-include $(deps)
endif

%.dep: %.c
	$(CC) $(CFLAGS) $< -MM > $@

clean:
	-rm -f $(objs) $(deps)
