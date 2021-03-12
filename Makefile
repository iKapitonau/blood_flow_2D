CC ?= gcc
LDFLAGS := $(LDFLAGS) -lconfig -lm
CFLAGS := $(CFLAGS) -Wall -Wextra -pedantic -Iinclude/

srcs = $(wildcard src/*.c)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.dep)

.PHONY: all clean

all: $(objs)

ifneq ($(MAKECMDGOALS), clean)
-include $(deps)
endif

%.dep: %.c
	$(CC) $(CFLAGS) $< -MM > $@

clean:
	-rm -f $(objs) $(deps)
