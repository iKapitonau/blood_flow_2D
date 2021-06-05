CC ?= gcc
LDFLAGS := $(LDFLAGS) -lconfig -lm -fopenmp
CFLAGS := $(CFLAGS) -Wall -Wextra -pedantic -Iinclude/

exec = blood_flow_2d
srcs = $(wildcard src/*.c)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.dep)

.PHONY: clean debug release parallel

release: CFLAGS += -O3
release: $(exec)

debug: CFLAGS += -g
debug: $(exec)

parallel: CFLAGS += -fopenmp -DWITH_OPENMP_
parallel: release

$(exec): $(objs)
	$(CC) $(objs) $(LDFLAGS) -o $@

ifneq ($(MAKECMDGOALS), clean)
-include $(deps)
endif

%.dep: %.c
	$(CC) $(CFLAGS) $< -MM > $@

clean:
	-rm -f $(objs) $(deps) $(exec)
