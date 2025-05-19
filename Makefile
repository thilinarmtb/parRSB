CC ?= mpicc
CFLAGS ?= -std=c99 -Wall -Wextra -Wpedantic -Wno-unused-function -Wno-c23-extensions
GSLIBPATH ?=
LDFLAGS ?=
INCFLAGS ?=
DEBUG ?= 0
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASFLAGS ?= -lblas
SHARED ?= 0

########################## Don't touch what follows ###########################
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib build>)
else
  LDFLAGS += -L$(GSLIBPATH)/lib -lgs
  INCFLAGS += -I$(GSLIBPATH)/include
endif

ifneq ($(DEBUG),0)
  CFLAGS += -g
else
  CFLAGS += -O2
endif

PP =
ifneq ($(SYNC_BY_REDUCTION),0)
  PP += -DPARRSB_SYNC_BY_REDUCTION
endif

ifneq ($(BLAS),0)
  PP += -DPARRSB_BLAS
  LDFLAGS += $(BLASFLAGS)
endif

LIB.ext = a
ifneq ($(SHARED),0)
  CFLAGS += -fPIC
  LIB.ext = so
endif

MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))
PROJDIR = $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRCDIR = $(PROJDIR)/src
BUILDDIR = $(PROJDIR)/build
INSTALLDIR = $(PROJDIR)/install
ifneq ($(strip $(DESTDIR)),)
  INSTALLDIR = $(DESTDIR)
endif
LIBDIR = $(INSTALLDIR)/lib
INCDIR = $(INSTALLDIR)/include

SRC.c = $(wildcard $(SRCDIR)/*.c)
SRC.o = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%.o,$(SRC.c))
EXAMPLE.c = $(wildcard $(PROJDIR)/example/*.c)
EXAMPLE.bin = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%,$(EXAMPLE.c))
TEST.c = $(wildcard $(PROJDIR)/test/*.c)
TEST.bin = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%,$(TEST.c))
LIB = $(BUILDDIR)/lib/libparRSB.$(LIB.ext)
LIB.a = $(BUILDDIR)/lib/libparRSB.a
LIB.so = $(BUILDDIR)/lib/libparRSB.so

LDFLAGS += -lm
INCFLAGS += -I$(SRCDIR)

.PHONY: all lib install example format clean

all: lib install example test

lib: $(LIB)

$(LIB.a): $(SRC.o) | $(BUILDDIR)
	@$(AR) cr $@ $?
	@ranlib $@

$(LIB.so): $(SRC.o)
	$(CC) -shared -o $@ $(SRC.o) $(LDFLAGS)

install: lib | $(INCDIR) $(LIBDIR)
	@cp -v $(LIB) $(INSTALLDIR)/lib
	@cp -v $(SRCDIR)/*.h $(INSTALLDIR)/include

example: $(EXAMPLE.bin) | lib install

test: $(TEST.bin) | lib install

format:
	find . -iname *.h -o -iname *.c | xargs clang-format -i

clean:
	@$(RM) -rf $(BUILDDIR)

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(BUILDDIR)/%.o: $(PROJDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) $(INCFLAGS) $(PP) -c $< -o $@

$(BUILDDIR)/%: $(PROJDIR)/%.c | lib install
	$(CC) $(INCFLAGS) -I$(SRCDIR)/test $< -o $@ -L$(LIBDIR) -lparRSB $(LDFLAGS)

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)/lib
	@mkdir -p $(BUILDDIR)/src
	@mkdir -p $(BUILDDIR)/example
	@mkdir -p $(BUILDDIR)/test

$(INCDIR):
	@mkdir -p $(INSTALLDIR)/include

$(LIBDIR):
	@mkdir -p $(INSTALLDIR)/lib
