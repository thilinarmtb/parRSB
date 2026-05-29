CC ?= mpicc
CFLAGS ?= -std=c99 -Wall -Wextra -Wpedantic -Wno-unused-function -Wno-c23-extensions
GSLIB_DIR ?=
LDFLAGS ?=
INCFLAGS ?=
DEBUG ?= 0
SYNC ?= reduction
BLAS ?= 0
BLASFLAGS ?= -lblas
SHARED ?= 0
VISIBILITY ?= hidden

########################## Don't touch what follows ###########################
ifeq ($(GSLIB_DIR),)
  $(error Specify GSLIB_DIR=<path to gslib build>)
else
  LDFLAGS += -L$(GSLIB_DIR)/lib -lgs
  INCFLAGS += -I$(GSLIB_DIR)/include
endif

ifneq ($(DEBUG),0)
  CFLAGS += -g
else
  CFLAGS += -O2
endif

PP =
ifeq ($(strip $(SYNC)),reduction)
  PP += -DPARRSB_SYNC_REDUCTION
else ifeq ($(strip $(SYNC)),barrier)
  PP += -DPARRSB_SYNC_BARRIER
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

ifneq ($(strip $(VISIBILITY)),)
  PP += -DPARRSB_VISIBILITY=$(strip $(VISIBILITY))
endif

MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))
PROJDIR = $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRCDIR = $(PROJDIR)/src
BUILDDIR = $(PROJDIR)/build
INSTALLDIR = $(PROJDIR)/install
ifneq ($(strip $(DESTDIR)),)
  INSTALLDIR = $(DESTDIR)
endif

SRC.c = $(wildcard $(SRCDIR)/*.c)
SRC.o = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%.o,$(SRC.c))
EXAMPLES.c = $(wildcard $(PROJDIR)/examples/*.c)
EXAMPLES.bin = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%,$(EXAMPLES.c))
TESTS.c = $(wildcard $(PROJDIR)/tests/*.c)
TESTS.bin = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%,$(TESTS.c))
LIB = $(BUILDDIR)/lib/libparRSB.$(LIB.ext)
LIB.a = $(BUILDDIR)/lib/libparRSB.a
LIB.so = $(BUILDDIR)/lib/libparRSB.so

LDFLAGS += -lm
INCFLAGS += -I$(SRCDIR)

.PHONY: all lib install examples format clean

all: lib install examples tests

lib: $(LIB)

$(LIB.a): $(SRC.o)
	@$(AR) cr $@ $?
	@ranlib $@

$(LIB.so): $(SRC.o)
	$(CC) -shared -o $@ $? $(LDFLAGS)

install: lib
	@mkdir -p $(INSTALLDIR)/include
	@mkdir -p $(INSTALLDIR)/lib
	@cp -v $(LIB) $(INSTALLDIR)/lib/
	@cp -v $(SRCDIR)/*.h $(INSTALLDIR)/include/

examples: $(EXAMPLES.bin) | lib install

tests: $(TESTS.bin) | lib install

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

$(BUILDDIR)/%.o: $(PROJDIR)/%.c | build
	$(CC) $(CFLAGS) $(INCFLAGS) $(PP) -c $< -o $@

$(BUILDDIR)/%: $(PROJDIR)/%.c | lib install
	$(CC) $(INCFLAGS) -I$(SRCDIR)/tests $< -o $@ -L$(INSTALLDIR)/lib -lparRSB $(LDFLAGS)

build:
	@mkdir -p $(BUILDDIR)/lib
	@mkdir -p $(BUILDDIR)/src
	@mkdir -p $(BUILDDIR)/examples
	@mkdir -p $(BUILDDIR)/tests
