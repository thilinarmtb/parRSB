CC ?= mpicc
CFLAGS ?= -std=c99 -Wall -Wextra -Wpedantic -Wno-unused-function -Wno-c23-extensions
LDFLAGS ?=
DEBUG ?= 0
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASFLAGS ?= -lblas -llapack
GSLIBPATH ?=

########################## Don't touch what follows ###########################
ifeq ($(GSLIBPATH),)
	$(error Specify GSLIBPATH=<path to gslib build>)
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

MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))

PROJDIR = $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRCDIR = $(PROJDIR)/src
EXAMPLEDIR = $(PROJDIR)/example
BUILDDIR = $(PROJDIR)/build
INSTALLDIR = $(PROJDIR)/install
ifneq ($(strip $(DESTDIR)),)
	INSTALLDIR = $(DESTDIR)
endif

SRC.c = $(wildcard $(SRCDIR)/*.c)
SRC.o = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%.o,$(SRC.c))
EXAMPLE.c = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLE.bin = $(patsubst $(PROJDIR)/%.c,$(BUILDDIR)/%,$(EXAMPLE.c))
LIB = $(BUILDDIR)/lib/libparRSB.a

INCFLAGS = -I$(SRCDIR) -I$(GSLIBPATH)/include
CCCMD = $(CC) $(CFLAGS) $(INCFLAGS) $(PP)
LDFLAGS += -L$(INSTALLDIR)/lib -lparRSB -L$(GSLIBPATH)/lib -lgs -lm

.PHONY: all lib install example format clean

all: lib install example

lib: $(SRC.o) | $(BUILDDIR)
	@$(AR) cr $(LIB) $?
	@ranlib $(LIB)

install: lib | $(INSTALLDIR)
	@cp -v $(LIB) $(INSTALLDIR)/lib 2>/dev/null
	@cp $(SRCDIR)/*.h $(INSTALLDIR)/include 2>/dev/null

example: $(EXAMPLE.bin) | lib install

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
	$(CCCMD) -c $< -o $@

$(BUILDDIR)/%: $(PROJDIR)/%.c | lib install
	$(CCCMD) $< -o $@ $(LDFLAGS)

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)/lib
	@mkdir -p $(BUILDDIR)/src
	@mkdir -p $(BUILDDIR)/example

$(INSTALLDIR):
	@mkdir -p $(INSTALLDIR)/lib 2>/dev/null
	@mkdir -p $(INSTALLDIR)/include 2>/dev/null
