DEBUG ?= 0
PAUL ?= 1
CC ?= mpicc
CFLAGS ?= -O2
## ALGO=0 (Lanczos),1 (RQI),2 (FMG)
ALGO ?= 0
MPI ?= 1

MKFILEPATH =$(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT_  ?=$(patsubst %/,%,$(dir $(MKFILEPATH)))
SRCROOT    =$(realpath $(SRCROOT_))

GSLIBDIR=$(GSLIBPATH)

SRCDIR    =$(SRCROOT)/src
SORTDIR   =$(SRCROOT)/src/sort
PRECONDDIR=$(SRCROOT)/src/precond
GENCONDIR =$(SRCROOT)/src/gencon
BUILDDIR  =$(SRCROOT)/build
EXAMPLEDIR=$(SRCROOT)/examples
TESTDIR   =$(SRCROOT)/tests

TARGET=parRSB
LIB=$(BUILDDIR)/lib/lib$(TARGET).a

INCFLAGS=-I$(SRCDIR) -I$(SORTDIR) -I$(PRECONDDIR) -I$(GENCONDIR) -I$(GSLIBDIR)/include
LDFLAGS:=-L$(BUILDDIR)/lib -l$(TARGET) -L $(GSLIBDIR)/lib -lgs -lm

SRCS       =$(wildcard $(SRCDIR)/*.c)
SORTSRCS   =$(wildcard $(SORTDIR)/*.c)
PRECONDSRCS=$(wildcard $(PRECONDDIR)/*.c)
GENCONSRCS =$(wildcard $(GENCONDIR)/*.c)
TESTSRCS   =$(wildcard $(TESTDIR)/*.c)
EXAMPLESRCS=$(wildcard $(EXAMPLEDIR)/*.c)

SRCOBJS =$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SRCS))
SRCOBJS+=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SORTSRCS))
SRCOBJS+=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(PRECONDSRCS))
SRCOBJS+=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(GENCONSRCS))
TESTOBJS=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%,$(TESTSRCS))
EXAMPLEOBJS=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%,$(EXAMPLESRCS))

PP=

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifeq ($(ALGO),0)
  PP += -DGENMAP_LANCZOS
else ifeq ($(ALGO),1)
  PP += -DGENMAP_RQI
endif

ifneq ($(PAUL),0)
  PP += -DGENMAP_PAUL
endif

ifneq ($(MPI),0)
  PP += -DMPI
endif

INSTALLDIR=
ifneq (,$(strip $(DESTDIR)))
	INSTALLDIR=$(realpath $(DESTDIR))
endif

.PHONY: default
default: check lib install

.PHONY: all
all: check lib tests examples install

.PHONY: install
install: lib
ifneq ($(INSTALLDIR),)
	@mkdir -p $(INSTALLDIR)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLDIR)/lib 2>/dev/null
	@mkdir -p $(INSTALLDIR)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SORTDIR)/*.h $(PRECONDDIR)/*.h $(INSTALLDIR)/include 2>/dev/null
endif

.PHONY: lib
lib: $(SRCOBJS)
	@mkdir -p $(BUILDDIR)/lib
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

.PHONY: check
check: 
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib>/build)
endif

$(BUILDDIR)/src/%.o: $(SRCROOT)/src/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: examples
examples: install $(EXAMPLEOBJS)

$(BUILDDIR)/examples/%: $(SRCROOT)/examples/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(LDFLAGS)

.PHONY: tests
tests: install $(TESTOBJS)
	@cp $(TESTDIR)/run-tests.sh $(BUILDDIR)/tests/
	@cd $(BUILDDIR)/tests && ./run-tests.sh

$(BUILDDIR)/tests/%: $(SRCROOT)/tests/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -rf $(BUILDDIR) $(EXAMPLE) $(EXAMPLE).o

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/src/sort)
$(shell mkdir -p $(BUILDDIR)/src/precond)
$(shell mkdir -p $(BUILDDIR)/src/gencon)
$(shell mkdir -p $(BUILDDIR)/tests)
$(shell mkdir -p $(BUILDDIR)/examples)
