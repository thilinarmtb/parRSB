CC ?= mpicc
CFLAGS ?=
LDFLAGS ?=
DEBUG ?= 0
MPI ?= 1
UNDERSCORE ?= 1
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASDIR ?=
BLASFLAGS ?= -lblas -llapack

########################## Don't touch what follows ###########################
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib build>)
endif

MKFILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT := $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRC := $(SRCROOT)/src
EXAMPLE := $(SRCROOT)/examples
BUILD := $(SRCROOT)/build
INSTALL := $(BUILD)/install
ifneq ($(strip $(DESTDIR)),)
  INSTALL = $(realpath $(DESTDIR))
endif

SRCS = $(patsubst $(SRCROOT)/%.c,$(BUILD)/%.o,$(wildcard $(SRC)/*.c))
EXAMPLES = $(patsubst $(SRCROOT)/%.c,$(BUILD)/%,$(wildcard $(EXAMPLE)/*.c))

LIB = $(BUILD)/lib/libparRSB.a

ifneq ($(DEBUG),0)
  PP += -DPARRSB_DEBUG
  CFLAGS += -g
else
  CFLAGS += -O2
endif

ifneq ($(MPI),0)
  PP += -DMPI
endif

ifneq ($(UNDERSCORE),0)
  PP += -DPARRSB_UNDERSCORE
endif

ifneq ($(SYNC_BY_REDUCTION),0)
  PP += -DPARRSB_SYNC_BY_REDUCTION
endif

ifneq ($(BLAS),0)
  PP += -DPARRSB_BLAS
  ifneq ($(BLASDIR),)
    LDFLAGS+= -L$(BLASDIR)
  endif
  LDFLAGS += $(BLASFLAGS)
endif

INCFLAGS = -I$(SRC) -I$(GSLIBPATH)/include
CCCMD = $(CC) $(CFLAGS) $(INCFLAGS) $(PP)
LDFLAGS += -lm

.PHONY: all lib install examples format clean

all: lib install

lib: $(SRCS)
	@mkdir -p $(BUILD)/lib
	@$(AR) cr $(LIB) $?
	@ranlib $(LIB)

install: lib
	@mkdir -p $(INSTALL)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALL)/lib 2>/dev/null
	@mkdir -p $(INSTALL)/include 2>/dev/null
	@cp $(SRC)/*.h $(INSTALL)/include 2>/dev/null

examples: install $(EXAMPLES)

format:
	find . -iname *.h -o -iname *.c -o -iname *.okl | xargs clang-format -i

clean:
	@$(RM) -rf $(BUILD)

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(BUILD)/src/%.o: $(SRC)/%.c | dir
	$(CCCMD) -c $< -o $@

$(BUILD)/examples/%: $(EXAMPLE)/%.c | lib install dir
	$(CCCMD) $< -o $@ -L$(INSTALL)/lib -lparRSB -L$(GSLIBPATH)/lib -lgs $(LDFLAGS)

dir:
	mkdir -p $(BUILD)/src $(BUILD)/examples
