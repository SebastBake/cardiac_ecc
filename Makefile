MAKEFLAGS= --no-builtin-rules --warn-undefined-variables

GLOBAL_ROOT :=${OPENCMISS_ROOT}/cm
OC_CM_GLOBAL_ROOT = $(GLOBAL_ROOT)

include $(OC_CM_GLOBAL_ROOT)/utils/MakefileCommon.inc

BASE_PROGRAM_NAME = cardiac_ecc
PROGRAM_ROOT = /home/vijay/Documents/heart/sims/opencmiss/cardiac_ecc
OBJECT_DIR := $(PROGRAM_ROOT)/object/$(ENVIRONMENT)$(BUILD_TYPE)$(MPI_TOOLCHAIN)
EXE_DIR := $(PROGRAM_ROOT)/bin/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)_$(COMPILER_VERSION)


MAIN_SOURCE_DIR := $(OC_CM_GLOBAL_ROOT)/src
BASE_EXE_NAME = $(BASE_PROGRAM_NAME)
BASE_LIB_NAME = OpenCMISS

PROGRAM_SOURCE_DIR = $(PROGRAM_ROOT)/src
MODULE_DIR := $(OBJECT_DIR)
INC_NAME := opencmiss.mod
INCLUDE := $(INC_DIR)/$(INC_NAME)
LIB_NAME := lib$(BASE_LIB_NAME)$(EXE_ABI_SUFFIX)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a
LIBRARY := $(LIB_DIR)/$(LIB_NAME)
EXE_NAME := $(BASE_EXE_NAME)$(EXE_ABI_SUFFIX)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
EXECUTABLE := $(EXE_DIR)/$(EXE_NAME)

C_INCLUDE_DIRS := $(INC_DIR) \
    $(EXAMPLE_SOURCE_DIR)
F_INCLUDE_DIRS := $(MODULE_DIR)

CPPFLAGS += $(addprefix -I, $(C_INCLUDE_DIRS))
FPPFLAGS += $(addprefix -I, $(F_INCLUDE_DIRS))

OPENCMISS_LIBRARY = $(LIBRARY)
OPENCMISS_INCLUDE_PATH = $(addprefix -I, $(INC_DIR))

FPPFLAGS += $(OPENCMISS_INCLUDE_PATH) $(EXTERNAL_INCLUDE_PATH)

ELFLAGS += $(EXTERNAL_LIB_PATH)

.SUFFIXES: .f90 .c

OBJECTS = $(OBJECT_DIR)/cardiac_ecc.o

main: preliminaries $(EXECUTABLE)
preliminaries: $(OBJECT_DIR)$(EXE_DIR)

$(OBJECT_DIR) :
    mkdir -p $@

$(EXE_DIR) :
    mkdir -p $@

$(EXECUTABLE) : $(OBJECTS) $(OPENCMISS_LIBRARY)
    $(EXE_LINK) -o $@ $(OBJECTS) $(OPENCMISS_LIBRARY) $(ELFLAGS) $(EXTERNAL_LIBRARIES)

$(OBJECT_DIR)/cardiac_ecc.o : $(PROGRAM_SOURCE_DIR)/cardiac_ecc.f90
( cd $(OBJECT_DIR) ; $(FC) $(FFLAGS) $(FPPFLAGS) -c $<)


$(OPENCMISS_LIBRARY):
    ( cd $(OC_CM_GLOBAL_ROOT)l $(MAKE))

clean:
    @echo "Cleaning house"
    rm -rf $(OBJECT_DIR) $(EXECUTABLE)

allclean:
    @echo "Cleaning house"
    rm -rf object/ bin/

####
#Aliases
####
include $(OCE_MAKEINC_ROOT)/Makefile_Aliases.inc