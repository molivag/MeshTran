# Compiler
FC = ifort

# Directories
SRC_DIR = src
MOD_DIR = mod
BIN_DIR = bin

# Target
TARGET = meshTran

# Source files
SRC = \
$(SRC_DIR)/mesh_config.f90 \
$(SRC_DIR)/mesh_entities.f90 \
$(SRC_DIR)/meshTranFemtic.f90

# Flags
FFLAGS = -g -traceback -check all -fpe0 -warn all -debug all \
         -mcmodel large -fp-stack-check -check noarg_temp_created \
         -diag-disable=10448 \
         -module $(MOD_DIR)

# Default target
all: $(TARGET)

$(MOD_DIR):
	@mkdir -p $(MOD_DIR)

# Link executable
$(TARGET): $(MOD_DIR) $(SRC)
	@$(FC) $(FFLAGS) $(SRC) -o $@

# Clean build
clean:
	@rm -rf $(MOD_DIR) 
	@rm meshTran

.PHONY: all clean
