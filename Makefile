# Compiler
FC = ifort

# Directories
SRC_DIR = src
MOD_DIR = mod
BIN_DIR = bin

# Target
TARGET = ./MeshTran

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

# Link executable
$(TARGET): 	$(SRC)
	@$(FC) $(FFLAGS) $^ -o $@

# Clean build
clean:
	@rm -rf $(MOD_DIR)/*.mod
	@rm -f *.o
	@rm MeshTran

.PHONY: all clean
