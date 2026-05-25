ROOT := ..
SRC_DIR := $(ROOT)/src

ifeq ($(origin CXX),default)
CXX := $(shell if [ -x /usr/bin/mpicxx ]; then echo /usr/bin/mpicxx; elif command -v mpicxx >/dev/null 2>&1; then command -v mpicxx; else echo mpicxx; fi)
endif
TARGET ?= bin/mbs_po_mpi
OBJDIR ?= build/obj

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ARCH_FLAGS ?= $(shell bash $(ROOT)/scripts/detect_arch_flags.sh "$(CXX)")
OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
CXXFLAGS += -DUSE_MPI
LDFLAGS ?= -lm -lgomp
DEPFLAGS := -MMD -MP

INCLUDES := -I$(SRC_DIR) -I$(SRC_DIR)/math -I$(SRC_DIR)/handler \
            -I$(SRC_DIR)/common -I$(SRC_DIR)/geometry \
            -I$(SRC_DIR)/geometry/intrinsic -I$(SRC_DIR)/geometry/sse \
            -I$(SRC_DIR)/particle -I$(SRC_DIR)/scattering \
            -I$(SRC_DIR)/tracer -I$(SRC_DIR)/splitting \
            -I$(SRC_DIR)/cuda -I$(SRC_DIR)/bigint

SOURCES := $(shell find $(SRC_DIR) -path '$(SRC_DIR)/cuda' -prune -o -name '*.cpp' -print) \
           $(shell find $(SRC_DIR)/bigint -name '*.cc' 2>/dev/null)
CPU_SOURCES := GpuSupportStub.cpp
OBJECTS := $(patsubst $(ROOT)/%,$(OBJDIR)/%.o,$(SOURCES)) \
           $(patsubst %.cpp,$(OBJDIR)/cpu/%.o,$(CPU_SOURCES))
OBJECTS := $(OBJECTS:.cpp.o=.o)
OBJECTS := $(OBJECTS:.cc.o=.o)
DEPS := $(OBJECTS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	@echo "CPU MPI/OpenMP build"
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXX: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(ROOT)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR)/%.o: $(ROOT)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR)/cpu/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build bin

.PHONY: all clean

-include $(DEPS)
