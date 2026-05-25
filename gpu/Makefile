ROOT := ..
SRC_DIR := $(ROOT)/src
CUDA_PATH ?= /usr/local/cuda

CXX ?= g++
NVCC ?= nvcc
GPU_PRECISION ?= float
GPU_FAST_MATH ?= 1
TARGET ?= bin/mbs_po_gpu_float_fast
FAST_SUFFIX := $(if $(filter 1,$(GPU_FAST_MATH)),_fast,)
OBJDIR ?= build/$(GPU_PRECISION)$(FAST_SUFFIX)/obj
GPU_ARCH ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.' | sed 's/[^0-9].*//')
GPU_ARCH := $(if $(GPU_ARCH),$(GPU_ARCH),86)

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ARCH_FLAGS ?= $(shell bash $(ROOT)/scripts/detect_arch_flags.sh "$(CXX)")
OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
CXXFLAGS += -DUSE_CUDA -DMBS_GPU_DEFAULT_ON -DMBS_GPU_BUILD_ARCH=$(GPU_ARCH) -I$(CUDA_PATH)/include
LDFLAGS ?= -lm -lgomp -L$(CUDA_PATH)/lib64 -lcufft -lcudart
DEPFLAGS := -MMD -MP

NVCCFLAGS ?= -O3 -std=c++11 -arch=sm_$(GPU_ARCH) -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700
ifeq ($(GPU_PRECISION),float)
override NVCCFLAGS += -DMBS_GPU_FLOAT
endif
ifeq ($(GPU_FAST_MATH),1)
override NVCCFLAGS += --use_fast_math
endif

INCLUDES := -I$(SRC_DIR) -I$(SRC_DIR)/math -I$(SRC_DIR)/handler \
            -I$(SRC_DIR)/common -I$(SRC_DIR)/geometry \
            -I$(SRC_DIR)/geometry/intrinsic -I$(SRC_DIR)/geometry/sse \
            -I$(SRC_DIR)/particle -I$(SRC_DIR)/scattering \
            -I$(SRC_DIR)/tracer -I$(SRC_DIR)/splitting \
            -I$(SRC_DIR)/cuda -I$(SRC_DIR)/bigint

SOURCES := $(shell find $(SRC_DIR) -name '*.cpp') \
           $(shell find $(SRC_DIR)/bigint -name '*.cc' 2>/dev/null) \
           $(shell find $(SRC_DIR) -name '*.cu')
OBJECTS := $(patsubst $(ROOT)/%,$(OBJDIR)/%.o,$(SOURCES))
OBJECTS := $(OBJECTS:.cpp.o=.o)
OBJECTS := $(OBJECTS:.cc.o=.o)
OBJECTS := $(OBJECTS:.cu.o=.o)
DEPS := $(OBJECTS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	@echo "CUDA GPU build"
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXX: $(CXX)"
	@echo "GPU_PRECISION: $(GPU_PRECISION)"
	@echo "GPU_FAST_MATH: $(GPU_FAST_MATH)"
	@echo "GPU_ARCH: sm_$(GPU_ARCH)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "NVCCFLAGS: $(NVCCFLAGS)"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(ROOT)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR)/%.o: $(ROOT)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR)/%.o: $(ROOT)/%.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -DUSE_CUDA -I$(CUDA_PATH)/include $(INCLUDES) -c $< -o $@

float:
	$(MAKE) GPU_PRECISION=float GPU_FAST_MATH=0 TARGET=bin/mbs_po_gpu_float all

float_fast:
	$(MAKE) GPU_PRECISION=float GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_float_fast all

double_fast:
	$(MAKE) GPU_PRECISION=double GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_double_fast all

clean:
	rm -rf build bin

.PHONY: all float float_fast double_fast clean

-include $(DEPS)
