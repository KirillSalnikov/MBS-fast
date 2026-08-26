ROOT := ..
SRC_DIR := $(ROOT)/src
CUDA_PATH ?= /usr/local/cuda
CUDA_LIB_PATH ?= $(if $(wildcard $(CUDA_PATH)/lib64/libcudart.so),$(CUDA_PATH)/lib64,$(CUDA_PATH)/targets/x86_64-linux/lib)

CXX_ORIGIN := $(origin CXX)
CXX ?= g++
NVCC ?= $(if $(wildcard $(CUDA_PATH)/bin/nvcc),$(CUDA_PATH)/bin/nvcc,nvcc)
CUDA_HOST_CXX ?= $(shell bash $(ROOT)/scripts/select_cuda_host_compiler.sh "$(NVCC)" "$(CXX)")
CUDA_HOST_COMPAT := $(abspath $(SRC_DIR)/cuda/CudaHostCompat.h)
ifeq ($(CXX_ORIGIN),default)
CXX := $(CUDA_HOST_CXX)
endif
USE_MPI ?= 0
GPU_PRECISION ?= double
GPU_FAST_MATH ?= 0
GPU_LTO ?= 1

ifeq ($(filter float fp32 double fp64,$(GPU_PRECISION)),)
$(error GPU_PRECISION must be one of: fp32, float, fp64, double)
endif
ifneq ($(filter 0 1,$(GPU_FAST_MATH)),$(GPU_FAST_MATH))
$(error GPU_FAST_MATH must be 0 or 1)
endif
ifneq ($(filter 0 1,$(GPU_LTO)),$(GPU_LTO))
$(error GPU_LTO must be 0 or 1)
endif
ifneq ($(filter 0 1,$(USE_MPI)),$(USE_MPI))
$(error USE_MPI must be 0 or 1)
endif

ifneq ($(filter float fp32,$(GPU_PRECISION)),)
GPU_PRECISION_CANON := float
GPU_PRECISION_LABEL := fp32
GPU_PRECISION_DEFINES := -DMBS_GPU_FP32 -DMBS_GPU_FLOAT
else
GPU_PRECISION_CANON := double
GPU_PRECISION_LABEL := fp64
GPU_PRECISION_DEFINES := -DMBS_GPU_FP64
endif

FAST_SUFFIX := $(if $(filter 1,$(GPU_FAST_MATH)),_fast,)
MPI_SUFFIX := $(if $(filter 1,$(USE_MPI)),_mpi,)
LTO_SUFFIX := $(if $(filter 1,$(GPU_LTO)),_lto,)
TARGET ?= bin/mbs_po_gpu$(MPI_SUFFIX)_$(GPU_PRECISION_CANON)$(FAST_SUFFIX)
GPU_ARCH ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.' | sed 's/[^0-9].*//')
GPU_ARCH := $(if $(GPU_ARCH),$(GPU_ARCH),86)
OBJDIR ?= build/$(GPU_PRECISION_LABEL)$(FAST_SUFFIX)$(MPI_SUFFIX)$(LTO_SUFFIX)_sm$(GPU_ARCH)/obj

ifeq ($(USE_MPI),1)
CXX := mpicxx
endif

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ARCH_FLAGS ?= $(shell bash $(ROOT)/scripts/detect_arch_flags.sh "$(CXX)")
GIT_DESCRIBE ?= $(shell git -C $(ROOT) describe --tags --always --dirty --long 2>/dev/null || echo unknown)
OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
CXXFLAGS += -DUSE_CUDA -DMBS_GPU_DEFAULT_ON -DMBS_GPU_BUILD_ARCH=$(GPU_ARCH) \
            -DMBS_GIT_DESCRIBE=\"$(GIT_DESCRIBE)\" -I$(CUDA_PATH)/include
EXTRA_CXXFLAGS ?=
CXXFLAGS += $(EXTRA_CXXFLAGS)
ifeq ($(USE_MPI),1)
CXXFLAGS += -DUSE_MPI
endif
LDFLAGS ?= -lm -lgomp -L$(CUDA_LIB_PATH) -lcufft -lcudart
GPU_LTO_FLAGS ?= -flto=auto
ifeq ($(GPU_LTO),1)
override CXXFLAGS += $(GPU_LTO_FLAGS)
override LDFLAGS += $(GPU_LTO_FLAGS)
endif
DEPFLAGS := -MMD -MP

NVCCFLAGS ?= -O3 -std=c++11 -arch=sm_$(GPU_ARCH) -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700
override NVCCFLAGS += -ccbin $(CUDA_HOST_CXX) -include $(CUDA_HOST_COMPAT)
override NVCCFLAGS += -Xcompiler -fopenmp
override NVCCFLAGS += $(GPU_PRECISION_DEFINES)
override CXXFLAGS += $(GPU_PRECISION_DEFINES)
ifeq ($(GPU_FAST_MATH),1)
override NVCCFLAGS += --use_fast_math
override CXXFLAGS += -DMBS_GPU_FAST_MATH
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

cuda_check:
	@$(NVCC) --version >/dev/null 2>&1 || { \
		echo "CUDA compiler not found: $(NVCC)" >&2; \
		echo "Set CUDA_PATH=/path/to/cuda or NVCC=/path/to/nvcc." >&2; \
		exit 1; \
	}
	@$(CUDA_HOST_CXX) --version >/dev/null 2>&1 || { \
		echo "CUDA-compatible host compiler not found: $(CUDA_HOST_CXX)" >&2; \
		echo "Install a GCC version supported by this CUDA toolkit or set CUDA_HOST_CXX." >&2; \
		exit 1; \
	}

$(OBJECTS): | cuda_check
$(TARGET): | cuda_check

$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	@echo "CUDA GPU build"
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXX: $(CXX)"
	@echo "CUDA_HOST_CXX: $(CUDA_HOST_CXX)"
	@echo "GPU_PRECISION: $(GPU_PRECISION_LABEL) ($(GPU_PRECISION_CANON) storage)"
	@echo "GPU_CRITICAL_PHASE: fp64"
	@echo "GPU_FAST_MATH: $(GPU_FAST_MATH)"
	@echo "GPU_LTO: $(GPU_LTO)"
	@echo "USE_MPI: $(USE_MPI)"
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

double:
	$(MAKE) GPU_PRECISION=double GPU_FAST_MATH=0 TARGET=bin/mbs_po_gpu_double all

float_fast:
	$(MAKE) GPU_PRECISION=float GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_float_fast all

double_fast:
	$(MAKE) GPU_PRECISION=double GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_double_fast all

fp32: float

fp64: double

fp32_fast: float_fast

fp64_fast: double_fast

precision_matrix: float double float_fast double_fast

double_fast_mpi:
	$(MAKE) USE_MPI=1 GPU_PRECISION=double GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_mpi_double_fast all

float_debug:
	$(MAKE) GPU_PRECISION=float GPU_FAST_MATH=1 \
		TARGET=bin/mbs_po_gpu_float_debug \
		OBJDIR=build/float_debug/obj \
		EXTRA_CXXFLAGS=-DMBS_DEBUG_HELP all

double_debug:
	$(MAKE) GPU_PRECISION=double GPU_FAST_MATH=1 \
		TARGET=bin/mbs_po_gpu_double_debug \
		OBJDIR=build/double_debug/obj \
		EXTRA_CXXFLAGS=-DMBS_DEBUG_HELP all

clean:
	rm -rf build bin

.PHONY: all cuda_check float double float_fast double_fast fp32 fp64 fp32_fast fp64_fast \
	precision_matrix double_fast_mpi float_debug double_debug clean

-include $(DEPS)
