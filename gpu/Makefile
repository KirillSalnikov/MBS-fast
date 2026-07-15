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
GPU_PRECISION ?= float
GPU_FAST_MATH ?= 1
TARGET ?= bin/mbs_po_gpu_float_fast
FAST_SUFFIX := $(if $(filter 1,$(GPU_FAST_MATH)),_fast,)
MPI_SUFFIX := $(if $(filter 1,$(USE_MPI)),_mpi,)
GPU_ARCH ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.' | sed 's/[^0-9].*//')
GPU_ARCH := $(if $(GPU_ARCH),$(GPU_ARCH),86)
OBJDIR ?= build/$(GPU_PRECISION)$(FAST_SUFFIX)$(MPI_SUFFIX)_sm$(GPU_ARCH)/obj

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
ifeq ($(USE_MPI),1)
CXXFLAGS += -DUSE_MPI
endif
LDFLAGS ?= -lm -lgomp -L$(CUDA_LIB_PATH) -lcufft -lcudart
DEPFLAGS := -MMD -MP

NVCCFLAGS ?= -O3 -std=c++11 -arch=sm_$(GPU_ARCH) -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700
override NVCCFLAGS += -ccbin $(CUDA_HOST_CXX) -include $(CUDA_HOST_COMPAT)
override NVCCFLAGS += -Xcompiler -fopenmp
ifeq ($(GPU_PRECISION),float)
override NVCCFLAGS += -DMBS_GPU_FLOAT
CXXFLAGS += -DMBS_GPU_FLOAT
endif
ifeq ($(GPU_FAST_MATH),1)
override NVCCFLAGS += --use_fast_math
CXXFLAGS += -DMBS_GPU_FAST_MATH
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
	@echo "GPU_PRECISION: $(GPU_PRECISION)"
	@echo "GPU_FAST_MATH: $(GPU_FAST_MATH)"
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

float_fast:
	$(MAKE) GPU_PRECISION=float GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_float_fast all

double_fast:
	$(MAKE) GPU_PRECISION=double GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_double_fast all

double_fast_mpi:
	$(MAKE) USE_MPI=1 GPU_PRECISION=double GPU_FAST_MATH=1 TARGET=bin/mbs_po_gpu_mpi_double_fast all

float_debug:
	$(MAKE) GPU_PRECISION=float GPU_FAST_MATH=1 \
		TARGET=bin/mbs_po_gpu_float_debug \
		OBJDIR=build/float_debug/obj \
		CXXFLAGS="$(CXXFLAGS) -DMBS_DEBUG_HELP" all

double_debug:
	$(MAKE) GPU_PRECISION=double GPU_FAST_MATH=1 \
		TARGET=bin/mbs_po_gpu_double_debug \
		OBJDIR=build/double_debug/obj \
		CXXFLAGS="$(CXXFLAGS) -DMBS_DEBUG_HELP" all

clean:
	rm -rf build bin

.PHONY: all cuda_check float float_fast double_fast double_fast_mpi float_debug double_debug clean

-include $(DEPS)
