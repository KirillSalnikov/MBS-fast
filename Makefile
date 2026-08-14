CXX_ORIGIN := $(origin CXX)
CXX ?= g++
CUDA_PATH ?= /usr/local/cuda
CUDA_LIBDIR ?= $(or $(firstword $(wildcard $(CUDA_PATH)/lib64 $(CUDA_PATH)/targets/x86_64-linux/lib)),$(CUDA_PATH)/lib64)
XELATEX ?= xelatex

ifeq ($(USE_CUDA),1)
NVCC ?= $(if $(wildcard $(CUDA_PATH)/bin/nvcc),$(CUDA_PATH)/bin/nvcc,nvcc)
CUDA_HOST_CXX ?= $(shell bash scripts/select_cuda_host_compiler.sh "$(NVCC)" "$(CXX)")
CUDA_HOST_COMPAT := $(abspath src/cuda/CudaHostCompat.h)
ifeq ($(CXX_ORIGIN),default)
CXX := $(CUDA_HOST_CXX)
endif
endif

ifeq ($(USE_MPI),1)
CXX := mpicxx
endif

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ARCH_FLAGS ?= $(shell bash scripts/detect_arch_flags.sh "$(CXX)")
GIT_DESCRIBE ?= $(shell git describe --tags --always --dirty --long 2>/dev/null || echo unknown)

OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
CXXFLAGS += -DMBS_GIT_DESCRIBE=\"$(GIT_DESCRIBE)\"
LDFLAGS ?= -lm -lgomp
DEPFLAGS = -MMD -MP

ifeq ($(USE_CUDA),1)
CXXFLAGS += -DUSE_CUDA -I$(CUDA_PATH)/include
LDFLAGS += -L$(CUDA_LIBDIR) -lcufft -lcudart
endif
ifeq ($(USE_MPI),1)
CXXFLAGS += -DUSE_MPI
endif

SRC_DIR = src
INCLUDES = -I$(SRC_DIR) -I$(SRC_DIR)/math -I$(SRC_DIR)/handler \
           -I$(SRC_DIR)/common -I$(SRC_DIR)/geometry \
           -I$(SRC_DIR)/geometry/intrinsic -I$(SRC_DIR)/geometry/sse \
           -I$(SRC_DIR)/particle -I$(SRC_DIR)/scattering \
           -I$(SRC_DIR)/tracer -I$(SRC_DIR)/splitting \
           -I$(SRC_DIR)/cuda -Isrc/bigint

SOURCES = $(shell find $(SRC_DIR) -name '*.cpp') \
          $(shell find src/bigint -name '*.cc' 2>/dev/null)
ifeq ($(USE_CUDA),1)
NVCCFLAGS ?= -O3 -std=c++11 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700
override NVCCFLAGS += -ccbin $(CUDA_HOST_CXX) -include $(CUDA_HOST_COMPAT)
GPU_PRECISION ?= double
GPU_FAST_MATH ?= 0
GPU_ARCH ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.' | sed 's/[^0-9].*//')
GPU_ARCH := $(if $(GPU_ARCH),$(GPU_ARCH),86)
ifeq ($(filter float fp32 double fp64,$(GPU_PRECISION)),)
$(error GPU_PRECISION must be one of: fp32, float, fp64, double)
endif
ifneq ($(filter 0 1,$(GPU_FAST_MATH)),$(GPU_FAST_MATH))
$(error GPU_FAST_MATH must be 0 or 1)
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
override NVCCFLAGS += $(GPU_PRECISION_DEFINES)
CXXFLAGS += $(GPU_PRECISION_DEFINES) -DMBS_GPU_BUILD_ARCH=$(GPU_ARCH)
ifeq ($(GPU_FAST_MATH),1)
override NVCCFLAGS += --use_fast_math
CXXFLAGS += -DMBS_GPU_FAST_MATH
endif
override NVCCFLAGS += -arch=sm_$(GPU_ARCH) -Xcompiler -fopenmp
SOURCES += $(shell find $(SRC_DIR) -name '*.cu')
endif
ifeq ($(USE_CUDA),1)
CUDA_FAST_SUFFIX := $(if $(filter 1,$(GPU_FAST_MATH)),_fast,)
CUDA_MPI_SUFFIX := $(if $(filter 1,$(USE_MPI)),_mpi,)
ROOT_CUDA_OBJDIR ?= build/root_cuda/$(GPU_PRECISION_LABEL)$(CUDA_FAST_SUFFIX)$(CUDA_MPI_SUFFIX)_sm$(GPU_ARCH)/obj
OBJECTS = $(addprefix $(ROOT_CUDA_OBJDIR)/,$(SOURCES))
else
OBJECTS = $(SOURCES)
endif
OBJECTS := $(OBJECTS:.cpp=.o)
OBJECTS := $(OBJECTS:.cc=.o)
OBJECTS := $(OBJECTS:.cu=.o)
DEPS = $(OBJECTS:.o=.d)

TARGET = bin/mbs_po
TARGET_FLOAT = bin/mbs_po_float
TARGET_FLOAT_FAST = bin/mbs_po_float_fast
TARGET_DOUBLE = bin/mbs_po_double
TARGET_DOUBLE_FAST = bin/mbs_po_double_fast
FFT_PROBE = bin/fft_aperture_probe
GPU_TRACE_PROBE = bin/gpu_trace_projection_probe
GPU_QUATERNION_PROBE = bin/gpu_quaternion_rotation_probe

all: $(TARGET)

ifeq ($(USE_CUDA),1)
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
endif

cpu:
	$(MAKE) -C cpu all

gpu:
	$(MAKE) -C gpu all

gpu_float:
	$(MAKE) -C gpu float

gpu_float_fast:
	$(MAKE) -C gpu float_fast

gpu_double:
	$(MAKE) -C gpu double

gpu_double_fast:
	$(MAKE) -C gpu double_fast

gpu_fp32: gpu_float

gpu_fp64: gpu_double

gpu_fp32_fast: gpu_float_fast

gpu_fp64_fast: gpu_double_fast

gpu_double_fast_mpi:
	$(MAKE) -C gpu double_fast_mpi

split: cpu gpu_float_fast

docs:
	@command -v $(XELATEX) >/dev/null 2>&1 || { \
		echo "XeLaTeX not found: $(XELATEX)" >&2; \
		echo "Install XeLaTeX and the packages/fonts imported by docs/MANUAL*.tex." >&2; \
		exit 1; \
	}
	@command -v fc-match >/dev/null 2>&1 || { \
		echo "fontconfig not found: fc-match is required to verify manual fonts." >&2; \
		exit 1; \
	}
	@fc-match 'CMU Serif' | grep -qi 'CMU Serif' || { \
		echo "CMU Serif not found. On Ubuntu: sudo apt install fonts-cmu" >&2; \
		exit 1; \
	}
	cd docs && $(XELATEX) -interaction=nonstopmode -halt-on-error MANUAL.tex
	cd docs && $(XELATEX) -interaction=nonstopmode -halt-on-error MANUAL.tex
	cd docs && $(XELATEX) -interaction=nonstopmode -halt-on-error MANUAL_RU.tex
	cd docs && $(XELATEX) -interaction=nonstopmode -halt-on-error MANUAL_RU.tex
	@if grep -Eq 'Overfull|LaTeX Error|Undefined control sequence|Emergency stop' \
		docs/MANUAL.log docs/MANUAL_RU.log; then \
		echo "Manual layout/error diagnostics found in docs/MANUAL*.log." >&2; \
		echo "Fix every Overfull box and LaTeX error before release." >&2; \
		exit 1; \
	fi

test_cli:
	tests/run_cli_tests.sh

test_release: cpu test_cli
	MBS=cpu/bin/mbs_po_mpi tests/run_release_cli_matrix.sh

test_adaptive: cpu
	MBS=cpu/bin/mbs_po_mpi tests/run_adaptive_tests.sh

test_regression: cpu
	MBS=cpu/bin/mbs_po_mpi SKIP_BUILD=1 tests/run_tests.sh

test_extinction: cpu
	MBS=cpu/bin/mbs_po_mpi tests/run_extinction_reference.sh

test_so3: cpu
	MBS=cpu/bin/mbs_po_mpi tests/run_so3_symmetry_test.sh

test_poles: cpu
	MBS_BIN=cpu/bin/mbs_po_mpi tests/run_pole_stability_test.sh

test_beam_topology: cpu
	MBS_BIN=cpu/bin/mbs_po_mpi tests/run_beam_topology_stability_test.sh

test_coherence: cpu
	MBS_BIN=cpu/bin/mbs_po_mpi tests/run_coherent_sum_test.sh

test_forward_depth:
	tests/run_forward_depth_clipping_test.sh

test_warnings:
	tests/run_warning_gate.sh

test_cuda_build:
	$(MAKE) -C gpu GPU_PRECISION=double GPU_FAST_MATH=0 \
		TARGET=bin/mbs_po_gpu_double all

test_cuda_profiles:
	tests/run_cuda_precision_profile_test.sh

test_cuda: test_cuda_build cpu
	MBS_CPU=$(CURDIR)/cpu/bin/mbs_po_mpi \
		MBS_GPU=$(CURDIR)/gpu/bin/mbs_po_gpu_double \
		tests/run_cuda_release_gate.sh

test_cuda_if_available:
	@if command -v nvcc >/dev/null 2>&1 \
		&& command -v nvidia-smi >/dev/null 2>&1 \
		&& nvidia-smi -L >/dev/null 2>&1; then \
		$(MAKE) test_cuda; \
	else \
		echo "CUDA release gate: SKIP (nvcc or an NVIDIA GPU is unavailable)"; \
	fi

test_sanitize:
	$(MAKE) -C cpu TARGET=bin/mbs_po_mpi_sanitize \
		OBJDIR=build/sanitize/obj \
		OPT_FLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' \
		LDFLAGS='-lm -lgomp -fsanitize=address,undefined' -j
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1:abort_on_error=1 \
		UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
		MBS=$(CURDIR)/cpu/bin/mbs_po_mpi_sanitize \
		MBS_RELEASE_TIMEOUT_SECONDS=90 tests/run_release_cli_matrix.sh

test: test_release test_adaptive test_regression test_extinction test_so3 test_poles test_beam_topology test_coherence test_forward_depth test_warnings test_cuda_profiles test_cuda_if_available

$(TARGET): $(OBJECTS)
	@mkdir -p bin
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "USE_CUDA: $(USE_CUDA)"
	@echo "USE_MPI: $(USE_MPI)"
ifeq ($(USE_CUDA),1)
	@echo "GPU_PRECISION: $(GPU_PRECISION_LABEL) ($(GPU_PRECISION_CANON) storage)"
	@echo "GPU_CRITICAL_PHASE: fp64"
	@echo "GPU_FAST_MATH: $(GPU_FAST_MATH)"
	@echo "GPU_ARCH: sm_$(GPU_ARCH)"
	@echo "CUDA_OBJECTS: $(ROOT_CUDA_OBJDIR)"
	@echo "CUDA_HOST_CXX: $(CUDA_HOST_CXX)"
endif
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

ifeq ($(USE_CUDA),1)
$(ROOT_CUDA_OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

$(ROOT_CUDA_OBJDIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

$(ROOT_CUDA_OBJDIR)/%.o: %.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -DUSE_CUDA -I$(CUDA_PATH)/include $(INCLUDES) -c $< -o $@
else
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -DUSE_CUDA -I$(CUDA_PATH)/include $(INCLUDES) -c $< -o $@
endif

clean:
	find $(SRC_DIR) -name '*.o' -delete
	find $(SRC_DIR) -name '*.d' -delete
	find src/bigint -name '*.o' -delete 2>/dev/null; true
	find src/bigint -name '*.d' -delete 2>/dev/null; true
	rm -rf build/root_cuda
	rm -f $(TARGET) $(TARGET_FLOAT) $(TARGET_FLOAT_FAST) $(TARGET_DOUBLE) $(TARGET_DOUBLE_FAST) $(FFT_PROBE) $(GPU_TRACE_PROBE)

clean_cuda_objects:
	find $(SRC_DIR)/cuda -name '*.o' -delete
	find $(SRC_DIR)/cuda -name '*.d' -delete

cuda_float:
	$(MAKE) clean_cuda_objects
	$(MAKE) USE_CUDA=1 GPU_PRECISION=float TARGET=$(TARGET_FLOAT) all
	$(MAKE) clean_cuda_objects

cuda_float_fast:
	$(MAKE) clean_cuda_objects
	$(MAKE) USE_CUDA=1 GPU_PRECISION=float GPU_FAST_MATH=1 TARGET=$(TARGET_FLOAT_FAST) all
	$(MAKE) clean_cuda_objects

cuda_double:
	$(MAKE) clean_cuda_objects
	$(MAKE) USE_CUDA=1 GPU_PRECISION=double GPU_FAST_MATH=0 TARGET=$(TARGET_DOUBLE) all
	$(MAKE) clean_cuda_objects

cuda_double_fast:
	$(MAKE) clean_cuda_objects
	$(MAKE) USE_CUDA=1 GPU_PRECISION=double GPU_FAST_MATH=1 TARGET=$(TARGET_DOUBLE_FAST) all
	$(MAKE) clean_cuda_objects

cuda_variants:
	$(MAKE) cuda_float
	$(MAKE) cuda_double
	$(MAKE) cuda_float_fast
	$(MAKE) cuda_double_fast

ifeq ($(USE_CUDA),1)
fft_probe: $(FFT_PROBE)
gpu_trace_probe: $(GPU_TRACE_PROBE)
gpu_quaternion_probe: $(GPU_QUATERNION_PROBE)

$(FFT_PROBE): tools/fft_aperture_probe.cu
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -I$(CUDA_PATH)/include $< -o $@ -L$(CUDA_PATH)/lib64 -lcufft -lcudart

$(GPU_TRACE_PROBE): tools/gpu_trace_projection_probe.cu
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -I$(CUDA_PATH)/include $< -o $@ -L$(CUDA_PATH)/lib64 -lcudart

$(GPU_QUATERNION_PROBE): tools/gpu_quaternion_rotation_probe.cu
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -I$(CUDA_PATH)/include $< -o $@ -L$(CUDA_PATH)/lib64 -lcudart
else
fft_probe:
	@echo "fft_probe requires USE_CUDA=1"
	@false
gpu_trace_probe:
	@echo "gpu_trace_probe requires USE_CUDA=1"
	@false
gpu_quaternion_probe:
	@echo "gpu_quaternion_probe requires USE_CUDA=1"
	@false
endif

.PHONY: all cuda_check cpu gpu gpu_float gpu_float_fast gpu_double gpu_double_fast \
	gpu_fp32 gpu_fp64 gpu_fp32_fast gpu_fp64_fast split docs \
	test test_cli test_release test_adaptive test_regression test_so3 test_poles \
	test_beam_topology test_coherence test_forward_depth test_warnings test_cuda_build test_cuda \
	test_cuda_profiles test_cuda_if_available test_sanitize clean \
	clean_cuda_objects cuda_float cuda_float_fast cuda_double cuda_double_fast cuda_variants \
	fft_probe gpu_trace_probe gpu_quaternion_probe

-include $(DEPS)
