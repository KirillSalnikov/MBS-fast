CXX_ORIGIN := $(origin CXX)
CXX ?= g++
CUDA_PATH ?= /usr/local/cuda
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
LDFLAGS += -L$(CUDA_PATH)/lib64 -lcufft -lcudart
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
ifeq ($(GPU_PRECISION),float)
override NVCCFLAGS += -DMBS_GPU_FLOAT
CXXFLAGS += -DMBS_GPU_FLOAT
endif
ifeq ($(GPU_FAST_MATH),1)
override NVCCFLAGS += --use_fast_math
CXXFLAGS += -DMBS_GPU_FAST_MATH
endif
override NVCCFLAGS += -Xcompiler -fopenmp
SOURCES += $(shell find $(SRC_DIR) -name '*.cu')
endif
OBJECTS = $(SOURCES:.cpp=.o)
OBJECTS := $(OBJECTS:.cc=.o)
OBJECTS := $(OBJECTS:.cu=.o)
DEPS = $(OBJECTS:.o=.d)

TARGET = bin/mbs_po
TARGET_FLOAT = bin/mbs_po_float
TARGET_FLOAT_FAST = bin/mbs_po_float_fast
TARGET_DOUBLE_FAST = bin/mbs_po_double_fast
FFT_PROBE = bin/fft_aperture_probe
GPU_TRACE_PROBE = bin/gpu_trace_projection_probe

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

gpu_double_fast:
	$(MAKE) -C gpu double_fast

gpu_double_fast_mpi:
	$(MAKE) -C gpu double_fast_mpi

split: cpu gpu_float_fast

docs:
	@command -v $(XELATEX) >/dev/null 2>&1 || { \
		echo "XeLaTeX not found: $(XELATEX)" >&2; \
		echo "Install XeLaTeX and the packages/fonts imported by docs/MANUAL*.tex." >&2; \
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

test_sanitize:
	$(MAKE) -C cpu TARGET=bin/mbs_po_mpi_sanitize \
		OBJDIR=build/sanitize/obj \
		OPT_FLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' \
		LDFLAGS='-lm -lgomp -fsanitize=address,undefined' -j
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1:abort_on_error=1 \
		UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
		MBS=$(CURDIR)/cpu/bin/mbs_po_mpi_sanitize \
		MBS_RELEASE_TIMEOUT_SECONDS=90 tests/run_release_cli_matrix.sh

test: test_release test_adaptive test_regression test_extinction

$(TARGET): $(OBJECTS)
	@mkdir -p bin
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "USE_CUDA: $(USE_CUDA)"
	@echo "USE_MPI: $(USE_MPI)"
ifeq ($(USE_CUDA),1)
	@echo "GPU_PRECISION: $(GPU_PRECISION)"
	@echo "CUDA_HOST_CXX: $(CUDA_HOST_CXX)"
endif
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -DUSE_CUDA -I$(CUDA_PATH)/include $(INCLUDES) -c $< -o $@

clean:
	find $(SRC_DIR) -name '*.o' -delete
	find $(SRC_DIR) -name '*.d' -delete
	find src/bigint -name '*.o' -delete 2>/dev/null; true
	find src/bigint -name '*.d' -delete 2>/dev/null; true
	rm -f $(TARGET) $(TARGET_FLOAT) $(TARGET_FLOAT_FAST) $(TARGET_DOUBLE_FAST) $(FFT_PROBE) $(GPU_TRACE_PROBE)

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

cuda_double_fast:
	$(MAKE) clean_cuda_objects
	$(MAKE) USE_CUDA=1 GPU_FAST_MATH=1 TARGET=$(TARGET_DOUBLE_FAST) all
	$(MAKE) clean_cuda_objects

cuda_variants:
	$(MAKE) cuda_float
	$(MAKE) cuda_float_fast
	$(MAKE) cuda_double_fast

ifeq ($(USE_CUDA),1)
fft_probe: $(FFT_PROBE)
gpu_trace_probe: $(GPU_TRACE_PROBE)

$(FFT_PROBE): tools/fft_aperture_probe.cu
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -I$(CUDA_PATH)/include $< -o $@ -L$(CUDA_PATH)/lib64 -lcufft -lcudart

$(GPU_TRACE_PROBE): tools/gpu_trace_projection_probe.cu
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -I$(CUDA_PATH)/include $< -o $@ -L$(CUDA_PATH)/lib64 -lcudart
else
fft_probe:
	@echo "fft_probe requires USE_CUDA=1"
	@false
gpu_trace_probe:
	@echo "gpu_trace_probe requires USE_CUDA=1"
	@false
endif

.PHONY: all cuda_check cpu gpu gpu_float gpu_float_fast gpu_double_fast split docs \
	test test_cli test_release test_adaptive test_regression test_sanitize clean \
	clean_cuda_objects cuda_float cuda_float_fast cuda_double_fast cuda_variants \
	fft_probe gpu_trace_probe

-include $(DEPS)
