CXX ?= g++
CUDA_PATH ?= /usr/local/cuda

ifeq ($(USE_MPI),1)
CXX := mpicxx
endif

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ARCH_FLAGS ?= $(shell bash scripts/detect_arch_flags.sh "$(CXX)")

OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
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
NVCC ?= nvcc
NVCCFLAGS ?= -O3 -std=c++11 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700
GPU_PRECISION ?= double
GPU_FAST_MATH ?= 0
ifeq ($(GPU_PRECISION),float)
override NVCCFLAGS += -DMBS_GPU_FLOAT
endif
ifeq ($(GPU_FAST_MATH),1)
override NVCCFLAGS += --use_fast_math
endif
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

$(TARGET): $(OBJECTS)
	@mkdir -p bin
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "USE_CUDA: $(USE_CUDA)"
	@echo "USE_MPI: $(USE_MPI)"
ifeq ($(USE_CUDA),1)
	@echo "GPU_PRECISION: $(GPU_PRECISION)"
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

.PHONY: all clean clean_cuda_objects cuda_float cuda_float_fast cuda_double_fast cuda_variants fft_probe gpu_trace_probe

-include $(DEPS)
