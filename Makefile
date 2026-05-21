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
ifeq ($(GPU_PRECISION),float)
NVCCFLAGS += -DMBS_GPU_FLOAT
endif
SOURCES += $(shell find $(SRC_DIR) -name '*.cu')
endif
OBJECTS = $(SOURCES:.cpp=.o)
OBJECTS := $(OBJECTS:.cc=.o)
OBJECTS := $(OBJECTS:.cu=.o)
DEPS = $(OBJECTS:.o=.d)

TARGET = bin/mbs_po
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
	rm -f $(TARGET) $(FFT_PROBE) $(GPU_TRACE_PROBE)

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

.PHONY: all clean fft_probe gpu_trace_probe

-include $(DEPS)
