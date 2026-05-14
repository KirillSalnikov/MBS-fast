CXX ?= g++

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ifneq (,$(findstring EPYC 7H12,$(CPU_MODEL)))
ARCH_FLAGS ?= -march=znver2 -mtune=znver2
else
ARCH_FLAGS ?= -march=native -mtune=native
endif

OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
LDFLAGS ?= -lm -lgomp

SRC_DIR = src
INCLUDES = -I$(SRC_DIR) -I$(SRC_DIR)/math -I$(SRC_DIR)/handler \
           -I$(SRC_DIR)/common -I$(SRC_DIR)/geometry \
           -I$(SRC_DIR)/geometry/intrinsic -I$(SRC_DIR)/geometry/sse \
           -I$(SRC_DIR)/particle -I$(SRC_DIR)/scattering \
           -I$(SRC_DIR)/tracer -I$(SRC_DIR)/splitting \
           -Isrc/bigint

SOURCES = $(shell find $(SRC_DIR) -name '*.cpp') \
          $(shell find src/bigint -name '*.cc' 2>/dev/null)
OBJECTS = $(SOURCES:.cpp=.o)
OBJECTS := $(OBJECTS:.cc=.o)

TARGET = bin/mbs_po

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p bin
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	find $(SRC_DIR) -name '*.o' -delete
	find src/bigint -name '*.o' -delete 2>/dev/null; true
	rm -f $(TARGET)

.PHONY: all clean
