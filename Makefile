CXX ?= g++

CPU_MODEL := $(shell lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)
ARCH_FLAGS ?= $(shell bash scripts/detect_arch_flags.sh "$(CXX)")

OPT_FLAGS ?= -O3 -funroll-loops
CXXFLAGS ?= $(OPT_FLAGS) $(ARCH_FLAGS) -std=gnu++11 -fopenmp
LDFLAGS ?= -lm -lgomp
DEPFLAGS = -MMD -MP

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
DEPS = $(OBJECTS:.o=.d)

TARGET = bin/mbs_po

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p bin
	@echo "CPU: $(CPU_MODEL)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	find $(SRC_DIR) -name '*.o' -delete
	find $(SRC_DIR) -name '*.d' -delete
	find src/bigint -name '*.o' -delete 2>/dev/null; true
	find src/bigint -name '*.d' -delete 2>/dev/null; true
	rm -f $(TARGET)

.PHONY: all clean

-include $(DEPS)
