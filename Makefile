# Choose -std=c++11 or -std=c++0x
CXXVERSION = $(shell $(CXX) -dumpversion | cut -b 1-3)
ifneq "$(filter g++,$(CXX))" ""
ifeq "$(CXXVERSION)" "4.6"
CPPSTD = -std=c++0x
endif
ifeq "$(CXXVERSION)" "4.4"
$(error Bad $(CXX) version $(CXXVERSION). Atomic operations are required)
endif
endif

ifeq "$(CPPSTD)" ""
CPPSTD = -std=c++11
endif

BIN_DIR     = ./bin
BIN_CPU_SANK = $(BIN_DIR)/sankoff
BIN_GPU_SANK = $(BIN_DIR)/cuda_sankoff

TARGET      = $(BIN_CPU_SANK) $(BIN_GPU_SANK)

SRC_DIR     = ./src
INC_DIR     = ./src
OBJ_DIR     = ./obj
CPPFLAGS   += -W -Wall -fopenmp $(CPPSTD)
LDFLAGS    += -fopenmp -lstdc++ -lm
CUDAFLAGS  += -lcuda -lcudart -L/usr/local/cuda/lib64/

ifndef DEBUG
    OPTIMIZE = yes
    LDFLAGS += -s
else
    CPPFLAGS += -g
endif

ifdef OPTIMIZE
    CPPFLAGS += -O3
else
    CPPFLAGS += -O0
endif

ifdef PROFILE_GENERATE
    CPPFLAGS += -fprofile-generate
    LDFLAGS  += -fprofile-generate
endif

ifdef PROFILE_USE
    CPPFLAGS += -fprofile-use
    LDFLAGS += -fprofile-use
endif

ifdef PROFILE_INFORMATION
    CPPFLAGS += -pg
    LDFLAGS += -pg
endif

CPU_SANK_SRCS += \
    $(SRC_DIR)/main.cpp \
    $(SRC_DIR)/Cost.cpp \
    $(SRC_DIR)/DPMatrix.cpp \
    $(SRC_DIR)/Foldalign.cpp \
    $(SRC_DIR)/Sankoff.cpp \

GPU_CPP_SANK_SRCS += \
    $(SRC_DIR)/sankoff_gpu_main.cpp \

GPU_CUDA_SANK_SRCS += \
    $(SRC_DIR)/Sankoff_GPU.cu \

INC_PATH += \
    -I$(INC_DIR) \

CPPFLAGS += \
    $(INC_PATH) \

CPU_SANK_OBJS = $(CPU_SANK_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
GPU_CPP_SANK_OBJS = $(GPU_CPP_SANK_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
GPU_CUDA_SANK_OBJS = $(GPU_CUDA_SANK_SRCS:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)

all:	$(TARGET) $(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(CPU_SANK_OBJS):	| $(OBJ_DIR)
$(GPU_CUDA_SANK_OBJS):	| $(OBJ_DIR)
$(GPU_CPP_SANK_OBJS):	| $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	nvcc $(GPUFLAGS) -c -o $@ $<

$(BIN_CPU_SANK):	$(CPU_SANK_OBJS) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS)
#TODO FIXME FIXME FIXME
$(BIN_GPU_SANK):	$(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS) obj/Cost.o | $(BIN_DIR)
	$(CXX) $^ -o $@ $(CUDAFLAGS) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(CPU_SANK_OBJS) $(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS)
