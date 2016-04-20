BIN_DIR     = ./bin
BIN_CPU_SANK = $(BIN_DIR)/sankoff
BIN_GPU_SANK = $(BIN_DIR)/cuda_sankoff

TARGET      = $(BIN_CPU_SANK) $(BIN_GPU_SANK)

SRC_DIR     = ./src
INC_DIR     = ./src
OBJ_DIR     = ./obj
CPPFLAGS   += -W -Wall -fopenmp
LDFLAGS    += -fopenmp -lstdc++ -lm

ifndef DEBUG
    OPTIMIZE = yes
    LDFLAGS += -s
else
    CPPFLAGS += -g
    GPUFLAGS += -g -G
endif

ifdef OPTIMIZE
    CPPFLAGS += -O3
    GPUFLAGS += -O3
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
    $(SRC_DIR)/DPMatrix.cpp \
    $(SRC_DIR)/Foldalign.cpp \
    $(SRC_DIR)/Sankoff.cpp \

GPU_CPP_SANK_SRCS += \
    $(SRC_DIR)/sankoff_gpu_main.cpp \

GPU_CUDA_SANK_SRCS += \
    $(SRC_DIR)/DPMatrix_GPU.cu \
    $(SRC_DIR)/Sankoff_GPU.cu \

INC_PATH += \
    -I$(INC_DIR)

CPPFLAGS += $(INC_PATH)

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
$(BIN_GPU_SANK):	$(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS) | $(BIN_DIR)
	nvcc $^ -o $@ $(GPULDFLAGS)

clean:
	rm -f $(TARGET) $(CPU_SANK_OBJS) $(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS)
