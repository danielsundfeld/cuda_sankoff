BIN_DIR     = ./bin
BIN_CPU_SANK = $(BIN_DIR)/sankoff
BIN_GPU_SANK = $(BIN_DIR)/cuda_sankoff

TARGET      = $(BIN_CPU_SANK) $(BIN_GPU_SANK)

SRC_DIR     = ./src
INC_DIR     = ./src -I./ViennaRNA-2.3.3/H/ -I./ViennaRNA-2.3.3/src/
OBJ_DIR     = ./obj

CPPFLAGS   += -W -Wall -fopenmp
LDFLAGS    += -fopenmp -lstdc++ -lm

GPUFLAGS   += -dc -arch sm_52
GPULDFLAGS   += -link -arch sm_52

VIENNA_LIB  = ./ViennaRNA-2.3.3/src/ViennaRNA/libRNA.a

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

COMMON_SRCS += \
    $(SRC_DIR)/Backtrace.cpp \
    $(SRC_DIR)/DPMatrix.cpp \
    $(SRC_DIR)/read_fasta.cpp \
    $(SRC_DIR)/sankoff_args.cpp \
    $(SRC_DIR)/Sequences.cpp \
    $(SRC_DIR)/TimeCounter.cpp \
    $(SRC_DIR)/bp_probs.cpp \

CPU_SANK_SRCS += \
    $(SRC_DIR)/sankoff_cpu_main.cpp \
    $(SRC_DIR)/Sankoff.cpp \

GPU_CPP_SANK_SRCS += \
    $(SRC_DIR)/sankoff_gpu_main.cpp \

# It is faster to include all .cu in the same .cu
GPU_CUDA_SANK_SRCS += \
    $(SRC_DIR)/Sankoff_GPU.cu \

INC_PATH += \
    -I$(INC_DIR)

CPPFLAGS += $(INC_PATH)

COMMON_OBJS = $(COMMON_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
CPU_SANK_OBJS = $(CPU_SANK_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
GPU_CPP_SANK_OBJS = $(GPU_CPP_SANK_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
GPU_CUDA_SANK_OBJS = $(GPU_CUDA_SANK_SRCS:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)

all:	$(TARGET) $(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(COMMON_OBJS): 	| $(OBJ_DIR)
$(CPU_SANK_OBJS):	| $(OBJ_DIR)
$(GPU_CUDA_SANK_OBJS):	| $(OBJ_DIR)
$(GPU_CPP_SANK_OBJS):	| $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	nvcc $(GPUFLAGS) -c -o $@ $<

$(BIN_CPU_SANK):	$(COMMON_OBJS) $(CPU_SANK_OBJS) | $(BIN_DIR)
	$(CXX) $^ $(VIENNA_LIB) -o $@ $(LDFLAGS)
$(BIN_GPU_SANK):	$(COMMON_OBJS) $(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS) | $(BIN_DIR)
	nvcc $^ $(VIENNA_LIB) -o $@ $(GPULDFLAGS)

clean:
	rm -f $(TARGET) $(COMMON_OBJS) $(CPU_SANK_OBJS) $(GPU_CUDA_SANK_OBJS) $(GPU_CPP_SANK_OBJS)
