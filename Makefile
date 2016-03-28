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
BIN         = $(BIN_DIR)/sankoff

TARGET      = $(BIN)

SRC_DIR     = ./src
INC_DIR     = ./src
OBJ_DIR     = ./obj
CPPFLAGS   += -W -Wall -fopenmp $(CPPSTD)
LDFLAGS    += -fopenmp -lstdc++ -lm

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

CPP_SRCS += \
    $(SRC_DIR)/main.cpp \
    $(SRC_DIR)/Cost.cpp \
    $(SRC_DIR)/DPMatrix.cpp \
    $(SRC_DIR)/Foldalign.cpp \
    $(SRC_DIR)/Sankoff.cpp \

INC_PATH += \
    -I$(INC_DIR) \

CPPFLAGS += \
    $(INC_PATH) \

OBJS = $(CPP_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

all:	$(TARGET)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<

$(OBJS):	| $(OBJ_DIR)

$(BIN):	$(OBJS) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET) $(OBJS)
