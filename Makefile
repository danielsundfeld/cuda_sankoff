# foldalign Makefile
BIN_DIR     = ./bin
BIN         = $(BIN_DIR)/foldalign

TARGET      = $(BIN)

SRC_DIR     = ./src
INC_DIR     = ./src
OBJ_DIR     = ./obj
CPPFLAGS   += -W -Wall -pthread
LDFLAGS    += -pthread -lstdc++ -lm

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
    $(SRC_DIR)/foldalign.cpp \

INC_PATH += \
    -I$(INC_DIR) \
    -I/usr/include \

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
