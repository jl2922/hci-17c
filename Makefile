# Default options.
CXX := mpic++
CXXFLAGS := -std=c++11 -Wall -Wextra -O3 -fopenmp
LDLIBS := -lboost_mpi -lboost_serialization
SRC_DIR := src
OBJ_DIR := build
EXE := hci.out
TEST_EXE := hci_test.out

# Host specific configurations.
HOSTNAME := $(shell hostname)
ifneq ($(findstring bridges, $(HOSTNAME)),)
include Makefile.config.bridges
endif
ifneq ($(findstring dft, $(HOSTNAME)),)
include Makefile.config.dft
endif

# Load Makefile.config if exists.
CONFIG_FILE := Makefile.config
ifneq ($(wildcard $(CONFIG_FILE)),)
include $(CONFIG_FILE)
endif

# Sources and intermediate objects.
SRCS := $(shell find $(SRC_DIR) \
		! -name "main.cc" ! -name "*_test.cc" -name "*.cc")
TESTS := $(shell find $(SRC_DIR) -name "*_test.cc")
HEADERS := $(shell find $(SRC_DIR) -name "*.h")
MAIN := $(SRC_DIR)/main.cc
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
TEST_OBJS := $(TESTS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

# GTest related.
GTEST_DIR := gtest/googletest
GTEST_CXXFLAGS := $(CXXFLAGS) -isystem $(GTEST_DIR)/include -pthread
GTEST_HEADERS := $(GTEST_DIR)/include/gtest/*.h \
		$(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS := $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h \
		$(GTEST_HEADERS)
GTEST_MAIN := $(OBJ_DIR)/gtest_main.a

.PHONY: all test clean

all: $(EXE)

test: $(TEST_EXE)
	./$(TEST_EXE)

clean:
	rm -rf $(OBJ_DIR)
	rm -f ./$(EXE)
	rm -f ./$(TEST_EXE)

# Main program.

$(EXE): $(OBJS) $(MAIN) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(MAIN) $(OBJS) -o $(EXE) $(LDLIBS)
	
$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	mkdir -p $(@D) && $(CXX) $(CXXFLAGS) -c $< -o $@
	
# Test specific.

$(TEST_EXE): $(TEST_OBJS) $(OBJS) $(GTEST_MAIN) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(TEST_OBJS) $(OBJS) $(GTEST_MAIN) \
			-o $(TEST_EXE) $(LDLIBS) -lpthread

$(TEST_OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	mkdir -p $(@D) && $(CXX) $(GTEST_CXXFLAGS) -c $< -o $@

$(GTEST_MAIN): $(OBJ_DIR)/gtest-all.o $(OBJ_DIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

$(OBJ_DIR)/gtest-all.o: $(GTEST_SRCS)
	mkdir -p $(@D) && $(CXX) -I$(GTEST_DIR) $(GTEST_CXXFLAGS) -c \
			$(GTEST_DIR)/src/gtest-all.cc -o $@

$(OBJ_DIR)/gtest_main.o: $(GTEST_SRCS)
	mkdir -p $(@D) && $(CXX) -I$(GTEST_DIR) $(GTEST_CXXFLAGS) -c \
			$(GTEST_DIR)/src/gtest_main.cc -o $@
