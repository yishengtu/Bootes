CC := nvc++
C_CUDA :=nvcc
CFLAGS := -I /home/vhatz6ed/hdf5/include/ -L /home/vhatz6ed/hdf5/lib/ -g -lhdf5 -lhdf5_cpp -cuda -gpu=nordc#-acc -Minfo=accel -gpu=managed,lineinfo,ptxinfo
CUDA_FLAGS := -ccbin=nvc++ 
EXE_DIR := bin/
EXECUTABLE := $(EXE_DIR)bootes.out

SRC_FILES := $(wildcard src/algorithm/*.cpp) \
	     $(wildcard src/algorithm/eos/*.cpp) \
	     $(wildcard src/algorithm/util/*.cpp) \
	     $(wildcard src/algorithm/boundary_condition/*.cpp) \
	     $(wildcard src/algorithm/gravity/*.cpp) \
	     $(wildcard src/algorithm/hydro/*.cpp) \
	     $(wildcard src/algorithm/hydro/srcterm/*cpp) \
	     $(wildcard src/algorithm/mesh/*.cpp) \
	     $(wildcard src/algorithm/timeadvance/*.cpp) \
	     $(wildcard src/algorithm/reconstruct/*.cpp) \
	     $(wildcard src/algorithm/time_step/*.cpp) \
	     $(wildcard src/algorithm/dust/*.cpp) \
	     $(wildcard src/algorithm/dust/graingrowth/*.cpp) \
	     $(wildcard src/algorithm/boundary_condition/dust/*.cpp) \
	     $(wildcard src/main.cpp)

OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR), $(notdir $(SRC_FILES:.cpp=.o))) obj/coagulation.o
#OBJ_FILES := $(addprefix $(OBJ_DIR), $(notdir $(SRC_FILES:.cpp=.o)), obj/coagulation.o)
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

.PHONY : all dirs clean

all : dirs $(EXECUTABLE)

objs: dirs $(OBJ_FILES)

dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXECUTABLE) : $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $(OBJ_FILES)

$(OBJ_DIR)%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

obj/coagulation.o : src/algorithm/dust/graingrowth/coagulation.cu
	$(C_CUDA) $(CUDA_FLAGS) -c $< -o $@


clean:
	rm obj/*
	rm $(EXECUTABLE)


