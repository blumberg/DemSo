NVCC      = nvcc

TARGET    = DemSo	
SRC_DIR   = src
OBJ_DIR   = obj

NVCCFLAGS	:= -arch=sm_11 -I includes
LIBS		:= -lglut -lpthread -lGL

CPP_FILES = $(wildcard $(SRC_DIR)/*.cpp)
CU_FILES  = $(wildcard $(SRC_DIR)/*.cu)

H_FILES   = $(wildcard $(SRC_DIR)/*.h)
H_FILES  += $(wildcard $(SRC_DIR)/*.hpp)
CUH_FILES = $(wildcard $(SRC_DIR)/*.cuh)

OBJ_FILES = $(addprefix $(OBJ_DIR)/,$(notdir $(CPP_FILES:.cpp=.o)))
CUO_FILES = $(addprefix $(OBJ_DIR)/,$(notdir $(CU_FILES:.cu=.cu.o)))

OBJS =  $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(CPP_FILES)))
OBJS += $(patsubst %.cu,$(OBJ_DIR)/%.cu.o,$(notdir $(CU_FILES)))

$(TARGET) : $(OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $(OBJS) $(LIBS)

$(OBJ_DIR)/%.cu.o : $(SRC_DIR)/%.cu $(CUH_FILES)
	$(NVCC) $(NVCCFLAGS) -c -o $@ $< $(LIBS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(H_FILES)
	$(NVCC) $(NVCCFLAGS) -c -o $@ $< $(LIBS)

.PRONY: all
all: $(TARGET)

.PRONY: distclean	
distclean:
	rm -f $(OBJS) $(TARGET) src/*~ *~

.PRONY: clean	
clean:
	rm -f $(TARGET)

parser: src/parser.cpp src/DEMSimulation.cpp
	$(NVCC) src/parser.cpp src/DEMParticles.cpp src/DEMSimulation.cpp -o parser -I./includes $(NVCCFLAGS) 
