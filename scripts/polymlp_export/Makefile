# Builds polymlp_export against the C++ core of pypolymlp (BSD-3-Clause), e.g.
#   make POLYMLP_SRC=$(python -c "import pypolymlp, os; print(os.path.dirname(pypolymlp.__file__))")/cxx/src
# The gtinv tables (*.bin) of pypolymlp are copied next to the executable, where the library looks for them.
POLYMLP_SRC ?= $(error set POLYMLP_SRC to the cxx/src directory of pypolymlp)
CXX ?= c++
CXXFLAGS ?= -std=c++17 -O2 -w
BUILD ?= build
SRC := $(filter-out %boost.cpp,$(wildcard $(POLYMLP_SRC)/polymlp/*.cpp))
OBJ := $(addprefix $(BUILD)/,$(notdir $(SRC:.cpp=.o))) $(BUILD)/polymlp_export.o

polymlp_export: $(OBJ)
	$(CXX) -o $@ $(OBJ)
	cp $(POLYMLP_SRC)/polymlp/*.bin .

$(BUILD)/%.o: $(POLYMLP_SRC)/polymlp/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(POLYMLP_SRC)/polymlp -I$(POLYMLP_SRC) -c $< -o $@

$(BUILD)/polymlp_export.o: polymlp_export.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(POLYMLP_SRC)/polymlp -I$(POLYMLP_SRC) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD) polymlp_export *.bin

.PHONY: clean
