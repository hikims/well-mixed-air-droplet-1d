# Makefile for the well-mixed-air single-droplet model
# Repository version for manuscript reproducibility.
#
# Typical usage:
#   make all                 # build all main executables
#   make figure1             # run Figure 1 simulation and plot
#   make figures             # reproduce all main figures
#   make supplement          # run supplementary robustness checks
#   make clean
#
# Notes:
# - Model defaults are read from src/params1.h when possible.
# - Threading Building Blocks (TBB) is enabled by default.
# - On macOS, if TBB is not found automatically, use for example:
#     make TBB_PREFIX=/opt/homebrew/opt/tbb
#   or, on Intel Homebrew:
#     make TBB_PREFIX=/usr/local/opt/tbb

SHELL := /bin/bash

# -----------------------------------------------------------------------------
# Compilers and flags
# -----------------------------------------------------------------------------
CXX ?= g++
USE_TBB ?= 1
THREADS ?= 8
TBB_PREFIX ?= $(shell brew --prefix tbb 2>/dev/null)

CXXFLAGS ?= -O3 -std=c++14 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable
CPPFLAGS ?= -I$(SRC_DIR)
LDFLAGS ?=
LDLIBS ?=

ifeq ($(USE_TBB),1)
  # Older project code used USE_REAL_TBB; USE_TBB is also defined for clarity.
  CPPFLAGS += -DUSE_REAL_TBB -DUSE_TBB
  ifneq ($(strip $(TBB_PREFIX)),)
    CPPFLAGS += -I$(TBB_PREFIX)/include
    LDFLAGS += -L$(TBB_PREFIX)/lib
  endif
  LDLIBS += -ltbb
endif

# -----------------------------------------------------------------------------
# Project directories
# -----------------------------------------------------------------------------
SRC_DIR := src
BIN_DIR := bin
DATA_DIR := data/generated
FIG_MAIN_DIR := figures/main
FIG_SUPP_DIR := figures/supplement
LOG_DIR := logs
SCRIPTS_DIR := scripts
PARAMS := $(SRC_DIR)/params1.h

COMMON_DEPS := $(SRC_DIR)/model_lu_tbb.hpp $(PARAMS)

# Read selected defaults from params1.h, with safe fallbacks.
DT_FROM_PARAMS := $(shell awk '/^\#define[[:space:]]+default_dt[[:space:]]/ {gsub(/[()]/,"",$$3); print $$3}' $(PARAMS) 2>/dev/null)
SAVE_EVERY_FROM_PARAMS := $(shell awk '/^\#define[[:space:]]+default_save_every[[:space:]]/ {gsub(/[()]/,"",$$3); print $$3}' $(PARAMS) 2>/dev/null)
S0_FROM_PARAMS := $(shell awk '/^\#define[[:space:]]+S0[[:space:]]/ {gsub(/[()]/,"",$$3); print $$3}' $(PARAMS) 2>/dev/null)

DT ?= $(if $(DT_FROM_PARAMS),$(DT_FROM_PARAMS),5.0e-4)
SAVE_EVERY ?= $(if $(SAVE_EVERY_FROM_PARAMS),$(SAVE_EVERY_FROM_PARAMS),20)
S0 ?= $(if $(S0_FROM_PARAMS),$(S0_FROM_PARAMS),0.125)

# -----------------------------------------------------------------------------
# Executables
# -----------------------------------------------------------------------------
FIG1_BIN := $(BIN_DIR)/figure1_phase_diagram
FIG2_BIN := $(BIN_DIR)/figure2_trajectories
FIG45_BIN := $(BIN_DIR)/figure4_surface_tension
FIG6_BIN := $(BIN_DIR)/figure5_cover_protocol
GRID_BIN := $(BIN_DIR)/grid_aligned_validation
SUPP_DEFENSE_BIN := $(BIN_DIR)/supplement_defense_batch
SUPP_CYCLE_BIN := $(BIN_DIR)/supplement_oscillatory_cycle

MAIN_BINS := $(FIG1_BIN) $(FIG2_BIN) $(FIG45_BIN) $(FIG6_BIN)
SUPP_BINS := $(SUPP_DEFENSE_BIN)

# -----------------------------------------------------------------------------
# Figure 1: phase diagram
# -----------------------------------------------------------------------------
PHASE_SMIN ?= 0.05
PHASE_SMAX ?= 0.20
PHASE_DS ?= 0.005
PHASE_TEND ?= 120.0
FIG1_CSV ?= $(DATA_DIR)/phase_diagram_long.csv
FIG1_OUT ?= $(FIG_MAIN_DIR)/figure1_phase_diagram.png

# -----------------------------------------------------------------------------
# Figure 2 and Figure 3: trajectories and representative profiles
# Figure 3 uses the data generated for Figure 2.
# -----------------------------------------------------------------------------
FIG2_OUT ?= $(DATA_DIR)/figure2_data
FIG2_TEND ?= 200.0
FIG2_SAVE_EVERY ?= 20
FIG2_CASES ?= 0.080 0.125 0.180

# -----------------------------------------------------------------------------
# Figures 4 and 5: oscillatory and rotating surface-tension imbalance
# -----------------------------------------------------------------------------
FIG45_OUT ?= $(DATA_DIR)/figure4_5_surface_tension_data
FIG45_TEND ?= 200.0
FIG45_SAVE_EVERY ?= 20
FIG45_CASES ?= 0.080 0.125 0.180

# -----------------------------------------------------------------------------
# Figure 6: cover protocol
# -----------------------------------------------------------------------------
FIG6_OUT ?= $(DATA_DIR)/figure6_cover_protocol_data
FIG6_S0 ?= 0.125
FIG6_OPEN_TIME ?= 40.0
FIG6_CLOSE_TIME ?= 40.0
FIG6_REOPEN_TIME ?= 40.0
FIG6_SAVE_EVERY ?= 20

# -----------------------------------------------------------------------------
# Grid-aligned validation (manuscript Fig. 8, Supplementary Figs. S4-S6)
# Protocol parameters live in scripts/run_grid_aligned_validation.sh.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Supplementary robustness checks
# -----------------------------------------------------------------------------
SUPP_OUT ?= $(DATA_DIR)/supplement_robustness_data

.PHONY: help show dirs all main figures supplement \
        figure1 figure1_data figure1_plot figure1_bg \
        figure2 figure2_data figure2_plot \
        figure3 figure3_plot \
        figures4_5 figures4_5_data figures4_5_plot \
        figure6 figure6_data figure6_plot \
        grid_aligned hopf hopf_scan hopf_timestep hopf_plots \
        supplement_defense supplement_cycle \
        check-python check-names clean clean-data distclean

help:
	@echo "Targets:"
	@echo "  make all                         build all main executables"
	@echo "  make figures                     reproduce all main figures"
	@echo "  make figure1                     reproduce phase diagram"
	@echo "  make figure2                     reproduce representative trajectories"
	@echo "  make figure3                     reproduce representative profiles"
	@echo "  make figures4_5                  reproduce surface-tension figures"
	@echo "  make figure6                     reproduce cover protocol"
	@echo "  make hopf                        odd-sector Hopf roots + IMEX time-step stability"
	@echo "  make grid_aligned                grid-aligned validation (run after make hopf)"
	@echo "  make supplement                  run supplementary robustness checks"
	@echo "  make check-python                syntax-check Python scripts"
	@echo "  make check-names                 check for old file-name references"
	@echo "  make clean / make clean-data / make distclean"

show:
	@echo "CXX          = $(CXX)"
	@echo "USE_TBB      = $(USE_TBB)"
	@echo "TBB_PREFIX   = $(TBB_PREFIX)"
	@echo "THREADS      = $(THREADS)"
	@echo "DT           = $(DT)"
	@echo "SAVE_EVERY   = $(SAVE_EVERY)"
	@echo "SRC_DIR      = $(SRC_DIR)"
	@echo "BIN_DIR      = $(BIN_DIR)"
	@echo "DATA_DIR     = $(DATA_DIR)"
	@echo "FIG_MAIN_DIR = $(FIG_MAIN_DIR)"
	@echo "FIG_SUPP_DIR = $(FIG_SUPP_DIR)"
	@echo "CPPFLAGS     = $(CPPFLAGS)"
	@echo "CXXFLAGS     = $(CXXFLAGS)"
	@echo "LDFLAGS      = $(LDFLAGS)"
	@echo "LDLIBS       = $(LDLIBS)"

dirs:
	mkdir -p $(BIN_DIR) $(DATA_DIR) $(FIG_MAIN_DIR) $(FIG_SUPP_DIR) $(LOG_DIR)

all: main

main: $(MAIN_BINS)

# -----------------------------------------------------------------------------
# Build rules
# -----------------------------------------------------------------------------
$(FIG1_BIN): $(SRC_DIR)/figure1_phase_diagram.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(FIG2_BIN): $(SRC_DIR)/figure2_trajectories.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(FIG45_BIN): $(SRC_DIR)/figure4_surface_tension.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(FIG6_BIN): $(SRC_DIR)/figure5_cover_protocol.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(GRID_BIN): $(SRC_DIR)/grid_aligned_validation.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(SUPP_DEFENSE_BIN): $(SRC_DIR)/supplement_defense_batch.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(SUPP_CYCLE_BIN): $(SRC_DIR)/supplement_oscillatory_cycle.cpp $(COMMON_DEPS) | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

# -----------------------------------------------------------------------------
# Main figures
# -----------------------------------------------------------------------------
figure1_data: $(FIG1_BIN) | dirs
	$(FIG1_BIN) $(FIG1_CSV) $(PHASE_SMIN) $(PHASE_SMAX) $(PHASE_DS) $(PHASE_TEND) $(DT) $(THREADS)

figure1_plot: | dirs
	python3 $(SCRIPTS_DIR)/plot_figure1_phase_diagram.py $(FIG1_CSV) $(FIG1_OUT)

figure1: figure1_data figure1_plot

figure1_bg: $(FIG1_BIN) | dirs
	@LOG=$(LOG_DIR)/figure1_phase_diagram_$$(date +%Y%m%d_%H%M%S).log; \
	echo "Starting Figure 1 phase sweep in background"; \
	echo "Log file: $$LOG"; \
	(nohup bash -lc '$(MAKE) figure1 THREADS=$(THREADS) PHASE_SMIN=$(PHASE_SMIN) PHASE_SMAX=$(PHASE_SMAX) PHASE_DS=$(PHASE_DS) PHASE_TEND=$(PHASE_TEND)' > $$LOG 2>&1 & echo $$! > $(LOG_DIR)/figure1.pid); \
	echo "PID saved to $(LOG_DIR)/figure1.pid"

figure2_data: $(FIG2_BIN) | dirs
	$(FIG2_BIN) $(FIG2_OUT) $(FIG2_TEND) $(DT) $(FIG2_SAVE_EVERY) $(THREADS) $(FIG2_CASES)

figure2_plot: | dirs
	python3 $(SCRIPTS_DIR)/plot_figure2_representative_trajectories.py $(FIG2_OUT) $(FIG_MAIN_DIR)

figure2: figure2_data figure2_plot

figure3_plot: figure2_data | dirs
	python3 $(SCRIPTS_DIR)/plot_figure3_representative_profiles.py $(FIG2_OUT) $(FIG_MAIN_DIR) 1.0

figure3: figure3_plot

figures4_5_data: $(FIG45_BIN) | dirs
	$(FIG45_BIN) $(FIG45_OUT) $(FIG45_TEND) $(DT) $(FIG45_SAVE_EVERY) $(THREADS) $(FIG45_CASES)

figures4_5_plot: | dirs
	python3 $(SCRIPTS_DIR)/plot_figure4_5_surface_tension.py $(FIG45_OUT) $(FIG_MAIN_DIR)

figures4_5: figures4_5_data figures4_5_plot

figure6_data: $(FIG6_BIN) | dirs
	$(FIG6_BIN) $(FIG6_OUT) $(FIG6_S0) $(FIG6_OPEN_TIME) $(FIG6_CLOSE_TIME) $(FIG6_REOPEN_TIME) $(DT) $(FIG6_SAVE_EVERY) $(THREADS)

figure6_plot: | dirs
	python3 $(SCRIPTS_DIR)/plot_figure6_cover_protocol.py $(FIG6_OUT) $(FIG_MAIN_DIR) 1.0

figure6: figure6_data figure6_plot

# Grid-aligned validation. Requires 'make hopf' first for the reference CSVs.
grid_aligned: $(GRID_BIN) | dirs
	bash $(SCRIPTS_DIR)/run_grid_aligned_validation.sh

# -----------------------------------------------------------------------------
# Odd-sector Hopf analysis (manuscript Fig. 9, Figs. S2-S3, Tables S1-S2)
# Pure Python/SciPy; outputs are written to results/.
# -----------------------------------------------------------------------------
hopf_scan: | dirs
	python3 $(SCRIPTS_DIR)/run_hopf_scan.py

hopf_timestep: | dirs
	python3 $(SCRIPTS_DIR)/run_time_step_stability.py

hopf_plots: | dirs
	python3 $(SCRIPTS_DIR)/plot_hopf_results.py

hopf: hopf_scan hopf_timestep hopf_plots

figures: figure1 figure2 figure3 figures4_5 figure6

# -----------------------------------------------------------------------------
# Supplementary material
# -----------------------------------------------------------------------------
supplement_defense: $(SUPP_DEFENSE_BIN) | dirs
	bash $(SCRIPTS_DIR)/run_supplement_robustness_checks.sh
	python3 $(SCRIPTS_DIR)/plot_supplement_robustness_checks.py $(SUPP_OUT) $(FIG_SUPP_DIR) || true

# Optional legacy/supplemental cycle calculation. This is not part of the main
# manuscript workflow unless explicitly called.
supplement_cycle: $(SUPP_CYCLE_BIN)

supplement: supplement_defense

# -----------------------------------------------------------------------------
# Checks and cleaning
# -----------------------------------------------------------------------------
check-python:
	python3 -m py_compile $(SCRIPTS_DIR)/*.py

check-names:
	@echo "Checking for old executable/source/script references..."
	@grep -R -nE "phase_sweep_lu|figure2_lu|figure4_lu|figure5_cover_lu|figure6_lu|defense_batch_lu|figure7_lu" $(SCRIPTS_DIR) README.md docs 2>/dev/null && { echo "Old C++ file-name references remain. Please update them."; exit 1; } || true
	@grep -R -nE "plot_figure2_trajectories|plot_figure3_profiles|plot_figure4_surface_tension|plot_figure5_cover_protocol|plot_figure6_growth_rate|plot_defense_batch" $(SCRIPTS_DIR) README.md docs 2>/dev/null && { echo "Old Python script-name references remain. Please update them."; exit 1; } || true
	@echo "No old-name references found."

clean:
	rm -rf $(BIN_DIR)

clean-data:
	rm -rf $(DATA_DIR) $(LOG_DIR)

distclean: clean clean-data
	rm -rf $(FIG_MAIN_DIR) $(FIG_SUPP_DIR)
