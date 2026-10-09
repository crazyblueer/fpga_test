# fpga_test
# Hybrid Beamformer Combiner on DE10-Standard (Cyclone V SoC)

Digital prototype of a hybrid (analog + digital) beamformer for a 64-element, 8-RF-chain array at 28 GHz.
The beam points at the user and places nulls on jammers.

**Status:** the combiner (`hybrid_bf.v`) is verified in simulation against a bit-exact Python model and compiles in Quartus.
It is **not yet connected to the ARM (HPS)** and **not yet run on the board**.

---

## 1. What this project does

The Python script computes hybrid beamforming weights `(F_RF, w_BB)` from the received data.
The system is split into two parts:

| Part | Runs | Where |
|---|---|---|
| **Combiner**: `y = w_hybrid^H x` for every snapshot | every sample, fast | **FPGA** (`hybrid_bf.v`) |
| **Weight computation**: R64 reconstruction, pilot angle, MVDR loop, phase quantization | once per update, slow | **ARM (HPS), in C** (to do) |

The ARM computes the weights and writes them into the FPGA registers. The FPGA applies them to every snapshot.

### Combiner structure

```
x[0..63] (complex, Q1.15)
   |  8 groups of 8 elements, phase-only coefficients c[k]  (analog stage, F_RF)
   v
z[0..7]  (8 complex values, Q4.15)
   |  8 complex weights b[r]                               (digital stage, w_BB)
   v
y        (1 complex value, 24-bit real + 24-bit imag)
```

- Architecture: **one serial complex MAC**, 64 cycles (analog stage) + 8 cycles (digital stage) = **72 cycles per snapshot**
- At 50 MHz: about 1.44 us per snapshot, about 0.69 M snapshots/s

---

## 2. Repository layout

```
repo/
  python/     beamformer script (float model + fixed-point model + golden-vector export)
  verilog/    hybrid_bf.v (design), tb_hybrid_bf.v (testbench, simulation only)
  vectors/    x_vec.hex, coef.hex, bb.hex, y_gold.hex
  quartus/    fpga_test.qpf, fpga_test.qsf, fpga_test.sdc
  c/          (to do) C port of the weight algorithm
  README.md
```

Do **not** commit Quartus generated folders (`db/`, `incremental_db/`, `output_files/`). Add them to `.gitignore`.

---

## 3. Tools and versions

| Tool | Version | Used for |
|---|---|---|
| Python 3 + numpy + matplotlib | (fill in) | Reference model, golden vectors |
| Icarus Verilog | v12 (20220611) | Simulation (no license needed) |
| Quartus Prime Lite | 24.1std (Build 1077) | Synthesis, fitting, timing |
| Target device | 5CSXFC6D6F31C6 (DE10-Standard) | |

Questa/ModelSim is **not** used: it needs a license file that is not set up.

---

## 4. How to verify (do this after every change)

```
cd verilog
iverilog -g2012 -o sim hybrid_bf.v tb_hybrid_bf.v
vvp sim
```

The `.hex` files must be in the directory where you run `vvp` (the testbench reads them with `$readmemh`).

**Expected output:**

```
DONE: 200 snapshots, 0 errors
```

**Rule: no change is accepted unless this line is printed.** The comparison is bit-exact against `y_gold.hex`.

---

## 5. Data formats (the contract between Python, Verilog and C)

All values are signed two's complement. Q1.15 means integer = value x 32768.

| Item | Format |
|---|---|
| Input sample `x[n]` | 16-bit re + 16-bit im, Q1.15. Python scales raw data by **1/16** to avoid clipping |
| Coefficient `c[n]` (analog stage) | 16-bit re + 16-bit im, Q1.15. Value is `exp(-j*2*pi*idx/64)`, **already conjugated**, `idx` = 6-bit phase index. Hardware computes `sum(c*x)` with no conjugation |
| Digital weight `b[r]` | 16-bit re + 16-bit im, Q1.15. Value is `conj(w_BB[r]) / sqrt(D)`, then scaled by `0.99 / max|b|`. The 1/sqrt(D) amplitude of F_RF is folded in here. Overall scale does not change the beam shape |
| Intermediate `z[r]` | 20-bit, Q4.15 (sum of 8 products, `>>15`) |
| Output `y` | 24-bit re + 24-bit im, accumulated then `>>15` |

### File layouts (plain text, one hex word per line)

| File | Lines | Word layout |
|---|---|---|
| `x_vec.hex` | 200 x 64 | `{re[15:0], im[15:0]}` (8 hex digits). Snapshot `t`, element `n` is line `t*64+n` |
| `coef.hex` | 64 | `{cr, ci}` (8 hex digits). Element order 0..63; chain `r` owns elements `8r..8r+7` |
| `bb.hex` | 8 | `{br, bi}` (8 hex digits) |
| `y_gold.hex` | 200 | `{yr[23:0], yi[23:0]}` (12 hex digits) |

These are **not** Quartus memory-initialization files. Do not add them to the Quartus project (if Quartus asks for a "word size", cancel).

---

## 6. Module interface (`hybrid_bf`)

| Port | Dir | Description |
|---|---|---|
| `clk`, `rst` | in | Clock, synchronous reset |
| `start` | in | One-cycle pulse: begin a snapshot |
| `x_we`, `x_addr[5:0]`, `x_re`, `x_im` | in | Write one input sample (64 entries) |
| `c_we`, `c_addr[5:0]`, `c_re`, `c_im` | in | Write one analog-stage coefficient (64 entries) |
| `b_we`, `b_addr[2:0]`, `b_re`, `b_im` | in | Write one digital weight (8 entries) |
| `y_re`, `y_im` | out | Result (24-bit signed) |
| `done` | out | High when `y` is valid |

FSM: `IDLE -> S1 (64 cycles) -> S2 (8 cycles) -> FIN -> IDLE`.

**Do not change these ports without agreeing on it first.** The HPS bridge wrapper and the C code depend on them.

---

## 7. Results so far

### Algorithm (Python, scenario: user 10 deg, jammers 4 deg and 16 deg, SNR 0 dB, JSR 20 dB each, 64 elements, 3000 snapshots)

| Beamformer | SINR | Null at 4 deg | Null at 16 deg |
|---|---|---|---|
| Plain steering vector (no nulling) | -3.43 dB | none | none |
| Hybrid, floating point | -0.36 dB | -48.6 dB | -48.5 dB |
| Hybrid, quantized (6-bit phase, Q1.15) | not measured | -50.8 dB | -48.0 dB |

- Pilot angle estimate: 9.9 deg (true 10 deg)
- Upper bound with this normalization is 0 dB SINR, so the hybrid result is within about 0.36 dB of it
- Quantization (6-bit phase, 16-bit weights) does not degrade the nulls

### Hardware (Quartus 24.1std, baseline compile with 50 MHz `.sdc`)

| Metric | Value |
|---|---|
| ALMs (needed) | about 2,350 / 41,910 (6%) |
| Registers | 4,810 |
| DSP blocks | 6 / 112 (5%) |
| Block memory bits | **0** |
| Pins | 166 (not assigned, not used on the board yet) |
| **Fmax (Slow 1100mV 85C)** | **61.72 MHz** |
| **Setup slack at 50 MHz** | **+3.798 ns** (50 MHz met) |
| Hold slack | (fill in) |
| Compile time | about 7.7 min (Fitter 7 min) |

Timing constraint file `fpga_test.sdc`:

```
create_clock -name clk -period 20.000 [get_ports clk]
derive_clock_uncertainty
set_false_path -from [get_ports {x_* c_* b_* start rst}]
set_false_path -to   [get_ports {y_* done}]
```

---

## 8. Known weaknesses (baseline to improve)

1. **Memories are registers, not M10K.** Quartus reports "RAM logic xr, xi, cr, ci, br, bi, zr, zi is uninferred due to asynchronous read logic". The arrays are read combinationally, so they became 4,810 registers plus large multiplexers.
2. **Low Fmax (61.7 MHz).** One clock cycle does array read, multiply and accumulate. Pipelining would raise it well above 100 MHz.
3. **Slow compile** (7 min) caused mostly by the register and multiplexer tangle.
4. **Stale pin assignments** (`LEDR`, `SW`) in the `.qsf` from the old LED counter project. They only produce warnings.
5. **Not connected to the ARM** and **no pin mapping** (the final design uses internal bridge signals, so no pins are needed for these ports).

---

## 9. Tasks

### Task A: FPGA optimization (owner: FPGA teammate)

Goal: improve area and Fmax **without changing results**.

1. Move `x`, `c`, `b`, `z` into M10K using synchronous reads (one-cycle read latency, extra FSM state to wait for data).
2. Pipeline the datapath: register after memory read, after the multiplier, after the accumulator.
3. Optionally add parallel multipliers if throughput matters.
4. After **every** change: rerun the Icarus testbench (must print `0 errors`), recompile in Quartus, and record ALMs, registers, M10K bits, DSPs, Fmax and slack against the baseline in section 7.

Constraints: keep the golden vectors unchanged. Keep module ports unchanged unless agreed. The testbench's wait for `done` must still work if latency changes.

### Task B: C port of the weight algorithm (owner: C teammate)

Goal: reproduce what the Python script does, in C (`double`, C99, `<complex.h>`), producing `coef` (64 phase indices) and `b` (8 complex Q1.15 values).

Pieces:
1. R64 reconstruction (sweep, cross-correlations `pv`, least-squares solve using the precomputed `pinv(Bm^H Bm) Bm^H`, which depends only on geometry)
2. Pilot-based user angle estimate
3. Initial F_RF phases, then the 15-iteration loop: `R_eff` (8x8), `w_BB` = MVDR (needs an 8x8 complex inverse), F_RF update from phases
4. Quantization: 6-bit phase index, Q1.15 for `b` with the same scaling as section 5
5. Output in the same format as `coef.hex` and `bb.hex`

Optional first experiment (in Python): replace the 64x64 LCMV initialization with the phases of the user's steering vector. If the nulls stay deep, the 64x64 inverse can be dropped from C entirely.

Done means: C output matches Python's `coef.hex` and `bb.hex` (indices equal or within one step on few elements, `b` within a few LSBs), **and** the resulting null depths are still about -48 dB.

Structure the code as functions (`reconstruct_R64`, `estimate_angle`, `design_weights`, `quantize`) plus a `main`, so it can later be called by the program that talks to the FPGA. Develop on a PC first; cross-compile for the ARM later. Input data for the C program is exported from Python into `vectors/` (to do).

### Task C: Board bring-up and bridge (owner: board holder)

1. Practice the full flow on the LED counter: pin assignments, `.sdc`, compile, program via USB-Blaster.
2. Boot Linux on the ARM from the Terasic SD card image, with serial or SSH access.
3. Platform Designer: HPS plus an Avalon wrapper around `hybrid_bf` on the lightweight HPS-to-FPGA bridge.
4. C test program: `mmap` the bridge, write `c`, `b`, `x`, pulse `start`, read `y`, compare with `y_gold.hex`.
5. Integrate with Task B (C weights written into the FPGA) and Task A (optimized combiner).

Register map: **to be defined** before Task C starts (one address each for `c[k]`, `b[r]`, `x[k]`, `start`/`done`, `y_re`, `y_im`).

---

## 10. Team rules

- One owner per file. Use Git branches and merge deliberately.
- Nothing is merged unless its check passes:
  - Verilog: Icarus prints `0 errors`
  - C: matches Python `coef.hex` and `bb.hex`, null depths still deep
- Any change to data formats, bit widths or module ports must be agreed first and updated in this README.
- Weekly short sync: each person shows one result (for example "C matches Python", "Fmax is 90 MHz").

---

## 11. Scope note

The DE10-Standard has no 28 GHz front end and no 64 real antennas. This is a **digital prototype** driven by generated or recorded test data. The analog phase shifters (F_RF) are modeled as complex multiplies by `exp(j*phi)`.
Confirm the intended scope with the advisor.