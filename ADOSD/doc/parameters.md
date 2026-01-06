Configuration via ldpc_parm.h

All code and decoder configurations are specified in the header file
ldpc_parm.h.
Users select one quantum code and one set of decoder parameters by
enabling the corresponding #define entries.

###1. Quantum Code Selection

The quantum code to be simulated is chosen by defining the macro ORI.
Only one ORI definition should be active at a time.

#define ORI 14412   // [[144,12,12]] BB code


Each ORI value corresponds to a pre-defined stabilizer/normalizer matrix
stored in the codebase.

Supported Code Families

Examples of available codes include:

Surface codes
Rotated surface codes (optionally XZZX-type)
Toric codes
Standard toric codes
Twisted / XZZX toric codes
Color codes
(4.8.8) color codes
Color codes on the torus
BB codes
Lifted-connected surface codes

Each entry is annotated using the standard notation [[n, k, d]].


Additional Geometry Parameters

Some code families require additional parameters, e.g.:
#define Surf_d  13   // surface code distance (odd only)
#define hex_d   5    // color code distance

These parameters must be consistent with the selected ORI entry.

###2. Decoder Configuration

The decoder is based on belief propagation (BP) and its variants, optionally
combined with ordered statistics decoding (OSD).

2.1 BP Variant
#define LLR_BP  0


0: Linear-domain BP
1: Log-likelihood-ratio (LLR) BP


2.2 OSD / ADOSD Configuration
#define OSDW   (-2)
#define ADOSDw (2)


OSDW = -1 : BP only (no OSD)
OSDW = -2 : AP-OSD / ADOSD (approximate probabilistic OSD)

ADOSDw specifies the effective OSD order (e.g., 2 ≈ OSD-2 complexity)


2.3 Reliability-Based Subset Reduction
#define RelThr  (0.999995)

Defines the reliability threshold used to identify highly reliable variable
nodes for subset reduction.


2.4 Iteration Control
#define MAX_ITER  (100)

Maximum number of BP iterations.


2.5 MBP / AMBP Step Size Parameters
#define AFP   160
#define AFP2  0
#define AFPINC  (0)


AFP (alpha prime): controls the damping / step size in MBP
(larger value → smaller update step)

Mode Selection:

MBP (fixed step size)

AFP = constant, AFP2 = 0, AFPINC = 0


AMBP (adaptive step size)

AFP    = starting value
AFP2   = ending value
AFPINC = increment per retry


Example:

AFP = 160, AFP2 = 30, AFPINC = -1



2.6 Message Update Schedule
#define RND_SCHE 0


0: Deterministic schedule
1: Randomized schedule at initialization
2: Randomized per check node
3: Randomized per iteration

Random schedules often improve convergence performance in AMBP.



2.7 Decoder Scheduling Mode
#define BY_DEC  20


Controls the decoding schedule and execution mode:

Value	Description
20	Parallel decoding (default)
24	Serial decoding