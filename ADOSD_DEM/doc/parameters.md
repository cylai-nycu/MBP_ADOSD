Simulation Configuration Explanation



Decoder and OSD settings:



\#define OSDW (-2) – Determines the type of OSD used:



-1: No OSD applied.



-2: Adaptive Partial OSD (APOSD).



Other values: specify the OSD layer directly.



\#define ADOSDw (2) – Sets the OSD complexity budget:



0: OSD-0 (minimal correction).



2: OSD-2 (higher complexity, roughly double the OSD-0 effort).



Reliability and iteration control:



\#define RelThr (0.99) – Reliability threshold used to classify bits as reliable or unreliable.



\#define MAX\_ITER (10) – Maximum number of iterations for belief propagation (BP) before invoking OSD.



Error and simulation parameters:



\#define P\_ERR 0.01 – Physical error rate for the circuit-level simulation.



\#define LE\_tar 1000 – Target number of logical error events to collect for statistics.



\#define Shot\_max 100000 – Maximum number of Monte Carlo samples per simulation point.



\#define ShowTime 1000 – Interval (in samples) at which progress and statistics are printed.



Code variant settings:



\#define rev 1 – Indicates whether additional X detectors are used in the rotated surface code:



0: With additional X detectors.



1: Without additional X detectors.



Surface code selections:



Rotated surface codes without additional X detectors (Z-only measurement):



\#define ORI 19451400 – Corresponds to the \[\[81,1,9]] code at the circuit level.



Other smaller codes (commented out) are available for \[\[9,1,3]], \[\[25,1,5]], and \[\[49,1,7]].



Rotated surface codes with additional X detectors (full circuit-level parity checks):



Codes from \[\[9,1,3]] up to \[\[81,1,9]] are provided but currently commented out.

