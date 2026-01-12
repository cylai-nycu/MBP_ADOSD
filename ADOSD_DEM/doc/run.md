\# BPOSD Project Folder Structure



\## 1. BPOSD\_CB\_proj



This is the original Code::Blocks project folder, containing the complete environment for compiling and running the program.



\### Folder Structure



BPOSD\_CB\_proj/

├── bin/ # Compiled executables or binary files

├── codes/ # Various source code files

├── obj/ # Object files generated during compilation (.o)

├── Results/ # Folder for program output and results

└ldpcQuan\_CB\_proj # Main program





ldpcQuan/

├── bp\_dec/ # BP decoder related code

├── lib\_math/ # Mathematical utility library

├── lib\_rand/ # Random number generator library

└── OSD/ # Ordered Statistics Decoding (OSD) module





Steps:

1  cd C:\\ADOSD\\ldpcQaun

2  gcc BPOSD.c bp\_dec/bp\_dec.c bp\_dec/bp\_llr.c lib\_rand/splitmix64.c lib\_rand/lib\_rand.c 	lib\_rand/xoshiro256starstar.c lib\_math/fast\_math.c OSD/OSD.c -o BPOSD
3  .\\BPOSD.exe










\### fast\_math.c / fast\_math.h



This file contains fast approximate functions for exp, log, and pow.

Based on Martin Ankerl's and Nic Schraudolph's algorithms, with some modifications.




