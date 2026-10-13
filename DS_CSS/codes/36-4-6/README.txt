[[36,4,6]] BB code with sparse original checks + redundant syndrome measurements.
Original BB: x^3=y^6=1, A=x+y^2+y^3, B=1+y+x^2.
X checks 23 (18 original weight-6 + 5 redundant), Z checks 23 (18 original weight-6 + 5 redundant).
46 measurement-error virtual qubits, 36 data qubits, 82 total columns, 46 DS rows.
Physical check ranks 16+16; physical code [[36,4,6]] (distance from known BB construction).
Pauli encoding: 0=I, 1=X, 2=Z, 3=Y.
DS_A4.txt: GF4 alist, 82 columns x 46 rows. Data columns 0..35; virtual 36..58 correspond to X check measurement faults (Z virtual); virtual 59..81 correspond to Z check measurement faults (X virtual).
DS_G4s.txt: GF4 alist, 82 columns x 40 rows, same physical normalizer as previous construction, zero virtual columns.
DS distance 5 under min_{nontrivial physical Pauli e} wt(e)+wt(H e), exactly checked by enumerating all pure X and Z errors of physical weights 1..4. Mixed Pauli faults cannot give a smaller value because each sector's contribution lower bounds total.
The redundant measurements must be physically implemented as additional stabilizer measurements; do not model them as independent measurements if they are merely computed by XOR from already noisy syndrome bits.
Dedicated C decoder dimensions: n_data=36, mX=23, mZ=23, N=82, M=46, k=4, M_G=40. Decoder code must support these dimensions.
