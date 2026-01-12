# =============================
# Compiler and flags
# =============================

CC      = gcc
CFLAGS  = -O2 -Wall
TARGET  = BPOSD

# =============================
# Source files
# =============================

SRC = \
	BPOSD.c \
	bp_dec/bp_dec.c \
	bp_dec/bp_llr.c \
	lib_rand/splitmix64.c \
	lib_rand/lib_rand.c \
	lib_rand/xoshiro256starstar.c \
	lib_math/fast_math.c \
	OSD/OSD.c

# =============================
# Build rules
# =============================

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
