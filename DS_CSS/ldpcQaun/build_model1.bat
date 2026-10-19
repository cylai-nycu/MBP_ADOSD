@echo off
setlocal
cd /d "%~dp0"
gcc -std=c99 -O2 -fopenmp -w DS_CSS_model1.c bp_dec\bp_dec.c bp_dec\bp_llr.c lib_rand\lib_rand.c lib_rand\xoshiro256starstar.c lib_rand\splitmix64.c lib_math\fast_math.c OSD\OSD.c -lm -o DS_CSS_model1.exe
if errorlevel 1 exit /b %errorlevel%
echo Built DS_CSS_model1.exe
