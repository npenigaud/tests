module load nvhpc/22.11

PREFIX_FIAT=/home/gmap/mrpm/marguina/fiat/pginvtx+mplacc-build

mpif90 \
  -acc -Mcuda -Minfo=all mpc9.F90 -o mpc9.x \
  -L/opt/softs/nvidia/hpc_sdk/Linux_x86_64/22.11/math_libs/11.8/targets/x86_64-linux/lib \
  -I$PREFIX_FIAT/module/fiat -L$PREFIX_FIAT/lib -lfiat -Wl,-rpath,$PREFIX_FIAT/lib

