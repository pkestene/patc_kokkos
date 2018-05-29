see slides
http://on-demand.gputechconf.com/gtc/2016/presentation/s6510-jeff-larkin-targeting-gpus-openmp.pdf

See also source code:
https://github.com/jefflarkin/GTC16-S6510


I you want to try this code on ouessant, e.g. using the clang compiler:
module load cuda/9.0
module load llvm/clang/17.10-1
export OLCF_CUDA_ROOT=$CUDA_ROOT
