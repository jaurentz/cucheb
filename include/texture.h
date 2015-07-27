#pragma once

texture<int2, 1> tex_double;

void cuda_bind_tex(size_t *off, double *x, size_t size) {
  CUDA_SAFE_CALL(
  cudaBindTexture(off, tex_double, x, size) );
}

void cuda_unbind_tex() {
  CUDA_SAFE_CALL(
  cudaUnbindTexture(tex_double) );
}

