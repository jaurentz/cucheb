#include <cucheb.h>

cuchebStatus_t cuchebSetGridBlocks(int n, dim3 *blockSize, dim3 *gridSize){

	// query device
	int dev;
	cudaDeviceProp prop;
	cuchebCheckError(cudaGetDevice(&dev),__FILE__,__LINE__);
	cuchebCheckError(cudaGetDeviceProperties(&prop,dev),__FILE__,__LINE__);

	// set blockSize
	*blockSize = dim3(prop.maxThreadsPerBlock,1,1);
	int num = blockSize->x*blockSize->y*blockSize->z;
	
	// set gridSize
	int powtwo;
	gridSize->x = (int)ceil((double)n/num);
	powtwo = (int)max(1,(int)floor(log2((double)prop.maxGridSize[0])));
	if(gridSize->x > pow(2,powtwo)){
		gridSize->y = (int)ceil((double)gridSize->x/pow(2,powtwo));
		gridSize->x = pow(2,powtwo);
	}
	powtwo = (int)max(1,(int)floor(log2((double)prop.maxGridSize[1])));
	if(gridSize->y > pow(2,powtwo)){
		gridSize->z = (int)ceil((double)gridSize->y/pow(2,powtwo));
		gridSize->y = pow(2,powtwo);
	}
	powtwo = (int)max(1,(int)floor(log2((double)prop.maxGridSize[2])));
	if(gridSize->z > pow(2,powtwo)){
		fprintf(stderr,"\nIn %s line: %d, GPU cannot support this many threads.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// return
	return CUCHEB_STATUS_SUCCESS;
}
