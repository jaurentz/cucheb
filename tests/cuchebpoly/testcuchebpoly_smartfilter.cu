#include <cucheb.h>

/* driver */
int main(){

  // cuhebpoly 1
  cuchebpoly ccp1;

  // initialize to identity and print
  cuchebpoly_init(&ccp1);
  cuchebpoly_print(&ccp1);

  // step_filter
  cuchebpoly_smartfilter(-1.0,1.0,-0.01,0.01,&ccp1);
  cuchebpoly_print(&ccp1);

  // free memory
  cuchebpoly_destroy(&ccp1);

  // return 
  return 0;

}
