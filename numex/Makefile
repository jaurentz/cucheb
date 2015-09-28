include ../make.inc

all: objects 

objects: FORCE
	@$(MAKE) -C ./parsec
	@$(MAKE) -C ./laplacian
	@$(MAKE) -C ./groundstates
	@$(MAKE) -C ./dimacs
	@$(MAKE) -C ./tables

FORCE:

clean:
	@$(MAKE) clean -C ./parsec
	@$(MAKE) clean -C ./laplacian
	@$(MAKE) clean -C ./groundstates
	@$(MAKE) clean -C ./dimacs
	@$(MAKE) clean -C ./tables
	
	
	
	
