include ../make.inc

all: objects 

objects: FORCE
	@$(MAKE) -C ./cuchebutils
	@$(MAKE) -C ./cuchebstats
	@$(MAKE) -C ./cuchebpoly
	@$(MAKE) -C ./cuchebmatrix
	@$(MAKE) -C ./cucheblanczos

FORCE:

clean:
	@$(MAKE) clean -C ./cuchebutils
	@$(MAKE) clean -C ./cuchebstats
	@$(MAKE) clean -C ./cuchebpoly
	@$(MAKE) clean -C ./cuchebmatrix
	@$(MAKE) clean -C ./cucheblanczos
	
	
	
	
