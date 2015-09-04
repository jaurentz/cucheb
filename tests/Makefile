include ../make.inc

all: tests

tests: FORCE
	@$(MAKE) -C ./cuchebpoly
	@$(MAKE) -C ./cuchebmatrix
	@$(MAKE) -C ./cucheblanczos
	@$(MAKE) -C ./cuchebblocklanczos

FORCE:

clean:
	@$(MAKE) clean -C ./cuchebpoly
	@$(MAKE) clean -C ./cuchebmatrix
	@$(MAKE) clean -C ./cucheblanczos
	@$(MAKE) clean -C ./cuchebblocklanczos
	
	
	
	
