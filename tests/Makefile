include ../make.inc

all: tests

tests: FORCE
	@$(MAKE) -C ./cuchebpoly
	@$(MAKE) -C ./cuchebmatrix
	@$(MAKE) -C ./cucheblanczos

FORCE:

clean:
	@$(MAKE) clean -C ./cuchebpoly
	@$(MAKE) clean -C ./cuchebmatrix
	@$(MAKE) clean -C ./cucheblanczos
	
	
	
	
