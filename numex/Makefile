include ../make.inc

all: objects 

objects: FORCE
	@$(MAKE) -C ./parsec
	@$(MAKE) -C ./parsec2
	@$(MAKE) -C ./groundstates
	@$(MAKE) -C ./dimacs
	@$(MAKE) -C ./tables

FORCE:

clean:
	@$(MAKE) clean -C ./parsec
	@$(MAKE) clean -C ./parsec2
	@$(MAKE) clean -C ./groundstates
	@$(MAKE) clean -C ./dimacs
	@$(MAKE) clean -C ./tables
	
	
	
	
