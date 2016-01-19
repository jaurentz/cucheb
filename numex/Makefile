include ../make.inc

all: objects 

objects: FORCE
	@$(MAKE) -C ./parsec
	@$(MAKE) -C ./graph
	@$(MAKE) -C ./engineering
	@$(MAKE) -C ./tables

FORCE:

clean:
	@$(MAKE) clean -C ./parsec
	@$(MAKE) clean -C ./graph
	@$(MAKE) clean -C ./engineering
	@$(MAKE) clean -C ./tables
	
	
	
	
