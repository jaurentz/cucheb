include ../make.inc

all: objects 

objects: FORCE
	@$(MAKE) -C ./parsec
	@$(MAKE) -C ./tables

FORCE:

clean:
	@$(MAKE) clean -C ./parsec
	@$(MAKE) clean -C ./tables
	
	
	
	
