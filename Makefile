# Wannier program's Makefile has been used in building this Makefile
# Ref:: (http://www.wannier.org/download/)

ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../src/Makefile.ether

all: ether     

default: ether

objdir: 
	@( cd $(ROOTDIR) && if [ ! -d obj ] ; \
		then mkdir obj ; \
	fi ) ;

ether: objdir
	(cd $(ROOTDIR)/obj && $(MAKE) -f $(REALMAKEFILE) ether )

clean:
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f ether  
	cd $(ROOTDIR) && rm -f src/*~
	@( cd $(ROOTDIR) && if [ -d obj ] ; \
		then cd obj && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf obj ; \
	fi );


.PHONY: ether default all clean objdir
