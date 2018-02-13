all:
	( cd src; $(MAKE) )
clean:
	( cd src; $(MAKE) clean )
veryclean:
	( cd src; $(MAKE) veryclean )
libraries:
	( cd src; $(MAKE) libraries )
