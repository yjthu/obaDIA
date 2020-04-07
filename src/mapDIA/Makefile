all: mapDIA

objs = Pre main data Est Option Post Module print
objs := $(addsuffix .o, $(objs))
objs := $(addprefix src/, $(objs))
CXXFLAGS = -pedantic-errors -Wall -O3 -Inlopt/include

$(objs): nlopt

mapDIA: nlopt $(objs)
	$(CXX) -o $@ $(objs) -Lnlopt/lib -lnlopt

nlopt:
	cd nlopt-2.4.2/;\
	    ./configure --prefix=$(CURDIR)/nlopt;\
	    $(MAKE) install

.PHONY : clean
clean :
	$(RM) mapDIA $(objs)
