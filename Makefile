
MAKER         = pgcc
OPTFLAGS    = -acc -Minfo=acc
DEBUG       = 
ACCFLAGS    = 
LIB_LIST    = -lm 

SUFFIX = c

PROGRAM = moldyn

SRCS  = main.c    input-parameters.c read-input.c \
        initialise-particles.c pseudorand.c loop-initialise.c force.c movout.c \
        movea.c moveb.c sum-energies.c hloop.c scalet.c tidyup.c check_cells.c output_particles.c force_ij.c

### End User configurable options ###

FFLAGS =  $(INCLUDE_DIR) $(OPTFLAGS)
LIBS = $(LIB_PATH) $(LIB_LIST)
FLIBS = $(FLIB_PATH) $(LIB_LIST)

OBJS = $(SRCS:.$(SUFFIX)=.o)

${PROGRAM}: ${OBJS}
	$(MAKER) $(OPTFLAGS) -o $(PROGRAM) $(OBJS) $(FLIBS)

${OBJS}:
	${MAKER} ${DEBUG} ${OPTFLAGS} ${INCLUDE_DIR} -c $*.${SUFFIX}

clean:
	rm -f ${OBJS}

Clean:
	rm -f ${PROGRAM} ${OBJS}
