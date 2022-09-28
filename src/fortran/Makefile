FC=f2py
FC_ARGS=-c

LIB=axis_hurr

SRC=hurr-vars.f90 \
	hurr-initial.f90 \
	hurr-fwd.f90 \
	hurr-tlm.f90 \
	hurr-adj.f90
${LIB}.so:	${SRC}
		${FC} ${FC_ARGS} -m ${LIB} ${SRC}
clean:
		rm ${LIB}.*.so
