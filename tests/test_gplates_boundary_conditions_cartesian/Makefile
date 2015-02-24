input: clean
	gmsh -2 -bin Mesh/spherical_shell.geo
	tar -xvjf gplates_data/GPlates_Velocities.tar.bz2

clean:
	rm -rf src/*.msh *.msh *.stat *.vtu *.node *.ele *.edge *.log-0 *.err-0 *checkpoint *.convergence \
	matrixdump matrixdump.info *.nc *.pyc *~*
