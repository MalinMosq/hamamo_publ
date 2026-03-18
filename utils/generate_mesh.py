# Malin Mosquera
import gmsh
import sys
import numpy as np

# ------------- Read from files -------------------------

# Read all the points to numpy array (matrix)
all_points = np.loadtxt('data_domain/points.txt')

# Read all the edges to numpy array (matrix)
all_edges = np.loadtxt('data_domain/edges.txt','int16')

# Read all the subdomains EDGE indices to list containing numpy arrays (vectors)
f = open('data_domain/subdomains_e.txt','r')
cycles_e = []
for count, line in enumerate(f):
	cycle = []
	for y in line.split():
		if int(y) > 0:
			cycle.append(int(y))
	cycles_e.append(cycle) #cycles_e.append(np.array([int(y) for y in line.split()]))
f.close()

# Read all the subdomains POINT indices to list containing numpy arrays (vectors)
f = open('data_domain/subdomains_p.txt','r')
cycles_p = []
for count, line in enumerate(f):
	cycle = []
	for y in line.split():
		if int(y) > 0:
			cycle.append(int(y))
	cycles_p.append(cycle) #np.array([int(y) for y in line.split()])
f.close()

# -------------- Gmesh construction ---------------------

gmsh.initialize()
gmsh.model.add("domain")
lc = 1

# Points
p_gmsh = []
for p in all_points:
	p_gmsh.append(gmsh.model.geo.add_point(p[0],p[1],0,lc))

# Lines
e_gmsh = []
for e in all_edges:
	e_gmsh.append(gmsh.model.geo.add_line(e[0], e[1]))

# Curve loops (boundaries of each subdomain) and plane surfaces (subdomains)
curve_loop = np.empty((len(cycles_e),1))
plane_surfaces = np.empty((len(cycles_e),1))
for current_cycle_index in np.arange(len(cycles_e)):
	current_cycle = cycles_e[current_cycle_index]
	connecting_value = cycles_p[current_cycle_index][0]
	string = "["
	for line in current_cycle:
		if all_edges[line-1,0] == connecting_value:
			string = string+"e_gmsh[% s]," %(line-1)
			connecting_value = all_edges[line-1,1]
		else:
			string = string+"-"+"e_gmsh[% s]," %(line-1)
			connecting_value = all_edges[line-1,0]
	string = string[:-1]+"]" # Remove last comma and add end paranthesis
	# Curve loops (boundaries of each subdomain)
	curve_loop[current_cycle_index] = gmsh.model.geo.add_curve_loop(eval(string))
	# Plane surfaces (subdomains)
	plane_surfaces[current_cycle_index] = gmsh.model.geo.add_plane_surface(curve_loop[current_cycle_index])

# Synchronise everything
gmsh.model.geo.synchronize()

# Meshes
gmsh.model.mesh.generate(1)
gmsh.model.mesh.generate(2)

# Write to file
# gmsh.write("data_mesh/mesh.msh")
gmsh.write("data_mesh/mesh_h0.m")








gmsh.finalize()


















