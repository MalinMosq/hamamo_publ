%% Matlab mesh
%% domain, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 7;
msh.POS = [
0 0.41780915 0;
1 0.55257784 0;
0 0 0;
1 0 0;
0 1 0;
1 1 0;
0.5 0.7522230824999999 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 2 0
 1 3 0
 1 5 0
 2 4 0
 2 6 0
 3 4 0
 5 6 0
];
msh.TRIANGLES =[
 2 4 1 0
 4 3 1 0
 5 1 7 0
 2 6 7 0
 1 2 7 0
 6 5 7 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
 5 0
 6 0
];
