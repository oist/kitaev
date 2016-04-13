loc='~/Prog/TQM/kitaev/diagonalization.h5';
info=h5info(loc);
g=info.Groups(1).Name;

E=h5read(loc,strcat(g,'/E'));
U=h5read(loc,strcat(g,'/U'));
V=h5read(loc,strcat(g,'/V'));