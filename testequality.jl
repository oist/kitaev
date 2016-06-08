using HDF5;

fid=h5open("diagonalization.h5","r")
h=fid["50/5"]

E=read(h["E"])
S=read(h["S"])
U=read(h["U"])
V=read(h["V"])

ze=S-U*diagm(E)*transpose(V);

println(maximum(ze))

close(fid)
