
using HDF5;

function MakeHexagon(Nx::Int)
      Npara=Nx^2;
      N=Npara*3;

      unit=transpose([cos(pi/6),sin(pi/6)])
      a=transpose([2*cos(π/6),0]);
      b=[cos(π/6),1+sin(π/6)];
      bM=transpose(repeat(b,outer=[1,Nx]))
      rotmat=[[cos(2*π/3) -sin(2*π/3)]
            [sin(2*π/3) cos(2*π/3)]]


      Xb=Array{Float64}(N,2);
      Xw=Array{Float64}(N,2);

      ## Creating the positions in the first parallelogram
      for i in 1:Nx
        Xb[i,:]=(i-1)*a;
        Xw[i,:]=unit+(i-1)*a;
      end

      for j in 2:Nx
        Xb[(Nx*(j-1)+1):(Nx*(j-1)+Nx),:]=Xb[1:Nx,:]+(j-1)*bM;
        Xw[(Nx*(j-1)+1):(Nx*(j-1)+Nx),:]=Xw[1:Nx,:]+(j-1)*bM;
      end

      subt=Xb[Npara,2]+1
      Xb[1:Npara,1]=Xb[1:Npara,1]-(Nx-1)*b[1];
      Xw[1:Npara,1]=Xw[1:Npara,1]-(Nx-1)*b[1];
      Xb[1:Npara,2]=Xb[1:Npara,2]-subt;
      Xw[1:Npara,2]=Xw[1:Npara,2]-subt;

      Xb[Npara+(1:Npara),1]=rotmat[1,1]*Xb[1:Npara,1]+rotmat[1,2]*Xb[1:Npara,2];
      Xw[Npara+(1:Npara),1]=rotmat[1,1]*Xw[1:Npara,1]+rotmat[1,2]*Xw[1:Npara,2];
      Xb[Npara+(1:Npara),2]=rotmat[2,1]*Xb[1:Npara,1]+rotmat[2,2]*Xb[1:Npara,2];
      Xw[Npara+(1:Npara),2]=rotmat[2,1]*Xw[1:Npara,1]+rotmat[2,2]*Xw[1:Npara,2];

      Xb[2*Npara+(1:Npara),1]=rotmat[1,1]*Xb[Npara+(1:Npara),1]+rotmat[1,2]*Xb[Npara+(1:Npara),2];
      Xw[2*Npara+(1:Npara),1]=rotmat[1,1]*Xw[Npara+(1:Npara),1]+rotmat[1,2]*Xw[Npara+(1:Npara),2];
      Xb[2*Npara+(1:Npara),2]=rotmat[2,1]*Xb[Npara+(1:Npara),1]+rotmat[2,2]*Xb[Npara+(1:Npara),2];
      Xw[2*Npara+(1:Npara),2]=rotmat[2,1]*Xw[Npara+(1:Npara),1]+rotmat[2,2]*Xw[Npara+(1:Npara),2];

      vb=sortperm(Xb[:,1]+Npara*Xb[:,2])
      vw=sortperm(Xw[:,1]+Npara*Xw[:,2])
      Xb=Xb[vb,:];
      Xw=Xw[vw,:];

      h["Xb"]=Xb;
      h["Xw"]=Xw;

  return Xb, Xw
end

function MakeBonds(Nx::Int,Xb::Array{Float64,2},Xw::Array{Float64,2})
  Npara=Nx^2;
  N=Npara*3;

  nrowb=[collect(Nx:(2*Nx));collect((2*Nx-1):-1:(Nx+1))]
  nroww=flipdim(nrowb,1)

  i0b=ones(nrowb)
  i0w=ones(nroww)
  for i in 2:length(nrowb)
      i0b[i]=i0b[i-1]+nrowb[i-1]
      i0w[i]=i0w[i-1]+nroww[i-1]
  end
  push!(i0b,N+1);
  push!(i0w,N+1);

  offsetb=zeros(N+1)
  offsetw=zeros(N+1)
  for i in 1:Nx
      offsetb[i0b[i]:i0b[i+1]]=i-1
      offsetb[i0b[i+Nx]:i0b[i+1+Nx]]=Nx-i
      offsetw[i0w[i]:i0w[i+1]]=i-1
      offsetw[i0w[i+Nx]:i0w[i+1+Nx]]=Nx-i
  end

  Bondsw=Array{Int64}(N,3);
  Bondsb=Array{Int64}(N,3);
  for i in 1:N
      Bondsb[i,2]=i+offsetb[i];
      Bondsw[i,2]=i-offsetw[i];
      Bondsb[i,1]=i+1+offsetb[i];
      Bondsw[i,1]=i-1-offsetw[i];
      Bondsb[i,3]=i-Nx;
      Bondsw[i,3]=i+Nx;
  end

  boundaryb=ones(Int64,N,3);
  boundaryw=ones(Int64,N,3);
  i=1
  while isapprox(Xb[1,2],Xb[i,2])
      Bondsb[i,3]=N-Nx+i
      boundaryb[i,3]=0
      Bondsw[N-Nx+i,3]=i
      boundaryw[N-Nx+i,3]=0
      #println(N-Nx+i,' ',i)
      i+=1
  end

  Bondsw[i0w[1:Nx],1]=i0b[Nx+1+(1:Nx)]-1;
  boundaryw[i0w[1:Nx],1]=0
  Bondsb[i0b[Nx+1+(1:Nx)]-1,1]=i0w[1:Nx];
  boundaryb[i0b[Nx+1+(1:Nx)]-1,1]=0

  Bondsw[i0w[1+(1:Nx)]-1,2]=i0b[Nx+(1:Nx)];
  boundaryw[i0w[1+(1:Nx)]-1,2]=0
  Bondsb[i0b[Nx+(1:Nx)],2]=i0w[1+(1:Nx)]-1;
  boundaryb[i0b[Nx+(1:Nx)],2]=0

  h["Bondsw"]=Bondsw;
  h["Bondsb"]=Bondsb;
  h["Bondaryw"]=boundaryw;
  h["Bondaryb"]=boundaryb;
  h["i0w"]=i0w;
  h["i0b"]=i0b;

  return Bondsb, Bondsw, boundaryb, boundaryw
end

function MakeS(N::Int, J::Array,Bondsw::Array{Int}, bw::Array{Int})

      S=spzeros(N, N);

      h["J"]=J;

      for ii in 1:N
          for jj in 1:3
              S[ii,Bondsw[ii,jj]]=J[ii,jj]#*bw[ii,jj];
          end
      end
      h["S"]=full(S);
      return S;
end

function DiagS(S::SparseMatrixCSC{Float64,Int64})

  @time F=svdfact(full(S))

  h["E"]=F[:S];
  h["U"]=F[:U];
  h["V"]=F[:V];
  h["Etot"]=sum(F[:S]);

  return F;
end

function LDOS(F::Base.LinAlg.SVD)
  site=round(Int,N/2);
  n=150;
  omega=collect(linspace(0,6,150));
  heights=zeros(omega);

  for  j in 1:N;
    k=ceil(Int,F[:S][j]*50);
    heights[k]+=.5*F[:U][j,site].^2;
  end

  h["omega"]=omega;
  h["heights"]=heights;

  return 1;
end

function Fluxes(i0::Int,nflip::Int)
  J=Array{Int64}(N,3);
  for l in 1:N
    for m in 1:3
        J[l,m]=1;
    end
  end
  for l in 0:(nflip-1)
    J[i0+l,3]=-1
  end
  return J
end


Nx=50;
N=3*Nx^2;
i=round(Int,N/2);
nflip=10;

fid=h5open("twoflux_periodic.h5","r+")
if exists(fid,"$Nx")
  g=fid["$Nx"];
else
  g=g_create(fid,"$Nx");
end

version=readall(`git rev-list --count HEAD`)[1:end-1];
if exists(g,version)
  k=g[version];
else
  k=g_create(g,version);
end

flipper="nflip $nflip"
if exists(k,flipper)
  error("You already ran this combo...")
else
  h=g_create(k,flipper)
end

Xb, Xw = MakeHexagon(Nx);
Bondsb, Bondsw, bb, bw= MakeBonds(Nx,Xb,Xw);
J=Fluxes(i-nflip,nflip);
S=MakeS(N,J,Bondsw, bw)

println("starting ED")
F=DiagS(S)
println("finished ED")

ret=LDOS(F);

close(fid);
