Nx=4
Ny=3


using Iterators;
using ASCIIPlots;
#help define the atoms
unit=transpose([cos(π/6),sin(π/6)]);
a=transpose([2*cos(π/6),0]);
b=[cos(π/6),1+sin(π/6)];

#number of atoms
N=Nx*Ny
Nt=Nx*Ny*2

function armod(x::UInt,y::UInt)
    if x>y || x==0 || x<0
        return mod(x-1,y)+1
    else
        return x
    end
end

#Initializing all the arrays
bM=transpose(repeat(b,outer=[1,Nx]))
function pos()
  Xb=Array{Float64}(N,2);
  Xw=Array{Float64}(N,2);

    for i in 1:Nx
      Xb[i,:]=(i-1)*a;
      Xw[i,:]=unit+(i-1)*a;
    end

    for j in 2:Ny
      Xb[(Nx*(j-1)+1):(Nx*(j-1)+Nx),:]=Xb[1:Nx,:]+(j-1)*bM;
      Xw[(Nx*(j-1)+1):(Nx*(j-1)+Nx),:]=Xw[1:Nx,:]+(j-1)*bM;
    end
  return Xb, Xw
end

function bonds()
  bw=Array{UInt16}(N,3);
  #bb=Array{UInt16}(N,3);

  for i in 1:N
    #bb[i,1]=i;
    bw[i,1]=i;
    #bb[i,2]=armod(i-1,Nt);
    bw[i,2]=armod(i+1,N);
    #bb[i,3]=armod(i-Nx,Nt);
    bw[i,3]=armod(i+Nx,N);
  end
  return bw#, bb
end

bw=bonds()



J=ones(UInt8,N,3);
#for i in 1:N
#    for j in 1:3
#        J[i,j]=1;
#    end
#end
#cb=collect(product(repeated(0:1,N)...));
#cw=collect(product(repeated(0:1,N)...));

u=ones(Int8,N,3);
function Wp(i::Int)
    j=armod(i+1,N);
    k=armod(i-Nx+N,N);
    m=armod(i-Nx+1+N,N);
    w=u[i,1];
    w*=u[i,3];
    w*=u[j,2];
    w*=u[j,3];
    w*=u[k,2];
    w*=u[m,1];
    return w
end

function D()


  return psinew
end

#for i=1:N
#  println(Wp(i))
#end

A=zeros(Int8,N,N);
for i in 1:N
  for j in 1:3
      A[i,bw[i,j]]= J[i,j]*u[i,j];
  end
end
#println(imagesc(A))
