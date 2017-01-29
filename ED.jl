using ASCIIPlots

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

      Xb[2*Npara+(1:Npara),1]=rotmat[1,1]*Xb[Npara+(1:Npara),1] +rotmat[1,2]*Xb[Npara+(1:Npara),2];
      Xw[2*Npara+(1:Npara),1]=rotmat[1,1]*Xw[Npara+(1:Npara),1] +rotmat[1,2]*Xw[Npara+(1:Npara),2];
      Xb[2*Npara+(1:Npara),2]=rotmat[2,1]*Xb[Npara+(1:Npara),1] +rotmat[2,2]*Xb[Npara+(1:Npara),2];
      Xw[2*Npara+(1:Npara),2]=rotmat[2,1]*Xw[Npara+(1:Npara),1] +rotmat[2,2]*Xw[Npara+(1:Npara),2];

      vb=sortperm(Xb[:,1]+Npara*Xb[:,2])
      vw=sortperm(Xw[:,1]+Npara*Xw[:,2])
      Xb=Xb[vb,:];
      Xw=Xw[vw,:];

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

  return Bondsb, Bondsw, boundaryb, boundaryw
end


nx=2

Xb, Xw= MakeHexagon(nx)
Bondsb, Bondsw, boundaryb, boundaryw=MakeBonds(nx,Xb,Xw);

n=3*nx^2
nstates=2^n;
psi=collect(0:(nstates-1));

println("about to make M")
M=spzeros(Int8,nstates^2,nstates^2);
println("Made M")
for ii in 1:nstates
    for jj in 1:nstates

        for kk in 1:n

            if ((psi[ii])&2^(kk-1)) * ((psi[jj])&2^(Bondsb[kk,3]-1))!=0
                M[(ii-1)*nstates+jj,(ii-1)*nstates+jj]+=1
            else
                M[(ii-1)*nstates+jj,(ii-1)*nstates+jj]-=1
            end

        end
    end
    println(ii)
end
println("finished z")
for kk in 0:2
    bb=Bondsb[kk+1,1]-1
    for ii in 1:nstates
        if ( psi[ii]&(2^kk ) )==0
            neiw=psi[ii]+2^kk
        else
            neiw=psi[ii]-2^kk
        end

        for jj in 1:nstates
            if ( psi[jj]&(2^bb) ) ==0
                neib=psi[jj]+2^bb
            else
                neib=psi[jj]-2^bb
            end
            M[(ii-1)*nstates+jj,neiw*nstates+neib+1]+=1
        end
    end
end
println("finished x")
zed=true
for kk in 0:2
    bb=Bondsb[kk+1,2]-1
    for ii in 1:nstates
        if ( psi[ii]&(2^kk) )==0
            neiw=psi[ii]+2^kk
            zed=true
        else
            neiw=psi[ii]-2^kk
            zed=false
        end

        for jj in 1:nstates

            if ( psi[jj]&(2^bb) ) ==0
                neib=psi[jj]+2^bb
                if zed==true
                    M[(ii-1)*nstates+jj,neiw*nstates+neib+1]+=1
                else
                    M[(ii-1)*nstates+jj,neiw*nstates+neib+1]+=-1
                end
            else
                neib=psi[jj]-2^bb
                if zed==false
                    M[(ii-1)*nstates+jj,neiw*nstates+neib+1]+=1
                else
                    M[(ii-1)*nstates+jj,neiw*nstates+neib+1]+=-1
                end
            end

        end
    end
end
println("finished y")
F=eigs(M,nev=100,which="SM")

print(lineplot(collect(1:nstates^2),F[:values]))
