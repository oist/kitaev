n=151;
d=6/n;

omega=linspace(0,6,n);
height=zeros(n,1);

for i = 1:length(E);
   j=ceil((2*E(i))*25);
   height(j)=height(j)+w3(i);
end