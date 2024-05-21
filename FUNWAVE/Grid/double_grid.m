function Znew = double_grid(Zold)

[n1 m1]=size(Zold);
m2=(m1-1)*2+1;
n2=(n1-1)*2+1;

for j=1:2:n2
for i=1:2:m2
Znew(j,i)=Zold((j+1)/2,(i+1)/2);
end
end
for j=1:2:n2
for i=2:2:m2-1
Znew(j,i)=0.5*(Znew(j,i-1)+Znew(j,i+1));
end
end
for j=2:2:n2-1
for i=1:1:m2
Znew(j,i)=0.5*(Znew(j-1,i)+Znew(j+1,i));
end
end
