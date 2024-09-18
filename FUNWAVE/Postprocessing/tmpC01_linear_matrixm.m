clear all

x2=[-76.291282 -76.270701];
y2=[36.947564 36.963229];

x1=[4984 7896];
y1=[830 2516];

a=diff(x2./y1)/diff(x1./y1);
b=diff(x2./x1)/diff(y1./x1);
c=diff(y2./y1)/diff(x1./y1);
d=diff(y2./x1)/diff(y1./x1);

X1d=[4984:7896];
Y1d=[830:2516];
[X1,Y1]=meshgrid(X1d,Y1d);


X2=a.*X1+b.*Y1;
Y2=c.*X1+d.*Y1;