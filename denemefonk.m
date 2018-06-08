function [eout] = denemefonk(K)

L=300e-6;               %Bobin
Rl=140e-3;              %Bobinin iç direnci.

C=470e-6;             %Kapasite
Rc=280e-3;              %Kapasitenin iç direnci.
Ryuk=10;               

a=Rc*Ryuk;
b=Ryuk/C;
c=L*Rc+L*Ryuk;
d=L/C+Rl*Rc+Rl*Ryuk+Rc*Ryuk;
e1=(Ryuk+Rl)/C;


xz1=0;
xz11=0;
xz2=0;

R=100;
ts=10e-6;
tsim=0.001;

c2=0;
et=0;


for i=0:ts:tsim
    
e=(R-c2);

%et=et+e^2;
%et=et+abs(e);
%et=et+i*abs(e);
et=et+i*e^4;

x1=(e)+xz1;
c1=ts*xz1*K(1);
xz1=x1;
u1=K(2)*e;
U=u1+c1;

x2=(U+2*c*xz11-c*xz2-d*xz11*ts+d*xz2*ts-e1*xz2*(ts^2))/c;
c2=(a*xz11*ts-a*xz2*ts+b*xz2*(ts^2));
xz2=xz11;
xz11=x2;
end
eout=et;
end

