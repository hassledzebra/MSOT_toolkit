


A = recon.model.ModelOperator;

bMod = magic(5);
B = recon.model.ModelOperator; B.Model = bMod;

cMod = magic(6);
cMod = cMod(1:5,1:3);
C = recon.model.ModelOperator; C.Model = cMod;

b1 = 5;
b2 = 1:100;

c1 = A*b1
c3 = A(b1)


c2 = A*b2
c4 = A(b2)



D = B*C

e1 = randn(5,1);
e2 = randn(3,1);

D*e2
e1.'*D*e2


% c2 = b1*A;


