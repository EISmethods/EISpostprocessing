
function[con,coneq,tau] = nlc_EIS_BMM(X,k,p,type);

[D,A,R,c]=unpackScheme(X,k,type);

%Stability Criteria

%If Optimizing for SSP
% r=-X(end);
% I=eye(s);
% P=inv(I+r*R)*(r*R);
% Q=inv(I+r*R)*(r*A);
% T=inv(I+r*R)*(D-r*A);
% con=-[P(:);Q(:);T(:)];

%If optimizing for Linear Stability
r=-X(end);
angle=[pi];
direction=(cos(angle)+sqrt(-1)*sin(angle));
scale=r*[linspace(eps,1,40)];
lam=kron(scale,direction); %Defines a Mesh of which the Stability function should be satisfied which is stretched by r

I=eye(k);
%V=@(z) inv(I-z*R)*(D+A*z) ;                     % Amplification/Iteration Matrix
EigV=@(z) eig(  inv(I-z*R)*(D+A*z) );            % Eigenvalues of Amplification Matrix

con=zeros(1,length(lam));

for kk=1:length(lam)
    con(kk)= max(abs( EigV(lam(kk))))-1;
end

con=[con(:)];
%Tau is the vector of the local truncation errors
tau=zeros(k,p);
e=ones(k,1);

for jj=1:p;
tau(:,jj)=(D*(c-e).^jj)/jj + A*(c-e).^(jj-1) + R*c.^(jj-1) - (c.^jj)/jj;
end

%Typically the order condition requires tau(1:p)=0, but EIS/super-convergence allows relaxed conditions

%coneq=[tau(:,1:p-1),D*tau(:,p)];  %p+1
coneq=[tau(:,1:p-2),D*tau(:,p-1:p),D*(R+A)*tau(:,p-1)];%p+2

end