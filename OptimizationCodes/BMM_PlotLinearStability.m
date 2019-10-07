function [] = BMM_PlotLinearStability(D,A,R,xa,ya);
I=eye(size(R));
V=@(z)  inv(I-z*R)*(D+A*z);
R1=@(z) max(abs(eig(V(z))));


 
%xa=-300;
xb=1;
%ya=-300;
yb=-ya; 
nptsx = 301; 
nptsy = 301;

%Uniform sampling
%x = linspace(xa,xb,nptsx);
%y = linspace(ya,yb,nptsy);
%Fine to Coarse sampling
xa=abs(xa);ya=abs(ya);
x=2-logspace(0,log10(xa+2),100);
y=[fliplr(1-logspace(0,log10(ya+1),50)),logspace(0,log10(ya+1),50)-1+eps]
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;
 
[imax,jmax]=size(Z);
for i=1:imax
    for j=1:jmax;
    Rval1(i,j)=R1(Z(i,j));
    end
end
%----------------------------------------------------------------------
% Plot S, region of absolute stability:
 
Rabs1 = abs(Rval1);
figure()
view(2)
colormap([.7 .7 1;1 1 1])
%colormap([.7 .7 .7;1 1 1])   % grayscale
caxis([.99 1.01])
%daspect([1 1 1])
hold on
%contour(x,y,Rabs1,[1 1],'r','linewidth',3)
contourf(x,y,Rabs1,[0,1],['r'])
 
 
box on
%hold off
title('Region of absolute stability','FontSize',15)
% plot axes:
hold on
plot([-xa xb],[0 0],'k')
plot([0 0],[-ya ya],'k')
end

 
 