function [S,revs,k]=Build_Spiral(spacing,Max_Radius,points,x0,y0)
revs=Max_Radius/spacing;                    %calc number of revolutions
k=spacing/(2*pi);                       
theta=linspace(0,revs*2*pi,points); %how many points and angle coord
r=k*theta;                       %radius coord
x=r.*cos(theta)+x0;                 %convert to x cartesian
y=r.*sin(theta)+y0;                 %convert to y cartesian
S=[x',y'];
end
