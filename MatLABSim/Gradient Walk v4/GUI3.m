function c=GUI3(revs,k,x,y,X,Y,Z)
%% Nice Spiral Graphing for Presentaions
theta=linspace(0,revs*2*pi,3000); %how many points
r=k*theta;                   %spiral
XX=r.*cos(theta);             %convert to x cartesian
YY=r.*sin(theta);             %convert to y cartesian

figure(4)
hold on
plot(x,y,'r.') 
plot(XX,YY) 
axis equal
grid on
axis off
%% Function with Spiral Overlay for Presentaions
figure(5)
hold on
surf(X,Y,Z)
axis off
plot(x,y,'r.')


shading interp

end