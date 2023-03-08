function c=GUI2(x,y,Intensity) %3D dots and lines of spiral mapped intensity,no axis
figure(2) 
hold on
plot3(x,y,Intensity,'r.')
plot3(x,y,Intensity)
%axis off
colormap('jet(200)')
title('Spiral Mapping')
xlabel('X Translation')
ylabel('Y Translation')
zlabel('Intensity')

%axis equal
shading interp
end