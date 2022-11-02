function c=GUI1(x,y,X,Y,Z)
figure(1)
subplot(2,1,1)
plot(x,y,'r.') 
axis equal
grid on
axis off
title('Sampling Spiral')

subplot(2,1,2)
hold on
surf(X,Y,Z)
axis off
plot(x,y,'r.')
title('Function and Sampling Spiral')
shading interp
end