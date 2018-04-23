function [U_2] = FraunProp(U_1,dx_1,dy_1,lambda,z)
% Fraunhofer Propagator
k = 2*pi/lambda;
lz = lambda*z;

[Mx_1,Ny_1] = size(U_1);
Lx_1 = dx_1*Mx_1;
Ly_1 = dy_1*Ny_1;


Mx_2 = Mx_1;
Ny_2 = Ny_1;
Lx_2 = lz*Mx_2/Lx_1;
Ly_2 = lz*Ny_2/Ly_1;
dx_2 = Lx_2/Mx_2;
dy_2 = Ly_2/Ny_2;
x_2 = -Lx_2/2:dx_2:Lx_2/2-dx_2;
y_2 = -Ly_2/2:dy_2:Ly_2/2-dy_2;
[X_2,Y_2] = meshgrid(x_2,y_2);

fraunFactor = (exp(1j*k*z)/(1j*lz)).*exp(1j*k/(2*z)*(X_2.^2 + Y_2.^2));
U_2 = fraunFactor.*ifftshift(fft2(fftshift(U_1)))*dx_1*dy_1;

end

