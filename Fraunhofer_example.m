pp=1.36e-6;
z=1;
lambda = 632.80e-09;

A = imread('DOE_big.bmp');
A= A(:,:,1)/84;
Phase=(double(A)-2) *pi/2;
Amp=ones(size(Phase))*50;
E0= Amp.* exp(1i.*Phase);

%E=fresnel_advance (E0, pp, pp, z, lambda);
E=FraunProp(E0,pp,pp,lambda,z);

figure, imshow(abs(E),'Border','tight');

