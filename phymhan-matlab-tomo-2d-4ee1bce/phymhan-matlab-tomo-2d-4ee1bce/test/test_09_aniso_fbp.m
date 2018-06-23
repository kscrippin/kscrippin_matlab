%TEST 09: Filtered Back Projection using anisotropic filter
clear
clc
fprintf('TEST 03: Filtered Back Projection using anisotropic filter\r');
angles = 0:1:179.999;
%Phantom
im = imtest('phan',256);
[M,N] = size(im);
D = ceil(sqrt(M^2+N^2));
M_pad = ceil((D-M)/2)+1;
N_pad = ceil((D-N)/2)+1;
figure('name','Original Image: Phantom')
imshow(im)
[projmat,~] = tomoproj2d(im,angles);
fprintf('ANISOFBP  method\r')
im_rec1 = tomo_recon_anisofbp(projmat,angles);
im_rec1 = im_rec1(M_pad:D-M_pad,N_pad:D-N_pad);
im_rec1 = uint8(imscale(im_rec1));
figure('name','ANISOFBP')
imshow(im_rec1);
fprintf('FBP method\r')
im_rec2 = tomo_recon_fbp(projmat,angles);
im_rec2 = im_rec2(M_pad:D-M_pad,N_pad:D-N_pad);
im_rec2 = uint8(imscale(im_rec2));
figure('name','FBP')
imshow(im_rec2);
