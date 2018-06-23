function [W, p, t] = build2(im, angles)
t_s = toc;
%Pad image
[im_pad,D] = impad(im,'diag',0);
sz = size(im);
m = sz(1); %rows
n = sz(2); %cols
n_proj = length(angles);
N = m*n; %number of unknown variables
M = D*n_proj; %number of equations
try
    W = zeros(M,N); %weighting factor matrix
catch expr
    fprintf([expr.message '\nGenerating a sparse...\r'])
    W = sparse(M,N);
end
%proj_info = zeros(M,2); %[theta No.]
p = zeros(M,1); %projection vector
ix = false(M,1);
rc = [m;n]/2+0.5;
%cnt = 0;
fprintf('Building Weight Matrix...\r')
[x,y] = ind2sub(sz,1:N);
%xy = [x;y]; %x,y coordinate
%x,y coordinate with respect to rotation center
xy_c = [x;y]-repmat(rc,1,N);
for kp = 1:n_proj
    im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
    pvec = sum(im_rot,1);
    t = -angles(kp)/180*pi;
    R = [cos(t) -sin(t);sin(t) cos(t)];
    fprintf('\nAngle No.%d(%d Degree)\r',kp,angles(kp))
    xy_rot = R*xy_c+D/2+0.5;
    idx = round(xy_rot(2,:));
    ixM = D*(kp-1)+idx; %corresponding indice in W and p matrix
    for kn = 1:N
        W(ixM(kn),kn) = 1;
        p(ixM(kn)) = pvec(idx(kn));
        ix(ixM(kn)) = true;
    end
end
% %
W = sparse(W);

%Delete all-zero rows in W
fprintf('\nDelete all-zero rows in W...\r')
%ix = sum(W,2)==0;
W(~ix,:) = [];
p(~ix) = [];
n_eq = size(W,1);
fprintf('\nFinish building A. \n%d equations in total.\r',n_eq);
figure('name','Sparsity pattern of matrix W');
spy(W)
t_e = toc;
t = t_e-t_s;
fprintf('\nTotal: %f seconds\n\r',t);
end
