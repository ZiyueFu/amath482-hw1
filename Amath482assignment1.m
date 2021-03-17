% Clean workspace
clear all; close all; clc

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata

L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); 
x = x2(1:n); 
y =x; 
z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; 
ks = fftshift(k);
u = sech(x);
Uave = zeros(n,n);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Uave = Uave + fftn(Un(:,:,:));
    
%     M = max(abs(Un),[],'all');
%     close all, 
%     isosurface(X,Y,Z,abs(Un)/M,0.7)
%     axis([-20 20 -20 20 -20 20]), grid on, drawnow
%     pause(1)
end
Uave = fftshift(Uave)./49;
[row,ind] = max(Uave(:));
[a,b,c] = ind2sub(size(Uave),ind);
xp = Kx(a,b,c);
yp = Ky(a,b,c);
zp = Kz(a,b,c);
center = [xp, yp, zp];

% Filter
tau = 30;
% k0 here is the center for each x, y, z.
filter = exp(tau*(((Kx-xp).^2)+((Ky-yp).^2)+((Kz-zp).^2)));

for i = 1:49
    Un(:,:,:)=reshape(subdata(:,i),n,n,n);
    Uft = fftn(Un(:,:,:)); 
    Usft = fftshift(Uft);
    Ufilt = filter.*Usft;
    Unfilt = ifftn(Ufilt);
    
    [row, ind] = max(Unfilt(:));
    [b,a,c] = ind2sub(size(Unfilt),ind);
    A(i) = a;
    B(i) = b;
    C(i) = c;
end

plot3(x(A),y(B),z(C))
xlabel('x')
ylabel('y')
zlabel('z')
% averageX = max((Uave(:)));
% x_coor = averageX(1,1);
% averageY = max(abs(Uave(1,:,1)));
% y_coor = averageY(1,1);
% averageZ = max(abs(Uave(1,1,:)));
% z_coor = averageZ(1);
% center = [x_coor, y_coor, z_coor];
% Filter
% tau = 30;
% k0 = 0;
% filter = exp(-tau*(k-k0).^2);
% % plot3(fftshift(ks), fftshift(filter),x,fftshift(ks), fftshift(filter),y,fftshift(ks), fftshift(filter),z);
% title('Trajectory of Submarine')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% set(gca, 'Fontsize',10)


% hold on
% plot3(fftshift(ks), fftshift(filter),y)
% plot3(fftshift(ks), fftshift(filter),z)




%     M = max(abs(Utave),[],'all');
%     isosurface(Kx,Ky,Kz,abs(Utave)/M,0.7)
% for j = 1:length(x)
%     Ux(j,:,:) = fft(Un(j,:,:));
%     Uy(:,j,:) = fft(Un(:,j,:));
%     Uz(:,:,j) = fft(Un(:,j,:));
%     Kp(j,:,:) = fftshift(Kx(j,:,:));
%     Up(j,:,:) = fftshift(Ux(j,:,:));  
    
%     avg = abs(Ux)./max(abs(Ux))
%     isosurface(Ux, Uy, Uz, abs(Un)/M, 0.7)
% end

% slice = 0:0.5:24
% for j=1:length(X)
%     Ux(j,:,:) = fft(X(j,:,:));
%     KX(j,:,:) = fftshift(Kx(j,:,:));
%     Ux(:,j,:) = fft(Y(:,j,:));
%     KX(:,j,:) = fftshift(Ky(:,j,:));
% end
% 
% % isosurface(K, Ux, abs(Ut)/max(abs(Ut(1,:))))
% xlabel('x-axis'), ylabel('y-axis'), zlabel('z-axis')
% Uave = Uave/length(slice);
% plot(x, abs(Uave)/max(abs(Uave)),'k')
% hold on
% plot(ks, abs(fftshift(Ut(1,:))/max(abs(Ut(1,:)))),'k','Linewidth',2)

% Without get the center frequency, I can only assume the tau and k0
% tau = 10;
% k0 = 0;
% filter = exp(-tau*(k-k0).^2);
% plot3(fftshift(Kx), abs(fftshift(ut))/max(abs(fftshift(ut))))




