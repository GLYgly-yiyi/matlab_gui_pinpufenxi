function[H,W]=dtft(h,N) 
N=fix(N); 
L=length(h);
h=h(:); 
W=(2*pi/N)*[0:N-1]'; 
mid=ceil(N/2)+1; 
W(mid:N)=W(mid:N)-2*pi; 
W=fftshift(W); 
H=fftshift(fft(h,N)); 

