function u = convolve_fft(t,F,params)
beta = params.beta;
A = params.A;
B = params.B;
G0 = params.G0;
N = length(t);
T = t(end)-t(1);

om = (2*pi/T)*[0:N/2-1 0 -N/2+1:-1]; 
% This form for omega depends on Matlab's conventions for the FFT 
% mode order - look up "spectral differentiation in MATLAB"
% Also, it has an extra zero in it that you wouldn't expect to suppress
% some noise when doing differentiation. This is only necessary for
% an even number of modes and might not be right. 

%Kom = 1./(G0*(1i*om/w0).^beta + G0); % G0 = 4e2 pN/micron^3
%Kom(abs(om)<100*eps) = 0; % remove the zero mode -- note 
                          % that this could alternately be handled
                          % by assuming G(omega) = G_0 + (1i*om)^beta
                          
G1 = @(w) A*cos(pi*beta/8)*w.^beta + B*cos(3*pi/8)*w.^0.75; % form given in hoffman crocker paper
G2 = @(w) A*sin(pi*beta/8)*w.^beta + B*sin(3*pi/8)*w.^0.75; % ^

G = @(w) G1(w) + 1i*G2(w); % ^ DONT TAKE MAGNITUDE

rescale = 160/norm(G(1000)); % fitting to experimental results

G = @(w) rescale*(G1(w) + 1i*G2(w)) + G0; % Adding in a small constant to account for zero mode, should be unobservable to experiments
Kom = 1./G(om);

u = ifft(Kom.*fft(F));