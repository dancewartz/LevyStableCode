function [Fxs,Fys,Fzs,ts] = tensor_sum(params)
%% Importing parameters from the params struct
tmax = params.tmax;
dt = params.dt;
Npoints = params.Npoints; 
L = params.L;
m_monopole = params.m_monopole;
d = params.d;
tau = params.tau; % sets our timescale -- the random switching time
ts = linspace(0,tmax,round(tmax/dt));
F = params.F;
poisson_ratio = params.poisson_ratio; % For incompressible field 

%% Setting up simulation

% Trajectory coordinates
Fxs = zeros(size(ts));
Fys = zeros(size(ts));
Fzs = zeros(size(ts));

% next time a force turns on, how long it stays on
nexton = -(tau)*log(rand(1,Npoints));
durations = -(tau)*log(rand(1,Npoints));

% positions of the centers of the dipoles 
xs = L*(rand(1,Npoints)-0.5);
ys = L*(rand(1,Npoints)-0.5);
zs = L*(rand(1,Npoints)-1-d/L); % check this TODO 
rs = sqrt(xs.^2+ys.^2+zs.^2);

% Creating the uniform orientation vectors b
phis = 2*pi*rand(1,Npoints);
thetas = acos((2*rand(1,Npoints)-1));
bx = sin(thetas).*cos(phis);
by = sin(thetas).*sin(phis);
bz = cos(thetas);

% Each index of ux_ contains the respone from the force in the _ direction
uxc = zeros(1,Npoints);
uyc = zeros(1,Npoints);
uzc = zeros(1,Npoints);

tic;

%% Simulatation Löõôop
for s = 1:length(ts)
turnon = (ts(s) - nexton <= dt) & (ts(s) - nexton > 0); % dipoles that are turning on this time step
turnoff = (ts(s) >= nexton+durations); % dipoles that are turning off on this time step
status = turnon - turnoff; % status change of dipole: 1 for turned on during this step, -1 for turned off during this step
ind_change = find(status); % find indicies of all the dipoles that changed
num_change = length(ind_change); % number of dipoles that switched on this step (either on or off)

% Positions of one force in the dipole
xs1 = xs(ind_change) + bx(ind_change)*d/2;
ys1 = ys(ind_change) + by(ind_change)*d/2;
zs1 = zs(ind_change) + bz(ind_change)*d/2;

% Positions of the other force in the dipole
xs2 = xs(ind_change) - bx(ind_change)*d/2;
ys2 = ys(ind_change) - by(ind_change)*d/2;
zs2 = zs(ind_change) - bz(ind_change)*d/2;


for k = 1:num_change
    i_k = ind_change(k); % index of kth dipole to change status
    F_vector = F*[bx(i_k) by(i_k) bz(i_k)]'; % force fector

    position = [0 0 0]'; % we measure response at the origin
    position_prime_1 = [xs1(k) ys1(k) zs1(k)]'; % position of first force
    position_prime_2 = [xs2(k) ys2(k) zs2(k)]'; % position of the second force

    % Greens function, G1_image is used to compare to plots in the Schwarz
    % paper of dipoles in an elastic half space
    [G1,G1_image] = Greens_Function(position,position_prime_1,params,true); 
    [G2,G2_image] = Greens_Function(position,position_prime_2,params,true);

    % Status tells us whether we add or sub. Calculates the response due to
    % the dipole.
    response = status(i_k)*(G2 - G1)*F_vector; % TODO check this

    % Update the random walk with the new response
    uxc(i_k) = uxc(i_k) + response(1);
    uyc(i_k) = uyc(i_k) + response(2);
    uzc(i_k) = uzc(i_k) + response(3);
end

% Add up contributions due to each dipole
Fxs(s) = sum(uxc);
Fys(s) = sum(uyc);
Fzs(s) = sum(uzc);
   
% update all the state of those that switch
% need to change x,y,z,r, Fmag, duration,nexton
usi = size(nexton(turnoff)); % update size
nexton(turnoff) = ts(s) + (-tau)*log(rand(usi)); % for the forces that turned off, find out when they turn on again
durations(turnoff) = (-tau)*log(rand(usi)); % Find out how long those next turn on events last
xs(turnoff) = L*(rand(usi)-0.5); % Redistribute forces that turned off
ys(turnoff) = L*(rand(usi)-0.5); % ^ 
zs(turnoff) = L*(rand(usi)-1-d/L); % ^^
rs(turnoff) = sqrt(xs(turnoff).^2+ys(turnoff).^2+zs(turnoff).^2);
phis = 2*pi*rand(usi); % Reorient forces that turned off
thetas = acos((2*rand(usi)-1));
bx(turnoff) = sin(thetas).*cos(phis);
by(turnoff) = sin(thetas).*sin(phis);
bz(turnoff) = cos(thetas);

% give an update
if(rem(s,1e5)==0)
    fprintf('s = %d, %3.3g percent done \n',s,100*s/length(ts));
    toc;
end


end

% this is a good heuristic check, expect mean is small ish.....
mean(Fxs)
mean(Fys)
mean(Fzs)