function varargout = NOZZLEfp(action,varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Problem : NOZZLEfp
%    ---------
%
%    The problem is function is a simulator for in impingement cooling system of a 
%    nozzle in a gas turbine. It is defined for a constrained Black-Box 
%    optimization where the objective is to maximize the heat transfer 
%    coefficient of the system by acting on the design variables of the 
%    impingement plate.
%
%    Source: 
%    Cocchi, L., Marini, F., Porcelli, M., & Riccietti, E. (2024) 
%    Black-box optimization for the design of a jet plate for impingement cooling.
%    https://optimization-online.org/?p=25814
%
%    S2MPJ input: Filippo Marini, 8 X 2024.
%
%    classification = 'N-MOOR0-RN-5-11'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = 'NOZZLEfp';
switch(action)
case 'setup'
    pb.name      = name;
    pb.n = 5;
    pb.m = 11;
    pb.xlower    = [5e-4;5e-4; 5e-4; 2e-4; -Inf];
    pb.xupper    = [1e-2;1e-2; 1e-2; 1e-2; Inf];
    pb.x0        = [ 5e-3; 4e-3; 3e-3; 1e-3; 1];
    pb.xtype     = ['r';'r';'r';'r';'b'];
    pb.objderlvl = 0;
    pb.pbclass   = 'N-OOR0-RN-5-11';
    varargout{1} = pb;
    
case 'fx'
    [varargout{1},~] = IMP_CONS_fixed_params(varargin{1});
case 'cx'
    [~,varargout{1}] = IMP_CONS_fixed_params(varargin{1});
otherwise
    disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HTC,CONS] = IMP_CONS_fixed_params(x)

% Derived from the NOZZLE Black-Box function defined in
% 
%    Cocchi, L., Marini, F., Porcelli, M., & Riccietti, E. (2024) 
%    Black-box optimization for the design of a jet plate for impingement cooling.
%    https://optimization-online.org/?p=25814
%
% This matlab function returns the Root Mean Square of the distribution 
% of the Heat Transfer Coefficient for an impingement jet plate starting 
% from the design variables. It also returns the value of all constraints
% for the constrained Black-Box Optimiztion problem.
% The parameters derived from the bounduary conditions are fixed.
% Full description below
%
%------------------------------INPUT---------------------------------------
%__________VARIABLES______________________________________________________
% x(1)      = hole spacing along the x direction
% x(2)      = hole spacing along the y direction
% x(3)      = distance between the jet plate and the target surface
%             (i.e. meatus width)
% x(4)      = diameter of the holes
% x(5)      = boolean variable for the hole pattern of the jet plate, 1
%             for in-line pattern and 0 for staggered pattern
%__________PARAMETERS FIXED INSIDE THIS FUNCTION (SEE BELOW)____________
% Params       = structure array containing the scalar parameters needed,
%                containing the following fields

% ---------parameters for HTC and temperatures-----------------------------

% Params.Lx      = [m] x-length of the jet plate
% Params.Ly      = [m] y-length of the jet plate
% Params.Tc      = [K] inlet temperature of the cooling air
% Params.Pc      = [Pa] inlet pressure of the cooling air
% Params.Tg      = [K] temperature of the external hot gas
% Params.hg      = [W/(m^2*K)] HTC of the external hot gas
% Params.M       = [Kg/s] Inlet mass flow of cooling air
% Params.Cd      = Discharge coefficient of every jet
% Params.kw      = [W/(m*K)] Thermal conductivity of the target surface
% Params.ds      = [m] Thickness of the target surface
% --------parameters for inner iteative procedures-------------------------
% Params.maxit_h = maximum number of iterations in the T_f estimation process
% Params.tol_h   = [W/(m^2*K)] Tolerance on the distance between two susequent
%                   HTC distributions
% Params.maxit_P = maximum number of bisection iteration to estimate the
%                  outlet pressure of the cooling air
% Params.tol_P   = [Pa] Tolerance  for bisection method
%
%----------parameters for the constraints----------------------------------
%
% Params.DTmax   = [K] upper bound for the temperature gradient between the
%                  internal and external wall of the target surface
% Params.Tmax    = [K] upper bound for the external wall temperature
% Params.dp_p_max = upper bound for pressure ratio

%----------------------------OUTPUT----------------------------------------
% HTC  = normalized norm-2 of vector of length N_x containing the streamwise 
%        distribution of the HTC of the coolant for every spanwise row
% CONS = array of lenght 11 containing the value of all the constrained
%        functions of the problem of optimiztion for a jet impingement 
%        cooling system
%--------------------------------------------------------------------------
% Check on design variables
if numel(x) ~= 5
    error('NOT ENOUGH INPUT ARGUMENTS')
elseif (x(5) ~= 1 && x(5) ~= 0)
    error('x(5) MUST BE EQUAL TO 1 OR 0.')
else
    xn = x(1); yn = x(2); zn = x(3); d = x(4); pattern = x(5);
end

% Fixing Parameters______________________________

Params.Lx = 5e-2;    % [m] Streamwise length of the impingement plate
Params.Ly = 5e-2;    % [m] Spanwise width of the impingement plate
Params.Tc = 773.15;  % [K] inlet temperature of the cooling air
Params.Pc = 1e6;     % [Pa] inlet pressure of the cooling air
Params.Tg = 1273.15; % [K] temperature of the external hot gas
Params.hg = 1000;    % [W/(m^2*K)] HTC of the external hot gas
Params.M = 1e-2;     % [Kg/s] total mass flow
Params.Cd = 0.85;    % Discharge coefficient
Params.kw = 20;      % [W/(m*K)] Thermal conductivity of the target surface
Params.ds = 3e-3;    % [m] Thickness of the target surface
Params.maxit_h = 100;% maximum number of inner iterations in the T_f 
                     % estimation process
Params.tol_h = 1e-8; % [W/(m^2*K)] Tolerance on the distance between two 
                     % susequent HTC distributions
Params.maxit_P = 100;% Maximum number of bisection iteration to estimate the
                     % outlet pressure of the cooling air
Params.tol_P = 1e-6; % [Pa] Tolerance  for bisection method
Params.DTmax = 60;   % [K] upper bound for the temperature gradient between
                     % the internal and external wall of the target surface
Params.Tmax  = 1073; % [K] upper bound for the external wall temperature
Params.dp_p_max = 4e-2; % upper bound for pressure ratio

Aux_arrays = get_aux_arrays; % get arrays for the interpolation of physical 
                             % quantities for the cooling air

Nx = floor(Params.Lx/xn); % Number of spanwise jet rows
Ny = floor(Params.Ly/yn); % Number of holes in every spanwise row

x_vec = (0.5*xn:xn:(Nx-0.5)*xn)'; % x coordinate of the centers of the jet holes
%-------------------------
% Flow rate distributions |
%-------------------------
% Definition of some auxiliary vectors
Ny_vec = Ny * ones(Nx,1);
d_vec = d * ones(Nx,1);
% Definition of the model by Florschuetz

Gj_mod = @(x,m) Florschuetz_Gj(x,xn,yn,zn,d,Params.Cd,Nx,Ny_vec,m);

% Evaluation of the streamwise distributions of the jet mass velocity for
% every spanwise row (Gj_vec) and of the ratio of the crossflow mass
% velocity to the jet mass velocity for every spanwise row (Gc_on_Gj_vec)

[Gj_vec, Gc_on_Gj_vec] = flow_rate_dist_2(Params.M,x_vec,d_vec,...
    Ny_vec,Params.Ly,zn,Gj_mod);

% Interpolation of the dynamic viscosity of the cooling air

mu = interp1(Aux_arrays.Tc,Aux_arrays.mu,Params.Tc,'pchip'); 

% Interpolation of the Prandtl number of the cooling air

Pr = interp1(Aux_arrays.Tc,Aux_arrays.Pr,Params.Tc,'spline'); 

% Streamwise distribution of the Reynolds number of the spanwise jet rows

Rej_vec = (d * Gj_vec)/mu;

% Streamwise distribution of the Nusselt number of the spanwise jet rows

Nu = Florschuetz_Nusselt(xn,yn,zn,d,Pr,Rej_vec,Gc_on_Gj_vec,pattern);

% Iterative estimation of the HTC of the cooling air
% Initialization
% Thermal conductivity of the cooling air

k_old = interp1(Aux_arrays.Tc,Aux_arrays.k,Params.Tc,'pchip'); 

% Streamwise distribution of the HTC of the coolant for every spanwise row

h_dist_old = Nu*k_old/d;
% Internal and external wall temperature

[Tint_dist_old,~] = Wall_Temp_dist(h_dist_old,Params.hg,xn,Params.Ly,Params.ds,...
    Params.kw,Params.Tg,Params.Tc);

for it = 1:Params.maxit_h
    % Film temperature
    T_f = 0.5*(Tint_dist_old + Params.Tc);
    % New thermal conductivity of the cooling air
    k = interp1(Aux_arrays.Tc,Aux_arrays.k,T_f,'pchip');
    % New HTC of the cooling air
    h_dist = Nu.*k/d;
    [Tint_dist,Text_dist] = Wall_Temp_dist(h_dist,Params.hg,xn,Params.Ly,Params.ds,...
        Params.kw,Params.Tg,Params.Tc);
    HTC_dist = norm(h_dist-h_dist_old)/norm(h_dist_old);

    if HTC_dist <= Params.tol_h
        break
    else
        h_dist_old = h_dist;
        Tint_dist_old = Tint_dist;
    end

end
% Objective value
HTC = norm(h_dist,2)/sqrt(numel(h_dist));

DeltaT = norm(Text_dist - Tint_dist)/sqrt(numel(Text_dist));
Text = norm(Text_dist)/sqrt(numel(Text_dist)); 

% Assembling the penalty function

% Evaluation of the outlet pressure of the cooling air
% Interpolation of the isentropic expansion factor of the cooling air 
gamma = interp1(Aux_arrays.Tc,Aux_arrays.gamma,Params.Tc,'pchip'); 
f_press = @(x) mass_flow_diff(x,Params.Pc,Params.Cd,Gj_vec(Nx),Params.Tc,gamma);
P_star = Params.Pc * (2/(gamma+1))^(gamma/(gamma-1));
[P_out,~,~,~] = my_bisec(f_press,P_star, Params.Pc,Params.tol_P,Params.maxit_P,0);


% LEQ CONSTRAINTS
CONS = zeros(11,1);                             
CONS(1) = zn - 3*d;                             % Validity constraint
CONS(2) = yn - 8*d;                             %    ''        '' 
CONS(3) = xn - 3.75*yn;                         %    ''        ''  
CONS(4) = xn - d*(15*pattern + 10*(1-pattern)); %    ''        '' 
CONS(5) = (DeltaT - Params.DTmax)/Params.DTmax; % Temperature constraint
CONS(6) = (Text - Params.Tmax)/Params.Tmax;     %    ''        '' 
%------------------------------------------------
% GEQ CONSTRAINTS
CONS(7) =  zn - d;                              % Validity constraint
CONS(8) =  yn - 4*d;                            %    ''        ''
CONS(9) =  xn - 0.625*yn;                       %    ''        ''
CONS(10) = xn - 5*d;                            %    ''        '' 
CONS(11) = (P_out - Params.Pc*Params.dp_p_max)/Params.Pc;% Pressure constraint

end

%_______________DEPENDENCIES________________________________________________

function [Gj]=Florschuetz_Gj(x,xn,yn,zn,d,C_d,N_x,N_y,M)

% Function that evaluates the jet mass velocity of a spanwise row (Gj)
% depending on the streamwise position of the centers of the holes
% of the spanwise row. Based on the model by Florschuetz et al. 
%
% "Jet array impingement with crossflow-correlation of streamwise 
% resolved flow and heat transfer distributions". 1981. NASA Tech. Report

% INPUT
% x   = streamwise position of centers of the holes of the spanwise row. Must
%       be a multiple of the half-partition xn/2. COULD BE ALSO A VECTOR!
% xn  = streamwise jet hole spacing
% yn  = spanwise jet hole spacing
% zn  = channel height or jet plate-to-impingement surface spacing
% d   = jet hole diameter
% C_d = jet plate discharge coefficent
% N_x = number of spanwise rows
% N_y = vector of lenght N_c containing the number of holes in every spanwise row
% M   = total flow rate provided to the cooling system
%
% OUTPUT
% Gj = streamwise jet mass velocity distribution

Aj = 0.25*pi*d^2;
if isscalar(N_y) 
    Aj_tot = Aj*N_x*N_y;
else
    Aj_tot = Aj*sum(N_y);
end
delta = Aj*C_d*sqrt(2)/(yn*zn);
gamma =delta*N_x*M./Aj_tot;
b = delta*x/xn;
Gj = gamma.*cosh(b)/sinh(delta*N_x);

end
%---------------------------------------------------------------------------

function [Gj_vec,Gc_on_Gj] = flow_rate_dist_2(m,x_h,d_vec,N_s,Ly,zn,Gj_mod)

% Function that evaluates the flow rate disribution among all the N_c spanwise
% rows of a jet impingement array. 
% The inputs we need are the geometrical parameters of the jet impingement 
% array, the initial guess for the flow rate at the most upstream spanwise
% row and a mathematical model for the jet mass velocity (Gj) for every spanwise row.
% INPUT:
% m      = total flow rate provided upstream of the jet impingement plate.
% x_h    = vector of length N_c containing the x-coordinates of the centers 
%          of the spanwise rows.  
% d_vec  = vector of length N_c containing the diameters of the holes of
%          every spanwise row.
% N_s    = vector of length N_c containing  the number of holes
%          for every sapnwise row.
% Ly     = y-direction width of the impingement plate.
% zn     = channel height or jet plate-to-impingement surface spacing.
% Gj_mod = mathematical model for the jet mass velocity (Gj) for every spanwise row. 
%         [FUNCTION HANDLE THAT MANAGES VECTORS VARIABLES]
%          Must be defined outside flow_rate_dist_2 as 
%          Gj_mod = @(x,m)Florschhuetz_Gj(x,xn,yn,zn,d,C_d,N_c,N_s,m)
% OUTPUT:
% Gj_vec   = vector of length N_c where the i-th element is the mass velocity
%            of a single jet hole on the i-th spanwise row. 
% Gc_vec   = vector of length N_c where the i-th element is the mass velocity
%            of the crossflow upstream of the i-th spanwise row. 
% Gc_on_Gj = ratio of the crossflow mass velocity (Gc) to the jet mass 
%            velocity (Gj) for every spanwise row
%--------------------------------------------------------------------------
% Creating a vector with the areas of the jet holes in every spanwise row
A_vec = 0.25*pi*d_vec.^2;
C_sec = Ly*zn; % Channel cross-section
% Initialization of vectorial outputs
N_c = length(x_h);

Gc_vec = zeros(N_c,1);
Gj_vec = Gj_mod(x_h,m);

% Evaluation of the distribution of the ratio Gc/Gj along x direction

for i = 2:N_c
    Gc_vec(i) = sum(Gj_vec(1:i-1).*A_vec(1:i-1).*N_s(1:i-1))/C_sec;
end
Gc_on_Gj = Gc_vec./Gj_vec;
end

%----------------------------------------------------------------------------------

function Nu = Florschuetz_Nusselt(xn,yn,zn,d,Pr,Re_j,Gc_on_Gj,pattern_flag)

% Function that evaluates the streamwise distribution of the Nusselt number
% following the correlation function developed by Florschuetz et al.
% INPUT:
%           xn = streamwise jet hole spacing (must be a scalar) 
%           yn = spanwise jet hole spacing (must be a scalar)
%           zn = channel height (jet plate-to-impingement surface spacing, must be a scalar)
%            d = jet hole diameter (scalar)
%           Pr = Prandtl number 
%         Re_j = jet Reynolds number streamwise distribution (Re_j = Gj*d/mu)
%                (vector)
%     Gc_on_Gj = streamwise distribution of the ratio of the crossflow mass velocity 
%                (Gc) to the jet mass velocity (Gj) (vector)
% pattern_flag = boolean variable for the hole pattern of the jet plate, 1
%                for in-line pattern and 0 for staggered pattern

% OUTPUT:
%           Nu = Nusselt number distibution for every spanwise row 

if nargin < 8
    pattern_flag = 1; % If not specified, the hole pattern is considered
end                   % in-line.

% Calculation of the coefficient A, B and of the exponents alpha and beta,
% always according to Florschuetz models.

switch pattern_flag
    case 0 % Staggerd hole pattern
        % 'C' coefficients C(1) is for A, C(2) is for alpha, C(3) is for B
        % and C(4) is for beta.
        C = [1.87 0.571 1.03 0.442];
        %---------------------------------
        % Exponents for the  four coeffcients A, alpha, B and beta
        A_exp = [-0.771 -0.999 -0.257];
        alpha_exp = [0.028 0.092 0.039];
        B_exp = [-0.243 -0.307 0.059];
        beta_exp = [0.098 -0.003 0.304];
        %---------------------------------
    case 1 % In-line hole pattern
        
        C = [1.18 0.612 0.437 0.092];
        
        A_exp = [-0.944 -0.642 0.169];
        alpha_exp = [0.059 0.032 -0.022];
        B_exp = [-0.095 -0.219 0.275];
        beta_exp = [-0.005 0.599 1.04];
        %---------------------------------
end

A = C(1)*prod([xn/d yn/d zn/d].^A_exp);
alpha = C(2)*prod([xn/d yn/d zn/d].^alpha_exp);
B = C(3)*prod([xn/d yn/d zn/d].^B_exp);
beta = C(4)*prod([xn/d yn/d zn/d].^beta_exp);

aux_1 = A*Re_j.^alpha;

aux_2 = 1 - B*((zn/d)*Gc_on_Gj).^beta;

Nu = aux_1.*aux_2*Pr^(1/3);
end

%-----------------------------------------------------------------------------


function [T_int, T_ext] = Wall_Temp_dist(h_c,h_g,xn,Ly,ds,k,T_h,T_c)

% Function that gives the streamwise distribution of the internal and the 
% external wall temperature. Starting from the distribution of the HTC of 
% the cooling air h_c we solve an heat equations modeling .
% INPUT:
% h_c = vector of length N_x containing the streamwise distributions of 
%       the HTC of the cooling air.
% h_g = HTC of the external hot gas.
% xn  = streamwise (along x-direction) pitch of the impingement holes.
% Ly  = spanwise width of the impingement plate.
% ds  = thickness of the target wall.
% k   = thermal conductivity of the target wall.
% T_g = temperature of the target wall.
% T_c = upstream temperature of the cooling air.

% OUTPUT:
% T_int = vector of length N_x containing the streamwise distribution 
%         of the internal wall temperature
% T_ext = vector of length N_x containing the streamwise distribution 
%         of the external wall temperature

% First we solve a discretized heat equation building a linear system 
% A * x = b

n = numel(h_c);
h_c = h_c(:);
S = Ly*xn; % Area of the surface of the 'element' of the nozzle wall in contact with 
           % the cooling air and the hot gas

Ar = Ly*ds; % Area of the contact surface between two elemnts of the nozzle wall

one_R_h = S/(1/h_g + 0.5*ds/k); % Inverse of convective thermal resistance 
                                % of the heat flow of the hot gas

one_R_c = S./(1./h_c + 0.5*ds/k); % Inverse of convective thermal resistance 
                                  % of the heat flow of the cooling air

one_R_i = (Ar*k)/xn; % Inverse of conduction thermal resistence

% Building the matrix

v = one_R_h + one_R_c;

v(1) = v(1) + one_R_i;
v(end) = v(end) + one_R_i;
v(2:end-1) = v(2:end-1) + 2*one_R_i;

w = one_R_i * ones(n,1);

A = spdiags([-w v -w],-1:1,n,n);

% Building rhs

b = one_R_h*T_h + one_R_c*T_c;

% Solving the linear system
T = A\b;

% T is the temperature distribution internal to the nozzle wall, now we use
% it to estimate the tempearture distributions onn the internal and
% external nozzle wall

a = 2*k/ds;

T_int = (a*T + h_c.*T_c)./(a + h_c);

T_ext = (a*T + h_g*T_h)./(a + h_g);

end
%--------------------------------------------------------------------------

function y = mass_flow_diff(P1,P0,Cd,G,T0,gamma)

% Auxiliary function that evaluates the difference between the outlet mass 
% flow velocity of the cooling air G derivated from the mass flow 
% distribution computation and the outlet mass flow velocity derivated from 
% the differnce P0-P1 through the isentropic mass flow formula; where P0 is
% the inlet pressure of the cooling air upstream the impingement plate and
% P1 is the outlet pressure of the cooling air downstream the impingement plate
%
% INPUT 
% P1 [Pa]      = outlet pressure of the cooling air 
% P0 [Pa]      = inlet pressure of the cooling air
% Cd []        = discharge coefficent of the impingement plate holes
% G [Kg/(m^2s)]= mass flow velocity of the cooling air at the last spanwise row of
%                holes
% T0 [K]       = inlet temperature of the cooling air
% gamma []     = isentropic exspansion factor of the cooling air
%
% OUTPUT
% y = difference of massflow

R = 8.314462; % [J/(K mol)] ideal gas constant 

alpha = G^2 * (gamma - 1) * R * T0;
alpha = alpha/(2 * Cd^2 * gamma);
alpha = sqrt(alpha);

a = 1/gamma;
b = 1 - a;
c = 0.5*b;
y = (P1.^a) .* sqrt(P0^b - P1.^b).*(P0.^c) - alpha;
end

%--------------------------------------------------------------------------

function [x,fval,k,flag] = my_bisec(f,x0,x1,tol,kmax,display)

% Function for very simple bisection method
flag = 0;
f0 = f(x0);
f1 = f(x1);
if f0*f1>0
    warning('Impossible to estimate correctly Pout')
    x = x1;
    fval = f1;
    k = 0;
    flag = -1;
    return
else

    for i = 1:kmax
        k = i;
        x2 = (x0 + x1)/2;
        f2 = f(x2);
        err = abs(x0-x1);
        if display == 1
            fprintf('x = %0.4e \t f(x1) = %0.4e \t err = %0.4e\n',x2,f2,err);
        end

        if f0*f2<0
            x1 = x2;
        else
            x0 = x2;
        end
        x2 = (x0 + x1)/2;
        err = abs(x2-x1);
        if err < tol
            x = x2;
            fval = f(x);
            flag = 1;

            if display == 1
                fprintf('x = %0.4e \t f(x1) = %0.4e \t err = %0.4e\n',x2,f2,err);
            end
            return

        end

    end
end

end

%--------------------------------------------------------------------------
% Auxiliary arrays for interpolation

function [Aux_arrays] = get_aux_arrays()

% This function returns a structure array containig the vectors needed for
% interpolation. There is no input, while the output is Aux_arrays with the
% following fields:
% Aux_arrays.Tc: Vector of length 100 containing the temperature values for
%                atmospheric air in Kelvin from 292 K (18.85°C) to 1700 K 
%                (1426.85°C) 

% Aux_arrays.mu: Vector of length 100 containing the values for the dynamic 
%                viscosity of the cooling air in correspondance with values 
%                of temperatures stored in Aux_arrays.Tc

% Aux_arrays.k:  Vector of length 100 containing the values for the thermal 
%                conductivity [W/(m*K)] of the cooling air in correspondance 
%                with values of temperatures stored in Aux_arrays.Tc

% Aux_arrays.Pr: Vector of length 100 containing the values for the Prandtl 
%                number of the cooling air in correspondance with values 
%                of temperatures stored in Aux_arrays.Tc

% Aux_arrays.cp: Vector of length 100 containing the values for the specific 
%                heat [J/K] of the cooling air in correspondance with values 
%                of temperatures stored in Aux_arrays.Tc

% Aux_arrays.gamma: Vector of length 100 containing the values for the isentropic 
%                expansion factor of the cooling air in correspondance with values 
%                of temperatures stored in Aux_arrays.Tc
% Aux_arrays.R: Ideal gas constant [J/(mol K)]

Aux_arrays.Tc = [292;306.222222222222;320.444444444444;334.666666666667;...
    348.888888888889;363.111111111111;377.333333333333;391.555555555556;...
    405.777777777778;420;434.222222222222;448.444444444444;462.666666666667;...
    476.888888888889;491.111111111111;505.333333333333;519.555555555556;...
    533.777777777778;548;562.222222222222;576.444444444445;590.666666666667;...
    604.888888888889;619.111111111111;633.333333333333;647.555555555556;...
    661.777777777778;676;690.222222222222;704.444444444445;718.666666666667;...
    732.888888888889;747.111111111111;761.333333333333;775.555555555556;...
    789.777777777778;804;818.222222222222;832.444444444445;846.666666666667;...
    860.888888888889;875.111111111111;889.333333333333;903.555555555556;...
    917.777777777778;932;946.222222222222;960.444444444445;974.666666666667;...
    988.888888888889;1003.11111111111;1017.33333333333;1031.55555555556;...
    1045.77777777778;1060;1074.22222222222;1088.44444444444;1102.66666666667;...
    1116.88888888889;1131.11111111111;1145.33333333333;1159.55555555556;...
    1173.77777777778;1188;1202.22222222222;1216.44444444444;1230.66666666667;...
    1244.88888888889;1259.11111111111;1273.33333333333;1287.55555555556;...
    1301.77777777778;1316;1330.22222222222;1344.44444444444;1358.66666666667;...
    1372.88888888889;1387.11111111111;1401.33333333333;1415.55555555556;...
    1429.77777777778;1444;1458.22222222222;1472.44444444444;1486.66666666667;...
    1500.88888888889;1515.11111111111;1529.33333333333;1543.55555555556;...
    1557.77777777778;1572;1586.22222222222;1600.44444444444;1614.66666666667;...
    1628.88888888889;1643.11111111111;1657.33333333333;1671.55555555556;...
    1685.77777777778;1700];

Aux_arrays.mu = [1.80789760363643e-05;1.87529738750714e-05;1.94118151297565e-05;...
    2.00563150105598e-05;2.06872276128228e-05;2.13052511093157e-05;...
    2.19110325706141e-05;2.25051723979384e-05;2.30882283738025e-05;...
    2.36607193478244e-05;2.42231285813986e-05;2.47759067777216e-05;...
    2.53194748242925e-05;2.58542262743400e-05;2.63805295922255e-05;...
    2.68987301860918e-05;2.74091522490983e-05;2.79121004286347e-05;...
    2.84078613410280e-05;2.88967049474897e-05;2.93788858054214e-05;...
    2.98546442077049e-05;3.03242072212604e-05;3.07877896349454e-05;...
    3.12455948257844e-05;3.16978155515523e-05;3.21446346768783e-05;...
    3.25862258392676e-05;3.30227540607636e-05;3.34543763103682e-05;...
    3.38812420218030e-05;3.43034935707190e-05;3.47212667150409e-05;...
    3.51346910017562e-05;3.55438901431251e-05;3.59489823649925e-05;...
    3.63500807296152e-05;3.67472934351847e-05;3.71407240940140e-05;...
    3.75304719911676e-05;3.79166323251478e-05;3.82992964320972e-05;...
    3.86785519948450e-05;3.90544832379992e-05;3.94271711101827e-05;...
    3.97966934544102e-05;4.01631251675138e-05;4.05265383494498e-05;...
    4.08870024432416e-05;4.12445843662546e-05;4.15993486334351e-05;...
    4.19513574730960e-05;4.23006709357821e-05;4.26473469967033e-05;...
    4.29914416521876e-05;4.33330090105655e-05;4.36721013778677e-05;...
    4.40087693386864e-05;4.43430618325251e-05;4.46750262259337e-05;...
    4.50047083807068e-05;4.53321527183995e-05;4.56574022813964e-05;...
    4.59804987907537e-05;4.63014827010169e-05;4.66203932522010e-05;...
    4.69372685191094e-05;4.72521454581527e-05;4.75650599518191e-05;...
    4.78760468509350e-05;4.81851400148490e-05;4.84923723496588e-05;...
    4.87977758445955e-05;4.91013816066711e-05;4.94032198936885e-05;...
    4.97033201457048e-05;5.00017110150364e-05;5.02984203948849e-05;...
    5.05934754466608e-05;5.08869026260742e-05;5.11787277080605e-05;...
    5.14689758106022e-05;5.17576714175056e-05;5.20448384001876e-05;...
    5.23305000385228e-05;5.26146790408006e-05;5.28973975628377e-05;...
    5.31786772262873e-05;5.34585391361875e-05;5.37370038977849e-05;...
    5.40140916326706e-05;5.42898219942611e-05;5.45642141826567e-05;...
    5.48372869589077e-05;5.51090586587147e-05;5.53795472055933e-05;...
    5.56487701235245e-05;5.59167445491177e-05;5.61834872433077e-05;...
    5.64490146026068e-05];

Aux_arrays.k = [0.0256170150115419;0.0267289925119716;0.0278215391967157;...
    0.0288953043864342;0.0299509270685847;0.0309890314187982;0.0320102237750685;...
    0.0330150906865435;0.0340041977563392;0.0349780890692354;0.0359372870477813;...
    0.0368822926193939;0.0378135856061693;0.0387316252709676;0.0396368509697758;...
    0.0405296828727826;0.0414105227260187;0.0422797546325820;0.0431377458379144;...
    0.0439848475077600;0.0448213954906015;0.0456477110587889;0.0464641016244153;...
    0.0472708614273918;0.0480682721942244;0.0488566037667968;0.0496361147010466;...
    0.0504070528358621;0.0511696558328414;0.0519241516877772;0.0526707592148839;...
    0.0534096885048821;0.0541411413581110;0.0548653116938668;0.0555823859371660;...
    0.0562925433841191;0.0569959565470718;0.0576927914806368;0.0583832080896937;...
    0.0590673604203925;0.0597453969351455;0.0604174607725437;0.0610836899930845;...
    0.0617442178115505;0.0623991728168296;0.0630486791799218;0.0636928568508351;...
    0.0643318217450288;0.0649656859200245;0.0655945577427647;0.0662185420482666;...
    0.0668377402900798;0.0674522506830305;0.0680621683386986;0.0686675853940520;...
    0.0692685911336315;0.0698652721056560;0.0704577122323968;0.0710459929151448;...
    0.0716301931340763;0.0722103895433033;0.0727866565613767;0.0733590664574945;...
    0.0739276894336504;0.0744925937029459;0.0750538455642722;0.0756115094735606;...
    0.0761656481117818;0.0767163224498692;0.0772635918107285;0.0778075139284849;...
    0.0783481450051133;0.0788855397645863;0.0794197515046673;0.0799508321464683;...
    0.0804788322818871;0.0810038012190284;0.0815257870257111;0.0820448365711557;...
    0.0825609955659414;0.0830743086003182;0.0835848191809524;0.0840925697661823;...
    0.0845976017998544;0.0850999557438083;0.0855996711090728;0.0860967864858351;...
    0.0865913395722386;0.0870833672020650;0.0875729053713497;0.0880599892639815;...
    0.0885446532763303;0.0890269310409486;0.0895068554493858;0.0899844586741562;...
    0.0904597721898979;0.0909328267937569;0.0914036526250308;0.0918722791841037;...
    0.0923387353507035];

Aux_arrays.Pr = [0.709693165548335;0.705826839396829;0.702323479788794;...
    0.699167273239829;0.696343958906374;0.693840076500763;0.691642434246302;...
    0.689737762617761;0.688112518164139;0.686752803144437;0.685644370053239;...
    0.684772684497571;0.684123024626764;0.683680599935976;0.683430676465925;...
    0.683358699053129;0.683450404295044;0.683691920298402;0.684069851131481;...
    0.684571345276265;0.685184148354364;0.685896641058923;0.686697863633507;...
    0.687577528458937;0.688526022390734;0.689534400473962;0.690594372581141;...
    0.691698284396996;0.692839094029568;0.694010345374028;0.695206139202586;...
    0.696421102807545;0.697650358888685;0.698889494252876;0.700134528784162;...
    0.701381885046515;0.702628358798538;0.703871090628778;0.705107538860859;...
    0.706335453828223;0.707552853577687;0.708758001028214;0.709949382585001;...
    0.711125688188479;0.712285792762007;0.713428739010215;0.714553721511514;...
    0.715660072042418;0.716747246067814;0.717814810329434;0.718862431464376;...
    0.719889865586172;0.720896948762417;0.721883588325031;0.722849754951871;...
    0.723795475461185;0.724720826263549;0.725625927419036;0.726510937250585;...
    0.727376047467727;0.728221478757927;0.729047476805793;0.729854308703362;...
    0.730642259717373;0.731411630382089;0.732162733888717;0.732895893744752;...
    0.733611441678786;0.734309715768322;0.734991058770024;0.735655816633568;...
    0.736304337181901;0.736936968942167;0.737554060112959;0.738155957654824;...
    0.738743006492072;0.739315548815064;0.739873923473079;0.740418465448781;...
    0.740949505406118;0.741467369304245;0.741972378070715;0.742464847327859;...
    0.742945087166770;0.743413401963919;0.743870090235800;0.744315444527523;...
    0.744749751331603;0.745173291033581;0.745586337881424;0.745989159975948;...
    0.746382019279789;0.746765171642649;0.747138866840831;0.747503348629206;...
    0.747858854803994;0.748205617274868;0.748543862145060;0.748873809798276;...
    0.749195674991353];

Aux_arrays.cp = [1005.60012020993;1006.02925331568;1006.58903307631;...
    1007.29626387805;1008.16540128977;1009.20809710156;1010.43294154241;...
    1011.84538292384;1013.44779544955;1015.23966215093;1017.21784023473;...
    1019.37687907824;1021.70936547451;1024.20627563052;1026.85731825628;...
    1029.65125748952;1032.57620819329;1035.61989927639;1038.76990313492;...
    1042.01383115740;1045.33949655826;1048.73504669572;1052.18906757655;...
    1055.69066353108;1059.22951512425;1062.79591830997;1066.38080768083;...
    1069.97576644974;1073.57302555070;1077.16545398335;1080.74654226327;...
    1084.31038058838;1087.85163309629;1091.36550937181;1094.84773417082;...
    1098.29451615478;1101.70251628065;1105.06881636113;1108.39088819919;...
    1111.66656360708;1114.89400554059;1118.07168051426;1121.19833240871;...
    1124.27295773758;1127.29478240578;1130.26323996288;1133.17795133325;...
    1136.03870598764;1138.84544450827;1141.59824249046;1144.29729571735;...
    1146.94290654042;1149.53547139654;1152.07546939133;1154.56345187953;...
    1157.00003297414;1159.38588091844;1161.72171025727;1164.00827474707;...
    1166.24636094717;1168.43678243793;1170.58037461485;1172.67799001072;...
    1174.73049410130;1176.73876155300;1178.70367287399;1180.62611143326;...
    1182.50696081447;1184.34710247432;1186.14741367756;1187.90876568270;...
    1189.63202215507;1191.31803778552;1192.96765709483;1194.58171340601;...
    1196.16102796767;1197.70640921354;1199.21865214435;1200.69853781947;...
    1202.14683294707;1203.56428956216;1204.95164478340;1206.30962063990;...
    1207.63892396034;1208.94024631738;1210.21426402095;1211.46163815470;...
    1212.68301465041;1213.87902439559;1215.05028337005;1216.19739280767;...
    1217.32093937980;1218.42149539732;1219.49961902842;1220.55585452976;...
    1221.59073248859;1222.60477007398;1223.59847129522;1224.57232726589;...
    1225.52681647207];

Aux_arrays.gamma = [1.39950058865879;1.39926213833470;1.39895152092704;...
    1.39855977444793;1.39807939772435;1.39750462009510;1.39683155099059;...
    1.39605821895242;1.39518451635508;1.39421206948837;1.39314405446856;...
    1.39198497840356;1.39074044300536;1.38941690496999;1.38802144435143;...
    1.38656154914530;1.38504492156625;1.38347930916266;1.38187236200845;...
    1.38023151574301;1.37856389916267;1.37687626435023;1.37517493690227;...
    1.37346578361797;1.37175419499100;1.37004507994838;1.36834287046500;...
    1.36665153391742;1.36497459129787;1.36331513967219;1.36167587751920;...
    1.36005913182492;1.35846688601916;1.35690080803075;1.35536227790102;...
    1.35385241453459;1.35237210128240;1.35092201014769;1.34950262448322;...
    1.34811426010946;1.34675708483173;1.34543113637114;1.34413633875139;...
    1.34287251720361;1.34163941166452;1.34043668895207;1.33926395370736;...
    1.33812075819327;1.33700661103982;1.33592098502395;1.33486332396812;...
    1.33383304883780;1.33282956311333;1.33185225750641;1.33090051408676;...
    1.32997370987882;1.32907121998413;1.32819242027952;1.32733668973739;...
    1.32650341240954;1.32569197911253;1.32490178884846;1.32413224999193;...
    1.32338278127072;1.32265281256495;1.32194178554665;1.32124915417983;...
    1.32057438509829;1.31991695787726;1.31927636521249;1.31865211301947;...
    1.31804372046361;1.31745071993117;1.31687265694967;1.31630909006515;...
    1.31575959068328;1.31522374287999;1.31470114318697;1.31419140035654;...
    1.31369413510984;1.31320897987196;1.31273557849694;1.31227358598528;...
    1.31182266819638;1.31138250155777;1.31095277277292;1.31053317852905;...
    1.31012342520628;1.30972322858904;1.30933231358081;1.30895041392290;...
    1.30857727191774;1.30821263815750;1.30785627125815;1.30750793759950;...
    1.30716741107140;1.30683447282630;1.30650891103831;1.30619052066880;...
    1.30587910323874];

Aux_arrays.R = 8.31446261815324;

end

%--------------------------------------------------------------------------
