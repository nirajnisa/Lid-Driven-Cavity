close all
clc;
%% Geometry and property values
L = 1;
N = 101;
h = L/(N-1);
Re = 100; %Change Reynolds number here
mu = 1/Re;
rho = 1;
%% Initialization
u_final(N,N)=0;
v_final(N,N)=0;
p_final(N,N)=1;

u_final(1,:)=1;

alpha_p =0.8;
%staggerd grid arrangement

u(N+1,N) = 0;
de(N+1,N) = 0;
v(N,N+1) = 0;
dn(N,N+1) = 0;
p(N+1,N+1) = 1;
u(1,:) = 2;

u_new(N+1,N) = 0;
v_new(N,N+1) = 0;
p_new(N+1,N+1) = 1;
u_new(1,:) = 2;

pc(N+1,N+1) = 0;
b(N+1,N+1) = 0;

%% SIMPLE algorithm
error = 1;
iter = 0;
figure(1);
while error>1e-6
    % x-momentum interior
    for i = 2:N
        for j = 2:N-1
            ue = 0.5*(u(i,j+1)+u(i,j));
            uw = 0.5*(u(i,j-1)+u(i,j));
            vn = 0.5*(v(i-1,j)+v(i-1,j+1));
            vs = 0.5*(v(i,j)+v(i,j+1));
            
            Fe = rho*h*ue;
            Fw = rho*h*uw;
            Fn = rho*h*vn;
            Fs = rho*h*vs;
            De = mu;
            Dw = mu;
            Dn = mu;
            Ds = mu;
            
            aE = De + max(-Fe,0);
            aW = Dw + max(Fw,0);
            aN = Dn + max(-Fn,0);
            aS = Ds + max(Fs,0);
            
            ae = aE + aW + aN + aS + (Fe - Fw + Fn - Fs);
            Ae = -h;
            de(i,j) = Ae/ae;
            
            u(i,j) = (aE*u(i,j+1) + aW*u(i,j-1) + aN*u(i-1,j) + aS*u(i+1,j))/ae + de(i,j)*(p(i,j+1)-p(i,j));
        end
    end
    
    %x-momentum boundary
    u(1,:) = 2 - u(2,:); %top boundary
    u(N+1,:) = - u(N,:); %bottom boundary
    u(2:N,1) = 0; %left boundary
    u(2:N,N) = 0; %right boundary
 
    
    %y-momentum interior
    for i = 2:N-1
        for j = 2:N
            ue = 0.5*(u(i,j)+u(i+1,j));
            uw = 0.5*(u(i,j-1)+u(i+1,j-1));
            vn = 0.5*(v(i-1,j)+v(i,j));
            vs = 0.5*(v(i,j)+v(i+1,j));
            
            Fe = rho*h*ue;
            Fw = rho*h*uw;
            Fn = rho*h*vn;
            Fs = rho*h*vs;
            De = mu;
            Dw = mu;
            Dn = mu;
            Ds = mu;
            
            aE = De + max(-Fe,0);
            aW = Dw + max(Fw,0);
            aN = Dn + max(-Fn,0);
            aS = Ds + max(Fs,0);
            
            an = aE + aW + aN + aS + (Fe - Fw + Fn - Fs);
            An = -h;
            dn(i,j) = An/an;
            
            v(i,j) = (aE*v(i,j+1) + aW*v(i,j-1) + aN*v(i-1,j) + aS*v(i+1,j))/an + dn(i,j)*(p(i,j)-p(i+1,j));
        end
    end
    
    %y-momentum boundary
    v(:,1) = -v(:,2); %left boundary
    v(:,N+1) = -v(:,N); %right boundary
    v(1,2:N) = 0; %top boundary
    v(N,2:N) = 0; %bottom boundary
    
    %continuity aka pressure correction
    pc(1:N+1,1:N+1)=0;
    
    for i=2:N
        for j=2:N
            Fe = rho*h*u(i,j);
            Fw = rho*h*u(i,j-1);
            Fn = rho*h*v(i-1,j);
            Fs = rho*h*v(i,j);
            
            aE = -rho*h*de(i,j);
            aW = -rho*h*de(i,j-1);
            aN = -rho*h*dn(i-1,j);
            aS = -rho*h*dn(i,j);
            aP = aE + aW + aN + aS;
            b(i,j)= -(Fe - Fw + Fn - Fs);
            
            pc(i,j) = (aE*pc(i,j+1) + aW*pc(i,j-1) + aN*pc(i-1,j) + aS*pc(i+1,j) + b(i,j))/aP;
            %coreecting pressure
            p_new(i,j) = p(i,j) + alpha_p*pc(i,j);
        end
    end
    
    %pressure correction at boundary
    p_new(1,:) = p_new(2,:);
    p_new(N + 1,:) = p_new(N,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,N + 1) = p_new(:,N);
    
    %correction velocities
    %u
    for i = 2:N
        for j = 2:N-1
            u_new(i,j) = u(i,j) + de(i,j)*(pc(i,j+1)-pc(i,j));
        end
    end
    %x-momentum boundary
    u_new(1,:) = 2 - u_new(2,:); %top boundary
    u_new(N+1,:) = - u_new(N,:); %bottom boundary
    u_new(2:N,1) = 0; %left boundary
    u_new(2:N,N) = 0; %right boundary

    
    %v
    for i = 2:N-1
        for j = 2:N
            v_new(i,j) = v(i,j) + dn(i,j)*(pc(i,j)-pc(i+1,j));
        end
    end
    v_new(:,1) = -v_new(:,2); %left boundary
    v_new(:,N+1) = -v_new(:,N); %right boundary
    v_new(1,2:N) = 0; %top boundary
    v_new(N,2:N) = 0; %bottom boundary
    
    %monitoring residual
    error = 0;
    for i = 2:N
        for j = 2:N
            error = error + abs(b(i,j));
        end
    end
    if(rem(iter,500)==0)
       figure(1);
       semilogy(iter, error, 'ro')
       hold on
       xlabel('Iterations')
       ylabel('Residual Error')
    end
    % Finishing the iteration
    u = u_new;
    v = v_new;
    p = p_new;
    iter = iter + 1;
end

% collocated variables
for i = 1:N
    for j = 1:N
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end

%% Centreline u velocity
figure(2);
y=0:h:L;
plot(u_final(:,(N+1)/2),1-y, 'LineWidth', 1)

xlabel('u')
ylabel('y')
