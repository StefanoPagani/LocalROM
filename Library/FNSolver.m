classdef FNSolver
    %FNSOLVER implements a solver for the FitzHugh-Nagumo equations
    
    properties
        L  % domain lenght
        epsilon % sqrt of the conductivity constant
        b % recovery coefficient
        gamma % recovery coefficient
        Nh % spatial discretization
        t0 % initial time
        tF % final time
        Nt % number of time-steps
    end
    
    methods
        function obj = FNSolver(param, Nh, t0, tF, Nt)
            %FNSOLVER Construct an instance of this class 
            %  inizialization of the model parameters:
            obj.L       = param(1);
            obj.epsilon = param(2);
            obj.b       = param(3);
            obj.gamma   = param(4);
            obj.Nh      = Nh;
            obj.t0      = t0;
            obj.tF      = tF;
            obj.Nt      = Nt;
        end
        
        function [u,w] = solveFOM(obj,newEpsilon)
            %SOLVEFOM this method provides the Full-order solution of the
            % FN problem
            
            % define the spatial and temporal mesh
            dx = 1/obj.Nh;
            x = [0:dx:1]';
            
            dt = (obj.tF-obj.t0)/obj.Nt;
            t = obj.t0:dt:obj.tF;
            
            % diffusion matrix assembling
            
            h = diff(x(1:2));
            A = sparse( diag(2./diff(x)) - diag(1./diff(x(1:end-1)),1) - diag(1./diff(x(1:end-1)),-1) ) ;
            A(1,1) = 1./h;
            A(end+1,end+1) = 1./h;  A(end-1,end) = -1./h;  A(end,end-1) = -1./h;
            
            
            A_ref    = newEpsilon.^2.*A;            

            % reaction matrix assembling            
            M_ref = sparse( diag( 4*diff(x)/6 ) + diag(1*diff(x(1:end-1))/6,1) + diag(1*diff(x(1:end-1))/6,-1) ) ;
            M_ref(1,1) = 2*h/6;
            M_ref(end+1,end+1) = 2*h/6;  M_ref(end-1,end) = h./6;  M_ref(end,end-1) = h./6;
            
            
            M_ref_dt = newEpsilon/dt*M_ref;
            
            
            
            % intial data
            u = zeros(obj.Nh+1,obj.Nt+1);
            w = u;

            
            % ionic term
            f = @(v)  v.*(v-0.1).*(v-1); 

            N = obj.Nh;
            
            c = 0;
            
            for j=1:obj.Nt

                % update recovery variable
                w(:,j+1) = 1/(1+dt*2)*( dt*c + w(:,j) + dt*obj.b*u(:,j) );
                    
               % left-hand side
               AFOM(1:obj.Nh+1,1:obj.Nh+1) = M_ref_dt + A_ref ; 

               % right-hand side
               FFOM(1:N+1,1) =  M_ref_dt*(u(:,j) ) ... %- A_ref*u(:,j+1)
                     - M_ref*f(u(:,j))  ...
                    -M_ref*w(:,j+1);

               % external stimulus 
               FFOM(1)       = FFOM(1) + 50000*(t(j+1))^3*exp(-15*t(j+1))*newEpsilon^2;

               % update voltage
               u(:,j+1) = AFOM \ FFOM;
   
            end
     end

    end    
end

