classdef FNSolver
    %FNSOLVER implements a solver for the FitzHugh-Nagumo equations
    
    %   Copyright (c) 2018, Politecnico di Milano 
    %   localROM - Stefano Pagani <stefano.pagani at polimi.it>
    
    properties
        L                % domain lenght
        epsilon          % sqrt of the conductivity constant
        b                % recovery coefficient
        gamma            % recovery coefficient
        Nh               % spatial discretization
        t0               % initial time
        tF               % final time
        Nt               % number of time-steps
        Xnorm            % metric for error compuation
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
            
            % H1-norm 
            
            % space discretization
            dx = 1/obj.Nh;
            x = [0:dx:1]';
            
            % diffusion matrix assembling            
            h = diff(x(1:2));
            A = sparse( diag(2./diff(x)) - diag(1./diff(x(1:end-1)),1) - diag(1./diff(x(1:end-1)),-1) ) ;
            A(1,1) = 1./h;
            A(end+1,end+1) = 1./h;  A(end-1,end) = -1./h;  A(end,end-1) = -1./h;

            % reaction matrix assembling            
            M_ref = sparse( diag( 4*diff(x)/6 ) + diag(1*diff(x(1:end-1))/6,1) + diag(1*diff(x(1:end-1))/6,-1) ) ;
            M_ref(1,1) = 2*h/6;
            M_ref(end+1,end+1) = 2*h/6;  M_ref(end-1,end) = h./6;  M_ref(end,end-1) = h./6;
            
            obj.Xnorm = A + M_ref;
            
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
            
            
            A_ref    = newEpsilon*newEpsilon.*A;            

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

            % space discretization lenght
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
                FFOM(1)       = FFOM(1) + 50000*(t(j+1))^3*exp(-15*t(j+1))*newEpsilon*newEpsilon;

                % update voltage
                u(:,j+1) = AFOM \ FFOM;
   
            end
        end

     function [u,w] = solveROM(obj,newEpsilon,LROMclass)
            %SOLVEFOM this method provides the reduced-order solution of the
            % FN problem
            
            V = LROMclass.V;
            
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
            
            % projection of full-order stiffness matrix in the global and the local
            % case
            if strcmp(LROMclass.clusterType,'Global')
                
                A_ROM = V'*(A*V);            
                A_ref_ROM    = newEpsilon*newEpsilon.*A_ROM;            

            else
                % loop over the clusters
                for iC = 1 : LROMclass.clusterNumber
                    A_ROM{iC} = V{iC}'*(A*V{iC});            
                    A_ref_ROM{iC}    = V{iC}'*((newEpsilon*newEpsilon.*A)*V{iC}) ;
                end
            end
            
            % reaction matrix assembling            
            M_ref = sparse( diag( 4*diff(x)/6 ) + diag(1*diff(x(1:end-1))/6,1) + diag(1*diff(x(1:end-1))/6,-1) ) ;
            M_ref(1,1) = 2*h/6;
            M_ref(end+1,end+1) = 2*h/6;  M_ref(end-1,end) = h./6;  M_ref(end,end-1) = h./6;
            
            % projection of full-order mass matrix + initial data definition in the global and the local
            % case
            if strcmp(LROMclass.clusterType,'Global')
            
                M_ref_ROM = V'*M_ref*V ;
                M_ref_halfROM = V'*M_ref ;                
                M_ref_w = V'*M_ref*LROMclass.Vw ; 
                M_ref_dt_ROM = newEpsilon/dt*M_ref_ROM ;
                
                % intial data
                u = zeros(obj.Nh+1,obj.Nt+1);
                uROM = V'*u;
                w = u; 
                wROM = LROMclass.Vw'*w;
                OnesP = LROMclass.Vw'*ones(size(w(:,1)));
                

            else
                
                for iC = 1:LROMclass.clusterNumber
                    
                    M_ref_ROM{iC} = V{iC}'*(M_ref*V{iC}) ;
                    M_ref_halfROM{iC} = V{iC}'*M_ref ;            
                    M_ref_dt_ROM{iC} = V{iC}'*((newEpsilon/dt*M_ref)*V{iC});                    
                    M_ref_dt_halfROM{iC} = V{iC}'*( newEpsilon/dt*M_ref )  ;  
                    M_ref_w{iC} = V{iC}'*(M_ref*LROMclass.Vw{iC}) ;

                end
                
                % intial data
                u = zeros(obj.Nh+1,obj.Nt+1);
                w = u;
                              

            end
            
            % ionic term
            f = @(v)  v.*(v-0.1).*(v-1); 

            % space discretization lenght
            N = obj.Nh;
            
            c = 0;
            
            for j=1:obj.Nt

                % global solution update
                if strcmp(LROMclass.clusterType,'Global')
                    
                    % update recovery variable
                    wROM(:,j+1) = 1/(1+dt*2)*( dt*c*OnesP + wROM(:,j) + dt*obj.b*(LROMclass.Vw'*u(:,j)) );
                    
                    % left-hand side
                    AROM = M_ref_dt_ROM + A_ref_ROM;
                    
                    % right-hand side
                    FROM =  M_ref_dt_ROM*(uROM(:,j) ) ... %- A_ref*u(:,j+1)
                         - M_ref_halfROM*f(u(:,j))  ...
                        -M_ref_w*wROM(:,j+1);
                    
                    % external stimulus 
                    FROM       = FROM  +  V(1,:)'* ( 50000*(t(j+1))^3*exp(-15*t(j+1))*newEpsilon*newEpsilon );

                    % update voltage
                    uROM(:,j+1)  = AROM \ FROM;
   
                    % storaging of the solution
                    u(:,j+1)    = V*uROM(:,j+1) ;
                    w(:,j+1)    = LROMclass.Vw*wROM(:,j+1) ;
                    
                else
                    
                    % cluster selection: case k-means State                    
                    if strcmp(LROMclass.clusterType,'kmeansState')
                        for iC = 1 : LROMclass.clusterNumber
                            distv(iC) = pdist( [ u(:,j) , LROMclass.centroids{iC} ]' , 'squaredeuclidean' );
                        end
                        [val,iSel] = min(distv);
                    end
                    % cluster selection: case PEBL
                    if strcmp(LROMclass.clusterType,'PEBL')
                        
                        % distance wrt node 1
                        EP1 = pdist( [ u(:,j) , LROMclass.centroids{1}*LROMclass.centroids{1}'*u(:,j) ]' , 'squaredeuclidean'  );
                        % distance wrt node 2
                        EP2 = pdist( [ u(:,j) , LROMclass.centroids{2}*LROMclass.centroids{2}'*u(:,j) ]' , 'squaredeuclidean'  );

                        if EP1>=EP2
                            iSel=2;
                        else
                            iSel=1;
                        end
                         
                        % binary tree searching 
                        for iTree = 2 : length(LROMclass.tree) 
                            
                            if iSel == LROMclass.tree{iTree}{1}
                                
                                i2 = LROMclass.tree{iTree}{2};
                                
                                % distance wrt node iSel
                                EP1 = pdist( [ u(:,j) , LROMclass.centroids{iSel}*LROMclass.centroids{iSel}'*u(:,j) ]' , 'squaredeuclidean'  );
                                % distance wrt node i2    
                                EP2 = pdist( [ u(:,j) , LROMclass.centroids{i2}*LROMclass.centroids{i2}'*u(:,j) ]' , 'squaredeuclidean'  );

                                if EP1>=EP2
                                    iSel=i2;
                                end
                            end
                                
                        end
                        
                    end
                    
                    % cluster selection: case k-means parameters
                    if strcmp(LROMclass.clusterType,'kmeansParam') && j==1
                        for iC = 1 : LROMclass.clusterNumber
                            distv(iC) = pdist( [ (newEpsilon) , LROMclass.centroids{iC} ]' , 'squaredeuclidean' );
                        end
                        [val,iSel] = min(distv);
                    end
                    
                    % cluster selection: case Time-based subdivision
                    if strcmp(LROMclass.clusterType,'Time') 
                        for iC = 1 : LROMclass.clusterNumber
                            distv(iC) = ( j>=LROMclass.centroids{iC}(1) && j<=LROMclass.centroids{iC}(2) );
                        end
                        [val,iSel] = max(distv);
                    end
                    
                    % initial data
                    if j==1
                        
                        uROM = zeros( size(M_ref_ROM{iSel}(:,1)) );
                        wROM = zeros( size( LROMclass.Vw{iSel}(1,:)' ) );
                        
                        OnesP = LROMclass.Vw{iSel}' * ones(size(w(:,1)));
                        
                    else
                        % projection on the current selected subspaces
                        if iSelOld~=iSel

                            wROM = LROMclass.Vw{iSel}' * ( w(:,j) );
                           
                            
                            OnesP = LROMclass.Vw{iSel}' * ones(size(w(:,j)));
                            
                        end
                    end

                    % solution updating
                    
                    % update recovery variable
                    wROM = 1/(1+dt*2)*( dt*c*OnesP + wROM + dt*obj.b*(LROMclass.Vw{iSel}'*u(:,j)) );
                    
                    % left-hand side
                    AROM = M_ref_dt_ROM{iSel} + A_ref_ROM{iSel};

                    % right-hand side
                    FROM =  M_ref_dt_halfROM{iSel}*(u(:,j) ) ... 
                         - M_ref_halfROM{iSel}*f(u(:,j))  ...
                         -M_ref_w{iSel}*wROM;
                   
                    % external stimulus 
                    FROM  = FROM  +  V{iSel}(1,:)'*( 50000*(t(j+1))^3*exp(-15*t(j+1))*newEpsilon*newEpsilon );  
                    
                    % update voltage
                    uROM = AROM \ FROM;
   
                    % storaging of the solution
                    u(:,j+1)    = V{iSel}*uROM;
                    w(:,j+1)    = LROMclass.Vw{iSel}*wROM;
                    
                    iSelOld = iSel;
                    
                end
           
            end
     end
     
     function [u,w] = solveROMHyperRed(obj,newEpsilon,LROMclass)
            %SOLVEFOM this method provides the Hyperreduced solution of the
            % FN problem
            
            V = LROMclass.V;
            
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
            
            % projection of full-order stiffness matrix in the global and the local
            % case
            if strcmp(LROMclass.clusterType,'Global')
                
                A_ROM       = V'*(A*V);            
                A_ref_ROM   = newEpsilon*newEpsilon.*A_ROM;            

            else
                % loop over the clusters
                for iC = 1 : LROMclass.clusterNumber
                    A_ROM{iC}        = V{iC}'*(A*V{iC});            
                    A_ref_ROM{iC}    = V{iC}'*((newEpsilon*newEpsilon.*A)*V{iC}) ;
                end
            end
            
            % mass matrix assembling            
            M_ref = sparse( diag( 4*diff(x)/6 ) + diag(1*diff(x(1:end-1))/6,1) + diag(1*diff(x(1:end-1))/6,-1) ) ;
            M_ref(1,1) = 2*h/6;
            M_ref(end+1,end+1) = 2*h/6;  M_ref(end-1,end) = h./6;  M_ref(end,end-1) = h./6;
            
            % projection of full-order mass matrix + initial data definition in the global and the local
            % case
            if strcmp(LROMclass.clusterType,'Global')
            
                M_ref_ROM = V'*M_ref*V ;
                %M_ref_halfROM = V'*M_ref ;                          
                M_ref_w = V'*M_ref*LROMclass.Vw ; 
                M_ref_dt_ROM = newEpsilon/dt*M_ref_ROM ;
                
                % intial data
                u = zeros(obj.Nh+1,obj.Nt+1);
                uROM = V'*u;
                w = u;             
                wROM = LROMclass.Vw'*w;
                OnesP = LROMclass.Vw'*ones(size(w(:,1)));
                
                % HyperRed structures             
                MD = (LROMclass.PMat'*LROMclass.VD);
                PDEIM = (LROMclass.VD) * ((MD)\eye(size(MD)));                
                M_ref_DEIM = V'*(M_ref * PDEIM) ;
                

            else
                %loop over the clusters                
                for iC = 1:LROMclass.clusterNumber
                    M_ref_ROM{iC} = V{iC}'*(M_ref*V{iC}) ;
                    M_ref_halfROM{iC} = V{iC}'*M_ref ;            
                    M_ref_dt_ROM{iC} = V{iC}'*((newEpsilon/dt*M_ref)*V{iC});                    
                    M_ref_dt_halfROM{iC} = V{iC}'*( newEpsilon/dt*M_ref )  ;  
                    M_ref_w{iC} = V{iC}'*(M_ref*LROMclass.Vw{iC}) ;
                    
                    % HyperReduction structures
                    MD = (LROMclass.PMat{iC}'*LROMclass.VD{iC});
                    PDEIM = (LROMclass.VD{iC}) * ((MD)\eye(size(MD)));
                    M_ref_DEIM{iC} = V{iC}'*(M_ref * PDEIM) ;
                
                end
                
                % intial data
                u = zeros(obj.Nh+1,obj.Nt+1);
                w = u;
         
            end
        
            % ionic term
            f = @(v)  v.*(v-0.1).*(v-1); 

            % space discretization dimension
            N = obj.Nh;
            
            c = 0;
            
            % loop over time
            for j=1:obj.Nt
                
                % solution update: global procedure
                if strcmp(LROMclass.clusterType,'Global')
                    
                    % update recovery variable
                    wROM(:,j+1) = 1/(1+dt*2)*( dt*c*OnesP + wROM(:,j) + dt*obj.b*(LROMclass.Vw'*u(:,j)) );
                    
                    % left-hand side
                    AROM = M_ref_dt_ROM + A_ref_ROM;
                    
                    % right-hand side
                    FROM =  M_ref_dt_ROM*(uROM(:,j) ) ... 
                        - M_ref_DEIM * f(u(LROMclass.iDEIM,j))   ... 
                        -M_ref_w*wROM(:,j+1);
                    
                    % external stimulus 
                    FROM       = FROM  +  V(1,:)'* ( 50000*(t(j+1))^3*exp(-15*t(j+1))*newEpsilon*newEpsilon );

                    % update voltage
                    uROM(:,j+1)  = AROM \ FROM;
   
                    u(:,j+1)    = V*uROM(:,j+1) ;
                    w(:,j+1)    = LROMclass.Vw*wROM(:,j+1) ;
                    
                else
                    
                    % cluster selection: case k-means state
                    if strcmp(LROMclass.clusterType,'kmeansState')
                        for iC = 1 : LROMclass.clusterNumber
                            distv(iC) = pdist( [ u(:,j) , LROMclass.centroids{iC} ]' , 'squaredeuclidean' );
                        end
                        [val,iSel] = min(distv);
                    end
                    
                    % cluster selection: case PEBL
                    if strcmp(LROMclass.clusterType,'PEBL')
                       
                        % distance wrt node 1
                        EP1 = pdist( [ u(:,j) , LROMclass.centroids{1}*LROMclass.centroids{1}'*u(:,j) ]' , 'squaredeuclidean'  );
                        % distance wrt node 2
                        EP2 = pdist( [ u(:,j) , LROMclass.centroids{2}*LROMclass.centroids{2}'*u(:,j) ]' , 'squaredeuclidean'  );

                        if EP1>=EP2
                            iSel=2;
                        else
                            iSel=1;
                        end
                         
                        % binary tree searching
                        for iTree = 2 : length(LROMclass.tree) 
                            
                            if iSel == LROMclass.tree{iTree}{1}
                                
                                i2 = LROMclass.tree{iTree}{2};
                                
                                % distance wrt node iSel
                                EP1 = pdist( [ u(:,j) , LROMclass.centroids{iSel}*LROMclass.centroids{iSel}'*u(:,j) ]' , 'squaredeuclidean'  );
                                % distance wrt node i2
                                EP2 = pdist( [ u(:,j) , LROMclass.centroids{i2}*LROMclass.centroids{i2}'*u(:,j) ]' , 'squaredeuclidean'  );

                                if EP1>=EP2
                                    iSel=i2;
                                end
                            end
                                
                        end
                        
                    end
                    
                    % cluster selection: case k-means parameters
                    if strcmp(LROMclass.clusterType,'kmeansParam') && j==1
                        for iC = 1 : LROMclass.clusterNumber
                            distv(iC) = pdist( [ (newEpsilon) , LROMclass.centroids{iC} ]' , 'squaredeuclidean' );
                        end
                        [val,iSel] = min(distv);
                    end
                    
                    % cluster selection: case time-based subdivision
                    if strcmp(LROMclass.clusterType,'Time') 
                        for iC = 1 : LROMclass.clusterNumber
                            distv(iC) = ( j>=LROMclass.centroids{iC}(1) && j<=LROMclass.centroids{iC}(2) );
                        end
                        [val,iSel] = max(distv);
                    end

                    % initial data
                    if j==1
                        uROM = zeros( size(M_ref_ROM{iSel}(:,1)) );
                        wROM = zeros( size( LROMclass.Vw{iSel}(1,:)' ) );
                        
                        OnesP = LROMclass.Vw{iSel}' * ones(size(w(:,1)));
                        
                    else
                        % projection over the current selected subspace
                        if iSelOld~=iSel

                            wROM = LROMclass.Vw{iSel}' * ( w(:,j) );
                            
                            OnesP = LROMclass.Vw{iSel}' * ones(size(w(:,j)));
                            
                        end
                    end

                    % update recovery variable
                    wROM = 1/(1+dt*2)*( dt*c*OnesP + wROM + dt*obj.b*(LROMclass.Vw{iSel}'*u(:,j)) );
                    
                    % left-hand side
                    AROM = M_ref_dt_ROM{iSel} + A_ref_ROM{iSel};

                    % right-hand side
                    FROM =  M_ref_dt_halfROM{iSel}*(u(:,j) ) ... 
                          - M_ref_DEIM{iSel} * f(u(LROMclass.iDEIM{iSel},j))  ... 
                         -M_ref_w{iSel}*wROM;
                   
                    % external stimulus 
                    FROM       = FROM  +  V{iSel}(1,:)'*( 50000*(t(j+1))^3*exp(-15*t(j+1))*newEpsilon*newEpsilon );  
                    
                    % update voltage
                    uROM = AROM \ FROM;
   
                    % solution storaging
                    u(:,j+1)    = V{iSel}*uROM;
                    w(:,j+1)    = LROMclass.Vw{iSel}*wROM;
                    
                    iSelOld = iSel;
                    
                end
                        
            end
     end
     
    end
    
end
     
