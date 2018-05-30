classdef localizedReduction
    %LOCALIZEDREDUCTION This class implements the localized reduction
    %strategies presented in the paper S. Pagani, A. Manzoni, A. Quarteroni. "Numerical approximation of
    %parametrized problems in cardiac electrophysiology by a local reduced
    %basis method", 2018
    
    %   Copyright (c) 2018, Politecnico di Milano 
%   localROM - Stefano Pagani <stefano.pagani at polimi.it>
    
    
    properties
        V             % struct containing the basis function related to the selected clusters
        centroids     % struct containing the centroids of the clusters   
        clusterType   % flag identifying the possible cluster technique
        clusterNumber % integer identifying the number of clusters
        tree          % PEBL structure
        VD            % struct containing non-linear term basis functions
        iDEIM         % struct containing DEIM indices
        PMat          % matrix
        Vw            % struct containing the basis function related to the selected clusters for the recovery variable
    end
    
    methods
        function obj = localizedReduction(inputType,inputNclust)
            %LOCALIZEDREDUCTION Construct an instance of this class
            % 
            obj.clusterType   = inputType;
            obj.clusterNumber = inputNclust; 
        end
        
        function obj = offline(obj,nTrain,PODtol,DEIMtol)
            %OFFLINE offline phase for the construction of the clusters and
            %the associated basis functions
           
            % model inputs
            param(1) = 1;      % domain lenght
            param(2) = 0.015;  % conducibility
            param(3) = 0.5;    % recovery parameter
            param(4) = 2;      % recovery parameter
            
            % number of time-steps
            Nt = 399;

            % solver constructor
            FNS = FNSolver(param, 1024, 0, 2, Nt );
            
            %Snapshots matrix
            S = [];
            Sw = [];
            
            % for reproducibility
            rng('default')
            
            % loop on the training set
            for iT = 1 : nTrain
                
                % equispaced points
                SParam(iT) = 0.003 + (iT-1)/(nTrain-1) *(0.05-0.003);
                
                % random realization
                %SParam(iT) = 0.003 + rand(1,1)*(0.05-0.003);
                
                % full-order solver
                [u,w] = FNS.solveFOM( SParam(iT) );
                
                % solution storaging
                S = [S , u];
                Sw = [Sw , w];
                
            end
            
            % non linear term snapshots
            f = @(v)  v.*(v-0.1).*(v-1); 
            SNL = f(S);
            
            % centroids computation: kmeans State
            if  strcmp(obj.clusterType,'kmeansState')
                
                % k-means algorithm
                [IDX, C] = kmeans(S', obj.clusterNumber, 'MaxIter', 1000 );
                
                % snapshots subdivision
                for iC = 1:obj.clusterNumber
                    
                    obj.centroids{iC} = C(iC,:)';
                    
                    % identify the snapshots related to the current cluster
                    selSnap = find(IDX==iC) ; 
                    SnapClust{iC} = S(: , selSnap + (mod(selSnap,Nt+1)~=0)  );                    
                    SnapClustW{iC} = Sw(: , selSnap + (mod(selSnap,Nt+1)~=0)  );
                    SnapClustNL{iC} = SNL(: , selSnap + (mod(selSnap,Nt+1)~=0)  );
                    
                    % compute POD
                    % voltage
                    [obj.V{iC}, ~ , sigma{iC}] = obj.POD(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    % recovery variable
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    % nonlinear terms for hyperreduction purposes
                    if nargin>3  
                    
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    
                end
                               
            end
                  
            % centroids computation: kmeans parameters
            if  strcmp(obj.clusterType,'kmeansParam')
                
                % k-means parameters
                [IDX, C] = kmeans((SParam'), obj.clusterNumber);
                
                % snapshots subdivision
                for iC = 1:obj.clusterNumber
                    obj.centroids{iC} = C(iC,:)';
                    
                    %  identify the snapshots related to the current cluster
                    selParam = find(IDX==iC) ; 
                    SnapClust{iC} = [];
                    SnapClustW{iC} = [];
                    SnapClustNL{iC} = [];
                    
                    for iTP = 1:length(selParam)

                        iT = selParam(iTP);
                        
                        SnapClust{iC} =  [ SnapClust{iC}  , S(:, 1 + (iT-1)*(Nt+1) : (iT)*(Nt+1) ) ];                        
                        SnapClustW{iC} =  [ SnapClustW{iC}  , Sw(:, 1 + (iT-1)*(Nt+1) : (iT)*(Nt+1) ) ];                       
                        SnapClustNL{iC} =  [ SnapClustNL{iC}  , SNL(:, 1 + (iT-1)*(Nt+1) : (iT)*(Nt+1) ) ];
                       
                    end
                    

                    % compute POD
                    % voltage
                    [obj.V{iC}, ~ , sigma{iC}] = obj.POD( uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    % recovery variable
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    % nonlinear terms for hyperreduction purposes
                    if nargin>3  
                    
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    
                end
                
            end

            % centroids computation: time-based
            if  strcmp(obj.clusterType,'Time')
                
                % snapshots subdivision
                for iC = 1:obj.clusterNumber
                    
                    obj.centroids{iC} = [1+(iC-1)*400/obj.clusterNumber,(iC)*400/obj.clusterNumber];
                    
                    SnapClust{iC} = [];
                    SnapClustW{iC} = [];
                    SnapClustNL{iC} = [];
                    
                    for iT = 1 : nTrain
                        
                        selSnap = [(iT-1)*(Nt+1)+obj.centroids{iC}(1) :  (iT-1)*(Nt+1)+obj.centroids{iC}(2) ];
                        SnapClust{iC} =  [ SnapClust{iC}  , S(:, selSnap  ) ];                    
                        SnapClustW{iC} =  [ SnapClustW{iC}  , Sw(:, selSnap  ) ];                       
                        SnapClustNL{iC} =  [ SnapClustNL{iC}  , SNL(:, selSnap  ) ];
                        
                    end

                    % compute POD
                    % voltage
                    [obj.V{iC}, ~, sigma{iC}] = obj.POD(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    %recovery variable
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    % nonlinear terms for hyperreduction purposes
                    if nargin>3  
                    
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    
                end
                                   
            end
            
            % centroids computation: PEBL based
            if  strcmp(obj.clusterType,'PEBL')
                
                % binary tree construction
                
                for iS = 1:size(S,2)
                    normvec(iS) = norm( S(:,iS) );
                end
                [val,ind] = max(normvec);
                
                obj.centroids{1} = S(:,ind)./norm(S(:,ind));
                
                SnapClust{1} = S;
                
                SnapClustW{1} = Sw;
                
                SnapClustNL{1} = SNL;
                
                % compute the distance wrt the first centroid
                for iS = 1:size(SnapClust{1},2)
                    distvec(iS) = pdist( [ SnapClust{1}(:,iS) , obj.centroids{1}*obj.centroids{1}'*SnapClust{1}(:,iS)  ]' , 'squaredeuclidean'  );
                end
                
                [val,ind] = max(distvec);
                
                obj.centroids{2} = SnapClust{1}(:,ind)./norm(SnapClust{1}(:,ind));
                
                indep1 = [];
                indep2 = [];
                
                
                for iS = 1:size(SnapClust{1},2)
                    
                    % distance wrt node 1
                    EP1 = pdist( [ SnapClust{1}(:,iS) , obj.centroids{1}*obj.centroids{1}'*SnapClust{1}(:,iS) ]' , 'squaredeuclidean'  );
                    % distance wrt node 2
                    EP2 = pdist( [ SnapClust{1}(:,iS) , obj.centroids{2}*obj.centroids{2}'*SnapClust{1}(:,iS) ]', 'squaredeuclidean'  );
                
                    if EP1>=EP2
                        indep2 = [indep2, iS];
                    else
                        indep1 = [indep1, iS];
                    end
                end
                
                % snapshots branching
                SnapClust{2} = SnapClust{1}(:,indep2 + (mod(indep2,Nt+1)~=0) );
                SnapClust{1} = SnapClust{1}(:,indep1 + (mod(indep1,Nt+1)~=0) );
                
                SnapClustW{2} = SnapClustW{1}(:,indep2 + (mod(indep2,Nt+1)~=0) );
                SnapClustW{1} = SnapClustW{1}(:,indep1 + (mod(indep1,Nt+1)~=0) );
                
                % non-linear term
                SnapClustNL{2} = SnapClustNL{1}(:,indep2 + (mod(indep2,Nt+1)~=0) );
                SnapClustNL{1} = SnapClustNL{1}(:,indep1 + (mod(indep1,Nt+1)~=0) );
                
                % nodes collection
                obj.tree{1} = {1,2};
                
                % binary tree branching
                for NC = (length(obj.tree)+2) : (obj.clusterNumber) 
                    
                    clear distvec
                    
                    % search maximum leaf error
                    for iNC = 1:NC-1
                        for iS = 1:size(SnapClust{iNC},2)
                            distvec{iNC}(iS) = pdist( [ SnapClust{iNC}(:,iS) , obj.centroids{iNC}*obj.centroids{iNC}'*SnapClust{iNC}(:,iS) ]' , 'squaredeuclidean'  );
                        end

                        [val(iNC),ind(iNC)] = max(distvec{iNC});
                    end

                    [~,indsplit] = max(val);

                    obj.centroids{NC} = SnapClust{indsplit}(:,ind(indsplit))./norm(SnapClust{indsplit}(:,ind(indsplit)));

                    indep1 = [];
                    indep2 = [];


                    for iS = 1:size(SnapClust{indsplit},2)
                        % distance wrt node indsplit
                        EP1 = pdist( [ SnapClust{indsplit}(:,iS) , obj.centroids{indsplit}*obj.centroids{indsplit}'*SnapClust{indsplit}(:,iS) ]' , 'squaredeuclidean'  );
                        % distance wrt node NC
                        EP2 = pdist( [ SnapClust{indsplit}(:,iS) , obj.centroids{NC}*obj.centroids{NC}'*SnapClust{indsplit}(:,iS) ]' , 'squaredeuclidean'  );

                        if EP1>=EP2
                            indep2 = [indep2, iS];
                        else
                            indep1 = [indep1, iS];
                        end
                    end

                    % snapshots branching                    
                    SnapClust{NC} = SnapClust{indsplit}(:,indep2  );
                    SnapClust{indsplit} = SnapClust{indsplit}(:,indep1 );

                    SnapClustW{NC} = SnapClustW{indsplit}(:,indep2  );
                    SnapClustW{indsplit} = SnapClustW{indsplit}(:,indep1 );
                    
                    SnapClustNL{NC} = SnapClustNL{indsplit}(:,indep2  );
                    SnapClustNL{indsplit} = SnapClustNL{indsplit}(:,indep1 );

                    % tree nodes update
                    obj.tree{NC-1} = {indsplit,NC};

                    
                end

                % loop over final clusters
                for iC = 1:obj.clusterNumber
                    
                    % compute POD
                    % solution
                    [obj.V{iC}, ~, sigma{iC}] = obj.POD(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    % recovery variable
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    % nonlinear terms for hyperreduction purposes
                    if nargin>3  
                        
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    
                end
                
                
            end
            
            
            
            if  strcmp(obj.clusterType,'Global')
                
                    % compute POD
                    % solution
                    [obj.V, ~, sigma] = obj.POD(uniquetol(S',5*1e-2,'Byrows',true)',PODtol);   
                    % recovery variable
                    [obj.Vw, ~, sigma] = obj.POD(uniquetol(Sw',5*1e-2,'Byrows',true)',PODtol);   
                    
                    % nonlinear terms for hyperreduction purposes
                    if nargin>3                    
                        
                        [obj.VD, ~ ,  ~] = obj.POD(uniquetol(SNL',5*1e-2,'Byrows',true)',DEIMtol);                    
                        [obj.iDEIM, obj.PMat, obj.VD] = DEIM(obj,obj.VD);
                    
                    end
                    
            end
            
        end
            
        function [U, V, sigma] = POD(obj,S,tolP)
            %POD proper orthogonal decomposition algorithm

            % svd method
            [U,Sigma,V] = svd(S,'econ');

            % singular values
            sigma = diag(Sigma);

            % truncation
            if tolP < 1
                sigmaC = cumsum(sigma.^2);
                [val,ind] = max( tolP^2>=(1-sigmaC./sigmaC(end)) ) ;
            else
                ind = min(tolP,size(V,2));
            end

            % final basis functions
            V = V(:,1:ind);
            U = U(:,1:ind);
                
            
        end


        function [selInd, PMat, PHI] = DEIM(obj,VD)
            % DEIM discrete emprical interpolation method
            
            % starting point
            [ ~ , selInd(1) ] = max( abs( VD(:,1) ) );
            
            PHI(:,1) = VD(:,1);
            PMat = zeros(size(VD));
            PMat(selInd(1),1) = 1;
            
            % greedy procedure
            for i = 2:size(VD,2)
            
                % interpolation
                coeffs = PHI(selInd,:)\VD(selInd,i);
                
                %residual
                residual = VD(:,i) - PHI*coeffs;
                
                % selection of new index
                [ ~ , selInd(i) ] = max( abs( residual ) );
            
                % update basis functions
                PHI(:,i) = VD(:,i);
                PMat(selInd(i),i) = 1;
            end          
            
        end    
        
    end
end

