classdef localizedReduction
    %LOCALIZEDREDUCTION This class implements the localized reduction
    %strategies presented in the paper "Numerical approximation of
    %parametrized problems in cardiac electrophysiology by a local reduced basis method"
    
    
    properties
        V % struct containing the basis function related to the selected clusters
        centroids % struct containing the centroids of the clusters   
        clusterType % flag identifying the possible cluster technique
        clusterNumber % integer identifying the number of clusters 
    end
    
    methods
        function obj = localizedReduction(inputType,inputNclust)
            %LOCALIZEDREDUCTION Construct an instance of this class
            %   Detailed explanation goes here
            obj.clusterType   = inputType;
            obj.clusterNumber = inputNclust;            
        end
        
        function obj = offline(obj,nTrain)
            %OFFLINE offline phase for the construction of the cluster and
            %the related basis functions
           
            % model inputs
            param(1) = 1; % domain lenght
            param(2) = 0.015;  % conducibility
            param(3) = 0.5;   % 
            param(4) = 2;
            
            Nt = 399;

            FNS = FNSolver(param, 1024, 0, 2, Nt );
            
            
            %Snapshots matrix
            S = [];
            
            
            for iT = 1 : nTrain
                SParam(iT) = 0.003 + rand(1,1)*(0.05-0.003);
                [u,w] = FNS.solveFOM( SParam(iT) );
                S = [S , u];
                %S
            end
            
            
            if  strcmp(obj.clusterType,'kmeansState')
                [IDX, C] = kmeans(S', obj.clusterNumber)

                
                for iC = 1:obj.clusterNumber
                    obj.centroids{iC} = C(iC,:)';
                    
                    SnapClust{iC} = S(: , find(IDX==iC));
                end
                
            end
            

            
            if  strcmp(obj.clusterType,'kmeansParam')
                [IDX, C] = kmeans(SParam', obj.clusterNumber)
                
                for iC = 1:obj.clusterNumber
                    obj.centroids{iC} = C(iC,:)';
                    
                    SnapClust{iC} = S(: , find(IDX==iC));
                end
                
            end
            
            if  strcmp(obj.clusterType,'Time')
                
                
                for iC = 1:obj.clusterNumber
                    obj.centroids{iC} = [1+(iC-1)*400/obj.clusterNumber,(iC)*400/obj.clusterNumber]
                    
                    %%
                    SnapClust{iC} = [];
                    
                    for iT = 1 : nTrain
                        [(iT-1)*(Nt+1)+obj.centroids{iC}(1) ,  (iT-1)*(Nt+1)+obj.centroids{iC}(2) ]
                        SnapClust{iC} =  [ SnapClust{iC}  , S(:, (iT-1)*(Nt+1)+obj.centroids{iC}(1): (iT-1)*(Nt+1)+obj.centroids{iC}(2) ) ];
                    end
                end
                                   
            end
            
            if  strcmp(obj.clusterType,'Global')
                    [U, V, sigma] = obj.POD(S,1e-1);
            end
                
        end
            
        function [U, V, sigma] = POD(obj,S,tolP)
            %POD proper orthogonal decomposition algorithm
            
            % correlation matrix method
%                 [ns, Nh] = size(S);  % we assume that the snapshots are organized in column
%                 if ns < Nh
%                     % create the correlation matrix
%                     C = S'*S;
%                     % solve eigenvalue problem
%                     [V, sigma, ~]     = svd(C);
%                 else
%                     % create the correlation matrix
%                     C = S*S';
%                     % solve eigenvalue problem
%                     [V, sigma, ~]     = svd(C);
%                 end

            % svd method
                [U,Sigma,V] = svd(S,'econ');
                
                sigma = diag(Sigma);
                
                if tolP < 1
                    sigmaC = cumsum(sigma.^2);
                    [val,ind] = max( tolP^2>=(1-sigmaC./sigmaC(end)) ) ;
                    keyboard
                else
                    ind = min(tolP,size(V,2));
                end
                
                V = V(:,1:ind);
                U = U(:,1:ind);
                
            
        end

        
    end
end

