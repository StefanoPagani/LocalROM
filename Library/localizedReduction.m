classdef localizedReduction
    %LOCALIZEDREDUCTION This class implements the localized reduction
    %strategies presented in the paper "Numerical approximation of
    %parametrized problems in cardiac electrophysiology by a local reduced basis method"
    
    
    properties
        V % struct containing the basis function related to the selected clusters
        centroids % struct containing the centroids of the clusters   
        clusterType % flag identifying the possible cluster technique
        clusterNumber % integer identifying the number of clusters
        tree % PEBL structure
        VD % struct containing non-linear term basis functions
        iDEIM % struct containing DEIM indices
        PMat % matrix
        Vw % struct containing the basis function related to the selected clusters for the recovery variable
    end
    
    methods
        function obj = localizedReduction(inputType,inputNclust)
            %LOCALIZEDREDUCTION Construct an instance of this class
            %   Detailed explanation goes here
            obj.clusterType   = inputType;
            obj.clusterNumber = inputNclust; 
        end
        
        function obj = offline(obj,nTrain,PODtol,DEIMtol)
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
            Sw = [];
            
            rng('default')
            
            for iT = 1 : nTrain
                
                SParam(iT) = 0.003 + (iT-1)/(nTrain-1) *(0.05-0.003);
                
                %SParam(iT) = 0.003 + rand(1,1)*(0.05-0.003);
                %SParam(iT) = 0.02 + rand(1,1)*(0.05-0.02);
                [u,w] = FNS.solveFOM( SParam(iT) );
                S = [S , u];
                Sw = [Sw , w];
                %S
            end
            
            % non linear term snapshots
            f = @(v)  v.*(v-0.1).*(v-1); 
            SNL = f(S);
            
%             plot(SParam,ones(size(SParam)),'*')
%             keyboard
            
            
            if  strcmp(obj.clusterType,'kmeansState')
                [IDX, C] = kmeans(S', obj.clusterNumber, 'MaxIter', 1000 );

                
%                 keyboard
                
                for iC = 1:obj.clusterNumber
                    obj.centroids{iC} = C(iC,:)';
                    
                    % identify the snapshots related to the current cluster
                    selSnap = find(IDX==iC) ; 
                    SnapClust{iC} = S(: , selSnap + (mod(selSnap,Nt+1)~=0)  );
                    
                    SnapClustW{iC} = Sw(: , selSnap + (mod(selSnap,Nt+1)~=0)  );
                    
                    SnapClustNL{iC} = SNL(: , selSnap + (mod(selSnap,Nt+1)~=0)  );
                    
%                     keyboard
                    % compute POD
                    [obj.V{iC}, ~ , sigma{iC}] = obj.POD(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    if nargin>3  
                    
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);
                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    
                end
                
                
                
            end
            

            
            if  strcmp(obj.clusterType,'kmeansParam')
                [IDX, C] = kmeans((SParam'), obj.clusterNumber);
                
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
                        
                         
                        %[X,Y] = meshgrid( linspace(0,FNS.L, FNS.Nh+1), linspace(FNS.t0,FNS.tF, FNS.Nt+1)  );

%                         surf( X, Y ,  S(:, 1 + (iT-1)*(Nt+1) : (iT)*(Nt+1) )' )
%                         shading interp
%                         iC
%                         pause(1)
                    end
                    
%                     keyboard
    
                    %[size(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',2), size(SnapClust{iC},2) ]

                    % compute POD
                    [obj.V{iC}, ~ , sigma{iC}] = obj.POD( uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    if nargin>3  
                    
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);
                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
%                     keyboard
                end
                
            end
            
            if  strcmp(obj.clusterType,'Time')
                
                
                for iC = 1:obj.clusterNumber
                    obj.centroids{iC} = [1+(iC-1)*400/obj.clusterNumber,(iC)*400/obj.clusterNumber];
                    
                    %%
                    SnapClust{iC} = [];
                    SnapClustW{iC} = [];
                    SnapClustNL{iC} = [];
                    
                    for iT = 1 : nTrain
%                         keyboard
                        selSnap = [(iT-1)*(Nt+1)+obj.centroids{iC}(1) :  (iT-1)*(Nt+1)+obj.centroids{iC}(2) ];
                        SnapClust{iC} =  [ SnapClust{iC}  , S(:, selSnap  ) ];
                    
                        SnapClustW{iC} =  [ SnapClustW{iC}  , Sw(:, selSnap  ) ];
                        
                        SnapClustNL{iC} =  [ SnapClustNL{iC}  , SNL(:, selSnap  ) ];
                        
                    end
                    
%                     keyboard
                    
                    
                    [obj.V{iC}, ~, sigma{iC}] = obj.POD(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    if nargin>3  
                    
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);
                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    %[obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                end
                                   
            end
            
            if  strcmp(obj.clusterType,'PEBL')
                
                for iS = 1:size(S,2)
                    normvec(iS) = norm( S(:,iS) );
                end
                [val,ind] = max(normvec);
                
                obj.centroids{1} = S(:,ind)./norm(S(:,ind));
                
                SnapClust{1} = S;
                
                SnapClustW{1} = Sw;
                
                SnapClustNL{1} = SNL;
                
                
                for iS = 1:size(SnapClust{1},2)
                    distvec(iS) = pdist( [ SnapClust{1}(:,iS) , obj.centroids{1}*obj.centroids{1}'*SnapClust{1}(:,iS)  ]' , 'squaredeuclidean'  );
                end
                
                [val,ind] = max(distvec);
                
                obj.centroids{2} = SnapClust{1}(:,ind)./norm(SnapClust{1}(:,ind));
                
                indep1 = [];
                indep2 = [];
                
                
                for iS = 1:size(SnapClust{1},2)
                    EP1 = pdist( [ SnapClust{1}(:,iS) , obj.centroids{1}*obj.centroids{1}'*SnapClust{1}(:,iS) ]' , 'squaredeuclidean'  );
                    EP2 = pdist( [ SnapClust{1}(:,iS) , obj.centroids{2}*obj.centroids{2}'*SnapClust{1}(:,iS) ]', 'squaredeuclidean'  );
                
                    if EP1>=EP2
                        indep2 = [indep2, iS];
                    else
                        indep1 = [indep1, iS];
                    end
                end
                
                SnapClust{2} = SnapClust{1}(:,indep2 + (mod(indep2,Nt+1)~=0) );
                SnapClust{1} = SnapClust{1}(:,indep1 + (mod(indep1,Nt+1)~=0) );
                
                SnapClustW{2} = SnapClustW{1}(:,indep2 + (mod(indep2,Nt+1)~=0) );
                SnapClustW{1} = SnapClustW{1}(:,indep1 + (mod(indep1,Nt+1)~=0) );
                
                % non-linear term
                SnapClustNL{2} = SnapClustNL{1}(:,indep2 + (mod(indep2,Nt+1)~=0) );
                SnapClustNL{1} = SnapClustNL{1}(:,indep1 + (mod(indep1,Nt+1)~=0) );
                
                
                obj.tree{1} = {1,2};
                
                
%                 keyboard
                for NC = (length(obj.tree)+2) : (obj.clusterNumber) 
                    
                    clear distvec
                    
                    % search maximum leaf error
                    for iNC = 1:NC-1
%                         keyboard
                        for iS = 1:size(SnapClust{iNC},2)
                            distvec{iNC}(iS) = pdist( [ SnapClust{iNC}(:,iS) , obj.centroids{iNC}*obj.centroids{iNC}'*SnapClust{iNC}(:,iS) ]' , 'squaredeuclidean'  );
                        end

                        [val(iNC),ind(iNC)] = max(distvec{iNC});
                    end

                    [~,indsplit] = max(val);
%                     keyboard
                    obj.centroids{NC} = SnapClust{indsplit}(:,ind(indsplit))./norm(SnapClust{indsplit}(:,ind(indsplit)));

                    indep1 = [];
                    indep2 = [];


                    for iS = 1:size(SnapClust{indsplit},2)
                        EP1 = pdist( [ SnapClust{indsplit}(:,iS) , obj.centroids{indsplit}*obj.centroids{indsplit}'*SnapClust{indsplit}(:,iS) ]' , 'squaredeuclidean'  );
                        EP2 = pdist( [ SnapClust{indsplit}(:,iS) , obj.centroids{NC}*obj.centroids{NC}'*SnapClust{indsplit}(:,iS) ]' , 'squaredeuclidean'  );

                        if EP1>=EP2
                            indep2 = [indep2, iS];
                        else
                            indep1 = [indep1, iS];
                        end
                    end

%                     keyboard
                    
                    SnapClust{NC} = SnapClust{indsplit}(:,indep2  );
                    SnapClust{indsplit} = SnapClust{indsplit}(:,indep1 );

                    SnapClustW{NC} = SnapClustW{indsplit}(:,indep2  );
                    SnapClustW{indsplit} = SnapClustW{indsplit}(:,indep1 );
                    
                    SnapClustNL{NC} = SnapClustNL{indsplit}(:,indep2  );
                    SnapClustNL{indsplit} = SnapClustNL{indsplit}(:,indep1 );

                    
                    obj.tree{NC-1} = {indsplit,NC};

                    
                end
                
%                 keyboard
                
                for iC = 1:obj.clusterNumber
                    
                    [obj.V{iC}, ~, sigma{iC}] = obj.POD(uniquetol(SnapClust{iC}',5*1e-2,'Byrows',true)',PODtol);
                    
                    [obj.Vw{iC}, ~, ~] = obj.POD(uniquetol(SnapClustW{iC}',5*1e-2,'Byrows',true)',PODtol);   
                    
                    if nargin>3  
                        
                        [obj.VD{iC}, ~ ,  ~] = obj.POD(uniquetol(SnapClustNL{iC}',5*1e-2,'Byrows',true)',DEIMtol);
                    
                        [obj.iDEIM{iC}, obj.PMat{iC}, obj.VD{iC}] = DEIM(obj,obj.VD{iC});
                    
                    end
                    
                end
                
                
            end
            
            
            
            if  strcmp(obj.clusterType,'Global')
                
                
                    [obj.V, ~, sigma] = obj.POD(uniquetol(S',5*1e-2,'Byrows',true)',PODtol);   
                    
                    [obj.Vw, ~, sigma] = obj.POD(uniquetol(Sw',5*1e-2,'Byrows',true)',PODtol);   
                    
                    if nargin>3                    
                        
                        [obj.VD, ~ ,  ~] = obj.POD(uniquetol(SNL',5*1e-2,'Byrows',true)',DEIMtol);
                    
                        [obj.iDEIM, obj.PMat, obj.VD] = DEIM(obj,obj.VD);
                    
                    end
%                     keyboard
                    
            end
                
            %obj.V
            
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
                else
                    ind = min(tolP,size(V,2));
                end
                
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

