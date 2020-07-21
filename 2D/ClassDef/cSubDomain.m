classdef cSubDomain
    properties
        Nx; Ny;
        dx; dy;
                
        NodesIdx;
        DisjointNodesIdx;
        FinalLayerNodesIdx;
        BCNodesIdx;
        GhostNodesIdx;
        ActiveNodesIdx;
        
        Node = cNode;
        
        LocalActiveBCNodesIdx;
        
        InterpolationMatrix;
        LocalOperatorMatrix;
        
        DisjointExtensionOperator;
        OverlappingExtensionOperator;
    end
    methods
        %% SETTING UP THE DOMAIN
        function obj = SETUP(obj, XGrid, YGrid, GlobalMask, LocalMask, NoOverlap, InterpolationOrder, c)
            obj = obj.BUILD_SUBDOMAIN(XGrid, YGrid, GlobalMask, LocalMask, NoOverlap, InterpolationOrder);
            obj = obj.CONSTRUCT_NODES(XGrid, YGrid);
            obj = obj.FIND_LOCAL_BC_NODES();
            obj = obj.CONSTRUCT_CPM_MATRICES(c);
            obj = obj.CONSTRUCT_DD_MATRICES(GlobalMask);
        end
        %% PARTITIONING
        function obj = BUILD_SUBDOMAIN(obj, XGrid, YGrid, GlobalMask, LocalMask, NoOverlap, InterpolationOrder)
            % USEFULL VARIABLES
            DisjointMask = LocalMask;
            obj.Nx = length(XGrid); obj.dx = XGrid(1, 2) - XGrid(1, 1);
            obj.Ny = length(YGrid); obj.dy = YGrid(1, 1) - YGrid(2, 1);
            InterpolationArray = (-InterpolationOrder/2:1:InterpolationOrder/2);
            InterpolationStencil = obj.Ny*ones(InterpolationOrder+1, 1)*InterpolationArray + InterpolationArray'*ones(1, InterpolationOrder+1);
            
            % OVERLAP
            for i = 1:NoOverlap
                LocalMask = (bwdist(LocalMask) <= sqrt(2)) .* (GlobalMask - LocalMask) + LocalMask;
            end
            % INTERFACES ALIGNMENT
            LocalMaskOld = zeros(size(LocalMask));
            InterfaceNodesIdx = find(LocalMask, 1);
            while norm(LocalMaskOld(InterfaceNodesIdx) - LocalMask(InterfaceNodesIdx), inf) ~= 0
                LocalMaskOld = LocalMask;
                NotLocalMask = GlobalMask - LocalMask;
                InterfaceNodes = ((bwdist(NotLocalMask) == 1) .* LocalMask) + ((bwdist(LocalMask) == 1) .* NotLocalMask);
                InterfaceNodesIdx = find(InterfaceNodes);
                for i = 1:length(InterfaceNodesIdx)
                    id = InterfaceNodesIdx(i);
                    [X, Y] = CLOSESTPOINT(XGrid(id), YGrid(id));
                    k = dsearchn([XGrid(:), YGrid(:)], [X, Y]);
                    LocalMask(id) = LocalMask(k);
                end
            end
            
            % COMPLETEING INTERPOLATION STENCIL
            NotLocalMask = GlobalMask - LocalMask;
            LocalFinalLayerNodesMask = ((bwdist(NotLocalMask) <= sqrt(2)) .* LocalMask);
            LocalFinalLayerNodesIdx = find(LocalFinalLayerNodesMask);            
            LocalMaskIdx = find(LocalMask);
            temp_X = XGrid(LocalMaskIdx); temp_Y = YGrid(LocalMaskIdx);
            for i = 1:length(LocalFinalLayerNodesIdx)
                id = LocalFinalLayerNodesIdx(i);
                [X, Y] = CLOSESTPOINT(XGrid(id), YGrid(id));
                k = dsearchn([temp_X, temp_Y], [X, Y]);
                GlobalIdx = LocalMaskIdx(k) + InterpolationStencil;
                LocalMask(GlobalIdx) = 0.5*(LocalMask(GlobalIdx) == 0) + LocalMask(GlobalIdx);
            end
            % ADDING GHOST NODES
            LocalMask = 0.5 * (bwdist(LocalMask == 0.5) == 1) .* (GlobalMask - LocalMask) + LocalMask;
            
            % PROPERTIES
            LocalGhostNodesMask = ((bwdist(~LocalMask) <= 1) .* LocalMask);
            LocalGhostNodesIdx = find(LocalGhostNodesMask);
            LocalActiveNodesMask = LocalMask - LocalGhostNodesMask;
            LocalActiveNodesIdx = find(LocalActiveNodesMask);
            LocalMaskIdx = [LocalActiveNodesIdx; LocalGhostNodesIdx];
            
            obj.NodesIdx = LocalMaskIdx;
            obj.DisjointNodesIdx = find(DisjointMask);
            obj.FinalLayerNodesIdx = LocalFinalLayerNodesIdx;
            obj.BCNodesIdx = find(LocalMask == 0.5);
            obj.GhostNodesIdx = LocalGhostNodesIdx;
            obj.ActiveNodesIdx = LocalActiveNodesIdx;
        end
        %% NODE CONSTRUCTIiON
        % TO-DO: Has to be compatible with different interpolation
        % orders.
        function obj = CONSTRUCT_NODES(obj, XGrid, YGrid)
            % USEFUL PARAMETERS
            Na = length(obj.ActiveNodesIdx);
            N = length(obj.NodesIdx);
            
            % NODES CONSTRUCTION
            obj.Node(N) = cNode;            
            % ACTIVE NODES
            for i = 1:Na
                id = obj.NodesIdx(i);
                obj.Node(i).ID = id;
                obj.Node(i).x = XGrid(id);
                obj.Node(i).y = YGrid(id);
                
                % FINDING NEIGHBOURS
                T = find(obj.NodesIdx == id - 1, 1);
                L = find(obj.NodesIdx == id - obj.Ny, 1);
                M = i;
                R = find(obj.NodesIdx == id + obj.Ny, 1);
                B = find(obj.NodesIdx == id + 1, 1);
                obj.Node(i).Neighbours = [T, L M R, B];
                obj.Node(i).Weights = [1/obj.dy^2, 1/obj.dx^2 -2/obj.dx^2-2/obj.dy^2  1/obj.dx^2, 1/obj.dy^2];
                obj.Node(i).Ghost = 0;
                
                % CLOSEST POINT REPRESENTATION
                [obj.Node(i).CPx, obj.Node(i).CPy] = CLOSESTPOINT(obj.Node(i).x, obj.Node(i).y);
            end
            % GHOST NODES
            for i = Na+1:N
                id = obj.NodesIdx(i);
                obj.Node(i).ID = id;
                obj.Node(i).x = XGrid(id);
                obj.Node(i).y = YGrid(id);
                
                % FINDING NEIGHBOURS
                obj.Node(i).Ghost = 1;
                
                % CLOSEST POINT REPRESENTATION
                [obj.Node(i).CPx, obj.Node(i).CPy] = CLOSESTPOINT(obj.Node(i).x, obj.Node(i).y);
            end
            
            % CLOSEST POINT REPRESENTATION MODIFICATION
            LocalFinalLayerNodesIdx = zeros(length(obj.FinalLayerNodesIdx), 1);
            LocalBCNodesIdx = zeros(length(obj.BCNodesIdx), 1);
            p = 1; q = 1;
            for i = 1:N
                id = obj.Node(i).ID;
                if ismember(id, obj.FinalLayerNodesIdx)
                    LocalFinalLayerNodesIdx(p) = i;
                    p = p + 1;
                end
                if ismember(id, obj.BCNodesIdx)
                    LocalBCNodesIdx(q) = i;
                    q = q + 1;
                end
            end
            XVals = [obj.Node(LocalFinalLayerNodesIdx).CPx];
            YVals = [obj.Node(LocalFinalLayerNodesIdx).CPy];
            for i = 1:length(LocalBCNodesIdx)
                id = LocalBCNodesIdx(i);
                k = dsearchn([XVals', YVals'], [obj.Node(id).x, obj.Node(id).y]);
                obj.Node(id).CPx = obj.Node(LocalFinalLayerNodesIdx(k)).CPx;
                obj.Node(id).CPy = obj.Node(LocalFinalLayerNodesIdx(k)).CPy;
            end
            
            % INTERPOLATION STENCIL
            XVals = [obj.Node.x]; YVals = [obj.Node.y];
            for i = 1:N
                k = dsearchn([XVals', YVals'], [obj.Node(i).CPx, obj.Node(i).CPy]);
                L = obj.Node(k).Neighbours(2); LL = obj.Node(L).Neighbours(2);
                M = k;
                R = obj.Node(k).Neighbours(4); RR = obj.Node(R).Neighbours(4);
                temp = ones(5, 1)*[LL L M R RR] + (-2:1:2)'*ones(1, 5);
                temp = temp(:);
                obj.Node(i).CPNeighbours = temp';
                X = obj.Node(i).CPx - obj.Node(k).x; Y = obj.Node(i).CPy - obj.Node(k).y;
                temp = (2:-1:-2)'*obj.dy;
                P = [temp.^4 temp.^3 temp.^2 temp ones(5, 1)];
                temp = (-2:1:2)*obj.dx;
                Q = [temp.^4; temp.^3; temp.^2; temp; ones(1, 5)];
                A1 = [Y^4 Y^3 Y^2 Y 1]/P; A2 = Q\[X^4; X^3; X^2; X; 1];
                temp = A1'*A2';
                temp = temp(:);
                obj.Node(i).CPWeights = temp';
            end
        end
        %% FINDING LOCAL BC NODES
        function obj = FIND_LOCAL_BC_NODES(obj)
            N = length(obj.NodesIdx);
            ActiveBCNodesIdx = setdiff([obj.BCNodesIdx; obj.FinalLayerNodesIdx], obj.GhostNodesIdx);
            obj.LocalActiveBCNodesIdx = zeros(length(ActiveBCNodesIdx), 1);   p = 1;
            for i = 1:N
                id = obj.Node(i).ID;
                if ismember(id, ActiveBCNodesIdx)
                    obj.LocalActiveBCNodesIdx(p) = i;
                    p = p + 1;
                end
            end
        end
        %% MATRICES CONSTRUCTION
        function obj = CONSTRUCT_CPM_MATRICES(obj, c)
            % USEFUL VARIABLES
            N = length(obj.NodesIdx);
            Na = length(obj.ActiveNodesIdx);
            
            % LOCAL OPERATOR MATRIX CONSTRUCTION
            E = zeros(N, Na);
            for i = 1:N
                E(i, obj.Node(i).CPNeighbours) = obj.Node(i).CPWeights;
            end
            obj.InterpolationMatrix = E;
            Delta_h = zeros(Na, N);
            for i = 1:Na
                Delta_h(i, obj.Node(i).Neighbours) = obj.Node(i).Weights;
            end
            temp = diag(diag(Delta_h));
            obj.LocalOperatorMatrix = (temp + (Delta_h - [temp zeros(Na, N - Na)])*E);
            obj.LocalOperatorMatrix = c*eye(Na) - obj.LocalOperatorMatrix;
            % TRANSMISSION CONDITION
            obj.LocalOperatorMatrix(obj.LocalActiveBCNodesIdx, :) = 0;
            obj.LocalOperatorMatrix(obj.LocalActiveBCNodesIdx, obj.LocalActiveBCNodesIdx) = eye(length(obj.LocalActiveBCNodesIdx));
        end
        %% DD MATRICES
        function obj = CONSTRUCT_DD_MATRICES(obj, GlobalMask)
            % USEFUL VARIABLES
            GhostNodesMask = ((bwdist(~GlobalMask) <= 1) .* GlobalMask);
            ActiveNodesMask = GlobalMask - GhostNodesMask;
            ActiveGlobalNodesIdx = find(ActiveNodesMask);
            NN = length(ActiveGlobalNodesIdx);
            
            % OVERLAPPING EXTENSION OPERATOR
            N = length(obj.ActiveNodesIdx);
            LocalActiveOverlappingNodesIdx = zeros(N, 1); p = 1;
            for i = 1:NN
                if ismember(ActiveGlobalNodesIdx(i), obj.ActiveNodesIdx)
                    LocalActiveOverlappingNodesIdx(p) = i;
                    p = p + 1;
                end
            end
            temp = eye(NN);
            obj.OverlappingExtensionOperator = temp(LocalActiveOverlappingNodesIdx, :);
            
            % DISJOINT EXTENSION OPERATOR
            ActiveDisjointNodesIdx = setdiff(obj.DisjointNodesIdx, obj.GhostNodesIdx);
            N = length(ActiveDisjointNodesIdx);
            LocalActiveDisjointNodesIdx = zeros(N, 1); p = 1;
            for i = 1:NN
                if ismember(ActiveGlobalNodesIdx(i), ActiveDisjointNodesIdx)
                    LocalActiveDisjointNodesIdx(p) = i;
                    p = p + 1;
                end
            end
            temp = eye(NN);
            temp1 = setdiff(LocalActiveOverlappingNodesIdx, LocalActiveDisjointNodesIdx);
            temp(temp1, :) = 0;
            obj.DisjointExtensionOperator = temp(LocalActiveOverlappingNodesIdx, :);
        end
    end
end