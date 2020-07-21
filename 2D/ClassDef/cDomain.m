classdef cDomain
    properties
        Nx; Ny;
        dx; dy;
                
        NodesIdx;
        GhostNodesIdx;
        ActiveNodesIdx;
        
        Node = cNode;
        
        InterpolationMatrix;
        OperatorMatrix;
        RHSMatrix;
    end
    methods
        %% SETTING UP THE DOMAIN
        function obj = SETUP(obj, XGrid, YGrid, GlobalMask, c, f)
            obj = obj.BUILD_DOMAIN(XGrid, YGrid, GlobalMask);
            obj = obj.CONSTRUCT_NODES(XGrid, YGrid);
            obj = obj.CONSTRUCT_CPM_MATRICES(c, f);
        end
        %% PARTITIONING
        function obj = BUILD_DOMAIN(obj, XGrid, YGrid, GlobalMask)
            obj.Nx = length(XGrid); obj.dx = XGrid(1, 2) - XGrid(1, 1);
            obj.Ny = length(YGrid); obj.dy = YGrid(1, 1) - YGrid(2, 1);
            LocalGhostNodesMask = ((bwdist(~GlobalMask) <= 1) .* GlobalMask);
            LocalGhostNodesIdx = find(LocalGhostNodesMask);
            LocalActiveNodesMask = GlobalMask - LocalGhostNodesMask;
            LocalActiveNodesIdx = find(LocalActiveNodesMask);
            LocalMaskIdx = [LocalActiveNodesIdx; LocalGhostNodesIdx];
            
            obj.NodesIdx = LocalMaskIdx;
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
        %% MATRICES CONSTRUCTION
        function obj = CONSTRUCT_CPM_MATRICES(obj, c, f)
            % USEFUL VARIABLES
            N = length(obj.NodesIdx);
            Na = length(obj.ActiveNodesIdx);
            
            % OPERATOR MATRIX CONSTRUCTION
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
            obj.OperatorMatrix = (temp + (Delta_h - [temp zeros(Na, N - Na)])*E);
            obj.OperatorMatrix = c*eye(Na) - obj.OperatorMatrix;
            
            % RHS-MATRIX CONSTRUCTION
            [theta, ~] = cart2pol([obj.Node(1:Na).x]', [obj.Node(1:Na).y]');
            obj.RHSMatrix = f(theta);
        end
    end
end