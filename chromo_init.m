function chromosome = chromo_init(N,M,V)
%% function chromosome = chromo_init(N,M,V,min_range,max_range)
% This function initializes the chromosomes. Each chromosome has the
% following at this stage
%       * set of decision variables
%       * objective function values
%
% where,
% N - Population size
% M - Number of objective functions
% V - Number of decision variables, in here V equals to number of nodes
% min_range - A vector of decimal values which indicate the minimum value
% for each decision variable.
% max_range - Vector of maximum possible values for decision variables.


%% Chromosome Format
% [[path],[length,corr_counts]]: where
% path: 1-by-613 array, which starts from 1, ends with 613, rest of path
% fulfilled with 0. eg: [1,34,526,83,92,613,0,...,0]
% length: scaler, indicates length of path, which is sum of each edge's
% distance.
% corr_counts: scaler, equal to number of inner nodes, eg: above path's
% corrcounts is 4.
chromosome=zeros(N,V+M);
% Nh=round(N/2);% number of hor nodes
% Nv=N-Nh;% number of ver nodes


air_info=[1,0,0,0,0];% currentNodeNo;err_hor;err_ver;long;node_counts;
global gf nn fpath;

vPath=struct('validpath',[]);
Id=1;
Index=1;
Treef=struct('Id',Id,'Parent',-1,'Child',[],'Dead',false,'Info',air_info);% Rootf
% ppath=1;
stopf=false;
N_tmp=N;

% DFS search tree
while ~stopf
    %     Ncs=Treef(end).Id;% ender  ID at bottom height
    %     childNumber=0;
    if Index<1
        stopf=true;
    end
    % current node is ender
    if Treef(Index).Info(1)==nn
        N_tmp=N_tmp-1;% dec number of path need to find
        %         stopf=true;% find one path will stop
        
        % find current path
        ipath=Index;
        ppath=Treef(Index).Info(1);% end node of current path
        while ~isequal(ppath(1),1)
            ipath=[Treef(ipath(1)).Parent,ipath];
            ppath=[Treef(ipath(1)).Info(1),ppath];
        end
        vPath=[vPath,struct('validpath',ppath)];
        
        % let last node of current path dead
        Index=Treef(Index).Parent;
        Treef(Index).Dead=true;
        Index=ipath(randi([2,length(ipath)-1]));
        
        % draw finded path
        set(fpath,'XData',gf(ppath,2),'YData',gf(ppath,3),'ZData',gf(ppath,4));
        drawnow;
    else % current node not ender
        % current node is dead
        if Treef(Index).Dead
            Index=Treef(Index).Parent;
        else % current node not dead
            % current node no child
            if isempty(Treef(Index).Child)
                %                 find current path
                ipath=Index;
                ppath=Treef(ipath).Info(1);% end node of current path
                while ~isequal(ppath(1),1)
                    ipath=Treef(ipath).Parent;
                    ppath=[Treef(ipath).Info(1),ppath];
                end
                
                % draw finded path
                set(fpath,'XData',gf(ppath,2),'YData',gf(ppath,3),'ZData',gf(ppath,4));
                drawnow;
                
                childLeafs=se_vnode(Treef(Index).Info,ppath);
                [childNumber,~]=size(childLeafs);
                
                % update tree
                for loop=1:childNumber
                    Id=Id+1;
                    Leaff=struct('Id',Id,'Parent',Index,'Child',[],'Dead',false,'Info',childLeafs(loop,:));
                    Treef(Id)=Leaff;
                    Treef(Index).Child=union(Treef(Index).Child,Id);
                end
            else
                
            end
            
            ChildNode=Treef(Index).Child;% ChildNode
            vChildNode=ChildNode(arrayfun(@(x)x.Dead==false,Treef(ChildNode)));% alive ChildNode
            
            if isempty(vChildNode)% all child dead
                Treef(Index).Dead=true;
                Index=Treef(Index).Parent;
            else % not all child dead
                vEndNodeNo=vChildNode(arrayfun(@(x)x.Info(1)==nn,Treef(vChildNode)));
                if ~isempty(vEndNodeNo)% find ender in range
                    Index=vEndNodeNo;
                else % ender is not in range
                    % sort vChildNode by distance (current to ender)
                    %                     vCNvec=zeros(length(vChildNode),3);% valid ChildNode vec
                    vCNc=zeros(length(vChildNode),1);% valid ChildNode distance
                    for loopi=1:length(vChildNode)
                        vCNc(loopi)=norm(gf(nn,2:4)-gf(Treef(vChildNode(loopi)).Info(1),2:4),2);
                        %                         vCNvec(loopi,:)=(gf(nn,2:4)-gf(Treef(vChildNode(loopi)).Info(1),2:4))/norm(gf(nn,2:4)-gf(Treef(vChildNode(loopi)).Info(1),2:4),2);
                        %                         vCNc(loopi)=vCNvec(loopi,:)*transpose(vec_ss);
                    end
                    [~,vCNi]=sort(vCNc,'ascend');
                    %                     [vCNsorted,vCNi]=sort(vCNc,'descend');
                    % choose distance least as next
                    %                     Index=vChildNode(vCNi(1));
                    % randomly choose one in half distance least as next
                    Index=vChildNode(vCNi(randi([1,round(length(vCNi)/2)])));
                    %                     Index=vChildNode(vCNi(randi([1,length(vCNi)])));
                end
                
            end
            
        end
        
    end
    
    if isequal(N_tmp,0)
        stopf=true;
    end
    
end


for loop=1:N
    chromosome(loop,1:length(vPath(loop+1).validpath))=vPath(loop+1).validpath;
    pathdata=result_out(vPath(loop+1).validpath);
    chromosome(loop,V+1:V+M)=pathdata(end,4:5);
end

end

