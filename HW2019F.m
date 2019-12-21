%% Setup the Import Options

clear
choice=input('Choose 1 or 2 problem to be solved:');
% opts = spreadsheetImportOptions("NumVariables", 6);
%
% % Specify sheet and range
% opts.Sheet = "data1";
% opts.DataRange = "A3:F615";
%
% % Specify column names and types
% opts.VariableNames = ["NodeNum", "Xm", "Ym", "Zm", "CorrType", "Mark"];
% opts.SelectedVariableNames = ["NodeNum", "Xm", "Ym", "Zm", "CorrType", "Mark"];
% opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
%
% % Import the data
% tab_f1 = readtable("./data/f1.xlsx", opts, "UseExcel", false);
%% Clear temporary variables

% clear opts
%% Read Raw Data

% save('./data/f2.txt','f2','-ascii');
load('./data/f1.txt','f1','-ascii');
load('./data/f2.txt','f2','-ascii');

% gf: global f array
% nn: number of nodes
global gf nn; 
if choice==1
    gf=f1;
    nn=length(gf(:,1));
elseif choice==2
    gf=f2;
    nn=length(gf(:,1));
else
    fprintf('Wrong input!\n');
    return
end
%% Parameters

global parmbnd theta delta;
if choice==1
    parmbnd=[25;15;20;25;];% alpha1=25;alpha2=15;beta1=20;beta2=25;
    theta=30;
    delta=0.001;
elseif choice==2
    parmbnd=[20;10;15;20];% alpha1=20;alpha2=10;beta1=15;beta2=20;
    theta=20;
    delta=0.001;
end
%% Draw Figure

if choice==1
    figure("Name","F1 3D Graph");
elseif choice==2
    figure("Name","F2 3D Graph");
end
[row_nan,col_nan]=find(isnan(gf)==1);
scatter3(gf(row_nan,2),gf(row_nan,3),gf(row_nan,4),'*k');
hold on
[row_hor,col_hor]=find(gf(:,5)==0);
scatter3(gf(row_hor,2),gf(row_hor,3),gf(row_hor,4),'.r');
[row_ver,col_ver]=find(gf(:,5)==1);
scatter3(gf(row_ver,2),gf(row_ver,3),gf(row_ver,4),'.b');
% hold off
%% Search Tree

air_info=[1,0,0,0,0];% currentNodeNo;err_hor;err_ver;long;node_counts;

% d1=zeros(length(gf));% distance bwn 2 nodes
% for lp=1:length(d1)
%     for lpi=(lp+1):length(d1)
%         d1(lp,lpi)=norm(gf(lp,2:4)-gf(lpi,2:4),2);
%     end
% end
% d1=d1+transpose(d1);

vPath=[];% vPath: valid Path set
N=10;% number of path need to find
Id=1;
Index=1;
hl=[1]; % num of leafs at height
vec_ss=(gf(nn,2:4)-gf(1,2:4))/norm(gf(nn,2:4)-gf(1,2:4),2);% vec of starter to end
Rootf=struct('Id',Id,'Parent',-1,'Child',[],'Dead',false,'Info',air_info);
Treef=[Rootf];
ppath=[1];
gpath=plot3(gf(ppath,2),gf(ppath,3),gf(ppath,4),'k');
stopf=false;

fprintf("parameters loaded, start searching\n");
% DFS search tree
while ~stopf
    %     Ncs=Treef(end).Id;% ender  ID at bottom height
    %     childNumber=0;
    
    % current node is ender
    if Treef(Index).Info(1)==nn
        N=N-1;% dec number of path need to find
        %         stopf=true;% find one path will stop
        
        % find current path
        ipath=Index;
        ppath=Treef(Index).Info(1);% end node of current path
        while ~isequal(ppath(1),1)
            ipath=Treef(ipath).Parent;
            ppath=[Treef(ipath).Info(1),ppath];
        end
        vPath(end+1).validpath=ppath;
        
        % let last node of current path dead
        Index=Treef(Index).Parent;
        Treef(Index).Dead=true;
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
                % update path graph now
                set(gpath,'XData',gf(ppath,2),'YData',gf(ppath,3),'ZData',gf(ppath,4));
                drawnow;
                
                childLeafs=se_vnode(Treef(Index).Info,ppath);
                [childNumber,tmp]=size(childLeafs);
                
                % update tree
                for loop=1:childNumber
                    Id=Id+1;
                    Leaff=struct('Id',Id,'Parent',Index,'Child',[],'Dead',false,'Info',childLeafs(loop,:));
                    Treef(Id)=Leaff;
                    Treef(Index).Child=[Treef(Index).Child,Id];
                end
            else
                
            end
            
            ChildNode=Treef(Index).Child;% ChildNode
            vChildNode=ChildNode(arrayfun(@(x)x.Dead==false,Treef(ChildNode)));% valid ChildNode
            
            if isempty(vChildNode)% all child dead
                Treef(Index).Dead=true;
                %                 Index=Treef(Index).Parent;
            else % not all child dead
                vEndNodeNo=vChildNode(arrayfun(@(x)x.Info(1)==nn,Treef(vChildNode)));
                if ~isempty(vEndNodeNo)% find ender in range
                    Index=vEndNodeNo;
                else % ender is not in range
                    % sort vChildNode by distance (current to ender)
%                     vCNvec=zeros(length(vChildNode),3);% valid ChildNode vec
                    vCNc=zeros(length(vChildNode),1);% valid ChildNode cos(angle)
                    for loopi=1:length(vChildNode)
                        vCNc(loopi)=norm(gf(nn,2:4)-gf(Treef(vChildNode(loopi)).Info(1),2:4),2);
%                         vCNvec(loopi,:)=(gf(nn,2:4)-gf(Treef(vChildNode(loopi)).Info(1),2:4))/norm(gf(nn,2:4)-gf(Treef(vChildNode(loopi)).Info(1),2:4),2);
%                         vCNc(loopi)=vCNvec(loopi,:)*transpose(vec_ss);
                    end
                    [vCNsorted,vCNi]=sort(vCNc,'ascend');
%                     [vCNsorted,vCNi]=sort(vCNc,'descend');
                    Index=vChildNode(vCNi(1));
                    %                     nextn=round(rand()*(length(vChildNode)-1))+1;
                    %                     Index=vChildNode(nextn);
                end
                
            end
            
        end
        
    end
    
    if isequal(N,0)
        stopf=true;
    end
    
end

% BFS search tree
% while ~stopf
%     Ncs=Treef(end).Id; % ender ID at bottom height
%     childLeaf=0;
%     for loop=hl(end)-1:-1:0 % leafs at bottom height
%         % Treef(Ncs-loop)  current leaf, start from left
%
%         % path of current leaf
%         pindex=Ncs-loop;% id of path node
%         ppath=[Treef(pindex).Info(1)];% current path
%         while ~isequal(ppath(1),1)
%             pindex=Treef(pindex).Parent;
%             ppath=[Treef(pindex).Info(1);ppath];
%         end
%
%         if ppath(end)==length(gf)
%             vPath=[vPath;struct('validpath',ppath)];
%             fprintf(vPath);
%         end
%
%         % search valid node of current leaf
%         vns=se_vnode(Treef(Ncs-loop).Info,ppath);
%         [r_vns,c_vns]=size(vns);
%         for vloop=1:c_vns
%             Id=Id+1;
%             Leaff=struct('Id',Id,'Parent',Ncs-loop,'Child',[],'Info',vns(2:end,vloop));
%             Treef=[Treef;Leaff];% add child of current leaf to Treef table
%             Treef(Ncs-loop).Child=[Treef(Ncs-loop).Child;Id];
%         end
%         childLeaf=childLeaf+c_vns;
%
%     end
%     hl=[hl;childLeaf];
%     fprintf(hl);
%
%     if childLeaf==0
%         stop=true;
%     end
% end
%% Finish

fprintf("Search Tree End!\n");

% save('vPathResult.mat','vPath');
%% Draw path

for loop=1:length(vPath)
%     subplot(1,length(vPath)-1,loop-1);
%     scatter3(gf(row_nan,2),gf(row_nan,3),gf(row_nan,4),'*k');
%     hold on
%     scatter3(gf(row_hor,2),gf(row_hor,3),gf(row_hor,4),'.r');
%     scatter3(gf(row_ver,2),gf(row_ver,3),gf(row_ver,4),'.b');
% 
%     plot3(gf(vPath(loop).validpath,2),gf(vPath(loop).validpath,3),gf(vPath(loop).validpath,4),'k');
    
    fprintf("#%d PATH:\n",loop);
    result_out(vPath(loop).validpath);
end

path=vPath(1).validpath;
set(gpath,'XData',gf(path,2),'YData',gf(path,3),'ZData',gf(path,4));
drawnow;
% result_out(path);
%%
%