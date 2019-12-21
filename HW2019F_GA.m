%% Setup the Import Options

clear
%% 读取源数据

choice=input('Choose 1 or 2 problem to be solved:');

global gf;
% f1=table2array(tab_f1);
% save('f1table.mat','f1')
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
%% 参数

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

% matrix: distance between two node
global distance;
distance=zeros(nn,nn);
for loop=1:nn
    for iloop=loop+1:nn
        distance(loop,iloop)=norm(gf(loop,2:4)-gf(iloop,2:4),2);
    end
end
distance=distance+transpose(distance);
%% 画出三维图

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
%% Prepare

% ParaP=gcp;

fprintf("All parameter loaded, start GA Algorithm...\n");
global fpath;
fpath=plot3(gf(row_nan,2),gf(row_nan,3),gf(row_nan,4),'k');% figure path: plot path in figure
%% NSGA-II

popu=round(nn/10);% population
gens=900;% generations

% Init chromosome
M=2; % number of obj functions
V=nn;% number of vars

pool=round(popu/2);
tour=2;

chromosome=chromo_init(popu,M,V);
chromosome=non_domination_sort_mod(chromosome,M,V);
fprintf("Chromosome Init Finished.\n");

for loop=1:gens
    % Selection
    parent_chromosome=tournament_selection(chromosome,pool,tour);
    %     parent_chromosome=parent_chromosome(:,1:M+V);
    
    % Perform crossover and mutation
    offspring_chromosome=genetic_operator(parent_chromosome,M,V);
    
    % intermediate_chromosome is a concatenation of current population and
    % offspring population
    % TODO: There are some intermediate_chromosomes same with each other, need
    % to IDENTIFY them, and kick duplicated ones out!
%     intermediate_chromosome=vertcat(chromosome(:,1:M+V),offspring_chromosome);
    intermediate_chromosome=union(chromosome(:,1:M+V),offspring_chromosome,'rows');
    
    % Non-domination-sort
    intermediate_chromosome=non_domination_sort_mod(intermediate_chromosome,M,V);
    
    % Perform selection
    chromosome=replace_chromosome(intermediate_chromosome,M,V,popu);
    
    % User defined operation: Figure
    for iloop=1:size(chromosome,1)
        tmp_path=chromosome(iloop,1:V);
%         tmp_path=chromosome(1,1:V);
        tmp_path=tmp_path(tmp_path~=0);
        set(fpath,'XData',gf(tmp_path,2),'YData',gf(tmp_path,3),'ZData',gf(tmp_path,4));
        drawnow;
    end
    fprintf("#%d generation finished...\n",loop);
    pause(0.1);
    
end

%% Result

fprintf("GA Algorithm end!\n");

for loop=1:popu
    
end
path1=chromosome(1,1:nn);
path1=path1(path1>0);
set(fpath,'XData',gf(path1,2),'YData',gf(path1,3),'ZData',gf(path1,4));
drawnow;
% plot3(gf(path1,2),gf(path1,3),gf(path1,4),'k');

%% Clear
% delete('ParaP');

%%
