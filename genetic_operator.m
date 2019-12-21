function f  = genetic_operator(parent_chromosome, M, V)

%% function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
%
% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed.
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
% mu - distribution index for crossover (read the enlcosed pdf file)
% mum - distribution index for mutation (read the enclosed pdf file)
% l_limit - a vector of lower limit for the corresponding decsion variables
% u_limit - a vector of upper limit for the corresponding decsion variables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector.

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

[N,m] = size(parent_chromosome);

clear m

child=[];
p = 1;
% Flags used to set if crossover and mutation were actually performed.
was_crossover = 0;
was_mutation = 0;

global distance;

for i = 1 : N
    % With 90 % probability perform crossover
    if rand(1) < 0.5
        % Initialize the children to be null vector.
        child_1 = [];
        child_2 = [];
        % Select the first parent
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        % Select the second parent
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        % Make sure both the parents are not the same.
        pcount=0;
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))&&pcount<100
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
            pcount=pcount+1;
        end
        if pcount<100
            
        else
            continue;
        end
        % Get the chromosome information for each randomnly selected
        % parents
        parent_1 = parent_chromosome(parent_1,1:M+V);
        parent_2 = parent_chromosome(parent_2,1:M+V);
        
        % perform crossover
        % get inner nodes of parent_1, parent_2
        inode1=parent_1(parent_1(1:V)~=0);
        inode1=inode1(2:end-1);% inner nodes of parent_1
        inode2=parent_2(parent_2(1:V)~=0);
        inode2=inode2(2:end-1);% inner nodes of parent_2
        % checked matrix, if not check path between 2 nodes, then 0; otherwise 1
        ichk=zeros(size(inode1,2),size(inode2,2));
        % if 2 inner nodes have same ID, then mark as checked
        internode=intersect(inode1,inode2);
        for loop=1:length(internode)
            ichk(find(internode(loop)==inode1),find(internode(loop)==inode2))=1;
            ichk(find(internode(loop)==inode2),find(internode(loop)==inode1))=1;
        end
        
        % global distance
        idistance=distance(inode1,inode2);
        W=1./idistance; % weight
        
        crvflag=1;
        while ~isempty(ichk==0)&&crvflag
            [Indr,Indc]=find(rand<W,1);
            Ind1=inode1(Indr);
            Ind2=inode2(Indc);
            
            % if select none,
            if isempty(Ind1)||isempty(Ind2)
                continue;
            end
            crov1=find(parent_1==Ind1,1);
            crov2=find(parent_2==Ind2,1);
            
            child_1t=[parent_1(1:crov1),parent_1(crov1+1:M+V)];
            child_2t=[parent_2(1:crov2),parent_2(crov2+1:M+V)];
            
            % mark (Ind1,Ind2) as checked
            ichk(Ind1,Ind2)=1;
            ichk(Ind2,Ind1)=1;
            
            % check if child is valid
            [c1data,vc1]=result_out(child_1t(child_1t(1:V)~=0));
            [c2data,vc2]=result_out(child_2t(child_2t(1:V)~=0));
            if ~vc1
                child_1=child_1t;
            end
            if ~vc2
                child_2=child_2t;
            end
            
            if vc1&&vc2
                % both are invlaid
            else
                crvflag=0;%
            end
        end
        
        % Set the crossover flag. When crossover is performed two children
        % are generate, while when mutation is performed only only child is
        % generated.
        was_crossover = 1;
        was_mutation = 0;
        % With 10 % probability perform mutation. Mutation is based on
        % polynomial mutation.
    else
        % Select at random the parent.
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        % Get the chromosome information for the randomnly selected parent.
        child_3 = parent_chromosome(parent_3,1:M+V);
        % Perform mutation of the selected parent
        % child_3 's path
        child_3p=child_3(child_3(1:V)~=0);
        % checked matrix, if one inner node checked, then mark =1; otherwise =0
        ichk=zeros(length(child_3p)-2,1);
        
        mflag=1;% mutation flag
        % Choose delete node or change node Randomly
        if rand(1)<0.2
            % change node
            while ~isempty(ichk==0)&&mflag
                Ind_3=find(ichk==0);
                Ind_3=Ind_3(randi(length(Ind_3)))+1;
                
                last_mndata=result_out(child_3p(1:Ind_3-1));% airInfo at mutation node's parent
                % valid nodes list for mutation node's parent
                vlist_last_mn=se_vnode(last_mndata(end,:),child_3p(1:Ind_3-1));
                %             list_last_mpt=list_last_mp(:,1);
                
                % delete current mutation node
                if ~isempty(vlist_last_mn)
                    vlist_last_mn=vlist_last_mn(vlist_last_mn(:,1)~=child_3p(Ind_3),:);
                end
                
                while ~isempty(vlist_last_mn)&&mflag
                    Ind_mnt=randi(length(vlist_last_mn(:,1)));% Index of mutation node tmp
                    child_3t=child_3;
                    child_3t(Ind_3)=vlist_last_mn(Ind_mnt);% tmp child_3
                    [child_3td,mflag]=result_out(child_3t(child_3t(1:V)~=0));
                    if mflag
                        vlist_last_mn(Ind_mnt,:)=[];
                    else
                        child_3=child_3t;
                        % Set the mutation flag
                        was_mutation = 1;
                        was_crossover = 0;
                    end
                end
                
                % mark current mutation node as checked
                ichk(Ind_3)=1;
                
                %                 if ~mflag
                %                     child_3=child_3t;
                %                     % Set the mutation flag
                %                     was_mutation = 1;
                %                     was_crossover = 0;
                %                 end
            end
        else
            % delete one node
            while isempty(ichk==0)&&mflag
                Ind_3=find(ichk==0);
                Ind_3=Ind_3(randi(length(Ind_3)))+1;
                child_3tp=[child_3p(1:Ind_3-1),child_3p(Ind_3+1:end)];
                [child_3td,mflag]=result_out(child_3tp);
                if mflag
                    ichk(Ind_3+1)=1;
                else
                    child_3(1:length(child_3p))=[child_3tp,0];
                    % Set the mutation flag
                    was_mutation=1;
                    was_crossover=0;
                end
                
            end
            
        end
        
    end
    
    % Keep proper count and appropriately fill the child variable with all
    % the generated children for the particular generation.
    if was_crossover
        if ~vc1
            child(p,:)=child_1;
            p=p+1;
        end
        if ~vc2
            child(p,:)=child_2;
            p=p+1;
        end
        was_crossover = 0;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
    
end

% evaluate obj function, concatenate chromosome with obj value
for loop=1:size(child,1)
    child_chtmp=child(loop,1:V);
    childdata=result_out(child_chtmp(child_chtmp~=0));
    child(loop,V+1:end)=childdata(end,end-1:end);
end

f = child;
