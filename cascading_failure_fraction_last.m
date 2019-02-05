function [N_failed,V_state,S_time] = cascading_failure_fraction_last(F0,Ax,Ay,Axy,Ayx,kx,ky,kxy,kyx,pmax,T) 
%****inputs***
%F0 is the set of initial failures that should be a vector(1*N)
%Ax is the adjacanecy matrix for the nodes in network X
%Ay is the adjacanecy matrix for the nodes in network Y
%Axy is the adjacancey matrix showing the connectivity of network X&Y(coulmn j shows the parents of node j of network Y)
%Ayx is the adjacancey matrix showing the connectivity of network X&Y(coulmn j shows the parents of node j of network X)
%kx and ky are respectively, are the minimum number of parents(intra-parents) that should be failed before node X or y can fail 
%kxy and kyx are respectively, are the minimum number of parents(inter-parents) that should be failed before node X or y can fail
%Kxy is used to calculate the probability of failure for the nodes of network Y based on their parents of network X
%pmax is the probability of a node being failed within one time step after all of its parents failed
%T is the time horizon that we want to see the behaviour of networks at the
%end of this
%****outputs
%N_failed is the number of failed nodes at the end of time T
%V_state is the vector of nodes state after time T
% if size(Axy)~=[length(Ax) length(Ay)] || size(Ayx)~=[length(Ay) length(Ax)] 
%     error('Invalid matrix size')
% end
% 
% if ~isvector(F0)
%     error('Input must be a vector')
% end

network_size=length(Ax)+length(Ay);
network_state=zeros(1,network_size);
for i=1:length(F0)
    f=F0(i);
    network_state(1,f)=1;  %we set the state of network according to the initial failure by putting one for the failed nodes
end

for j=1:T
    node_index=zeros(1,1);
    for h=1:network_size
        s1=network_state(1,h);
        sum_failed_parent_intra_X=0;
        sum_failed_parent_inter_YX=0;
        sum_failed_parent_intra_Y=0;
        sum_failed_parent_inter_XY=0;
        if s1==0
            if h<=length(Ax)
                parents_intra_h=Ax(:,h);
                parents_inter_h=Ayx(:,h);
                index_parents_intra_h=find(parents_intra_h==1);
                index_parents_inter_h=find(parents_inter_h==1);
                for i=1:length(index_parents_intra_h)
                    sum_failed_parent_intra_X=sum_failed_parent_intra_X+network_state(index_parents_intra_h(i));
                end
                for i=1:length(index_parents_inter_h)
                    sum_failed_parent_inter_YX=sum_failed_parent_inter_YX+network_state(index_parents_inter_h(i)+length(Ax));
                end
                
                number_parent_intra=length(index_parents_intra_h);
                number_parent_inter=length(index_parents_inter_h);
                
                if number_parent_intra~=0 && number_parent_inter~=0
                  kx1=ceil(kx*length(index_parents_intra_h));
                  kyx1=ceil(kyx*length(index_parents_inter_h));
                elseif number_parent_intra~=0 && number_parent_inter==0
                  kx1=ceil(kx*length(index_parents_intra_h));
                  kyx1=10000;
                elseif number_parent_intra==0 && number_parent_inter~=0
                  kx1=10000;
                  kyx1=ceil(kyx*length(index_parents_inter_h));  
                else
                   kx1=10000;
                   kyx1=10000;
                end
                if sum_failed_parent_intra_X>=kx1 && sum_failed_parent_inter_YX<kyx1
                    p_failed_intra_x=(sum_failed_parent_intra_X-kx1+1)*(pmax/(length(index_parents_intra_h)-kx1+1));
                    chance_number=rand(1,1); 
                    if chance_number<=p_failed_intra_x
                       node_index=[node_index h];
                    end
                elseif sum_failed_parent_intra_X<kx1 && sum_failed_parent_inter_YX>=kyx1
                    p_failed_inter_yx=(sum_failed_parent_inter_YX-kyx1+1)*(pmax/(length(index_parents_inter_h)-kyx1+1));
                    chance_number=rand(1,1); 
                    if chance_number<=p_failed_inter_yx
                        node_index=[node_index h];
                    end
                elseif sum_failed_parent_intra_X>=kx1 && sum_failed_parent_inter_YX>=kyx1
                     p_failed_inter_yx=(sum_failed_parent_inter_YX-kyx1+1)*(pmax/(length(index_parents_inter_h)-kyx1+1));
                     p_failed_intra_x=(sum_failed_parent_intra_X-kx1+1)*(pmax/(length(index_parents_intra_h)-kx1+1));
                     p_failed=(p_failed_inter_yx*(1-p_failed_intra_x))+(p_failed_intra_x*(1-p_failed_inter_yx))+(p_failed_inter_yx*p_failed_intra_x);
                     chance_number=rand(1,1); 
                    if chance_number<=p_failed
                        node_index=[node_index h];
                    end
                end
            elseif h>length(Ax)
                parents_intra_h=Ay(:,h-length(Ax));
                parents_inter_h=Axy(:,h-length(Ax));
                index_parents_intra_h=find(parents_intra_h==1);
                index_parents_inter_h=find(parents_inter_h==1);
                for i=1:length(index_parents_intra_h)
                    sum_failed_parent_intra_Y=sum_failed_parent_intra_Y+network_state(index_parents_intra_h(i)+length(Ax));
                end
                for i=1:length(index_parents_inter_h)
                    sum_failed_parent_inter_XY=sum_failed_parent_inter_XY+network_state(index_parents_inter_h(i));
                end
                number_parent_intra=length(index_parents_intra_h);
                number_parent_inter=length(index_parents_inter_h);
                if number_parent_intra~=0 && number_parent_inter~=0
                  ky1=ceil(ky*length(index_parents_intra_h));
                kxy1=ceil(kxy*length(index_parents_inter_h));
                elseif number_parent_intra~=0 && number_parent_inter==0
                  ky1=ceil(ky*length(index_parents_intra_h));
                  kxy1=10000;
                elseif number_parent_intra==0 && number_parent_inter~=0
                 ky1=10000;
                 kxy1=ceil(kxy*length(index_parents_inter_h));
                else
                   ky1=10000;
                   kxy1=10000;
                end

                
                    if sum_failed_parent_intra_Y>=ky1 && sum_failed_parent_inter_XY<kxy1
                    p_failed_intra_y=(sum_failed_parent_intra_Y-ky1+1)*(pmax/(length(index_parents_intra_h)-ky1+1));
                    chance_number=rand(1,1); 
                    if chance_number<=p_failed_intra_y
                        node_index=[node_index h];
                    end
                    elseif sum_failed_parent_intra_Y<ky1 && sum_failed_parent_inter_XY>=kxy1
                    p_failed_inter_xy=(sum_failed_parent_inter_XY-kxy1+1)*(pmax/(length(index_parents_inter_h)-kxy1+1));
                    chance_number=rand(1,1); 
                    if chance_number<=p_failed_inter_xy
                        node_index=[node_index h];
                    end
                    elseif sum_failed_parent_intra_Y>=ky1 && sum_failed_parent_inter_XY>=kxy1
                        p_failed_inter_xy=(sum_failed_parent_inter_XY-kxy1+1)*(pmax/(length(index_parents_inter_h)-kxy1+1));
                        p_failed_intra_y=(sum_failed_parent_intra_Y-ky1+1)*(pmax/(length(index_parents_intra_h)-ky1+1));
                        p_failed=(p_failed_inter_xy*(1-p_failed_intra_y))+(p_failed_intra_y*(1-p_failed_inter_xy))+(p_failed_inter_xy*p_failed_intra_y);
                        chance_number=rand(1,1); 
                        if chance_number<=p_failed
                        node_index=[node_index h];
                        end
                    end
            end
        end
        end
        
    node_index=node_index(1,2:length(node_index));
    for d=1:length(node_index)
        w=node_index(d);
        network_state(1,w)=1;
    end
   state_over_time{j}=network_state; 
end

N_failed=length(find(network_state==1));
V_state=network_state;
S_time=state_over_time;    


end


