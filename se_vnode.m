function [vlist]=se_vnode(air_info,p)
% search available nodes for node air_info
% p is nodesNo already passed
global nn;
nnv=ones([nn,1]);
nnv(p)=0;
vlist=[];
for lp=1:nn
    if isPositiveIntegerValuedNumeric(nnv(lp))
        [ee,m]=check_err(air_info,lp);
        if ~ee
            vlist=vertcat(vlist,m);
        end
    end
end
end