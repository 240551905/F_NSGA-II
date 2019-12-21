function [pdata,t] = result_out(path)
%RESULT_OUT Summary of this function goes here
%   Detailed explanation goes here
l=length(path);
pdata=zeros(l,5);
pdata(:,1)=path;
for loop=2:l
    [t,D]=check_err(pdata(loop-1,:),pdata(loop,1));
    if t
%         fprintf("invalid path!\n");
        return;
    end
    pdata(loop,:)=D;
end

for loop=1:l
    fprintf("Node %d: horizon error: %f, vertical error: %f, long: %f, correctness counts: %f.\n",pdata(loop,1),pdata(loop,2),pdata(loop,3),pdata(loop,4),pdata(loop,5));
end

end

