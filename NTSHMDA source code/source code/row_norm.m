function [result]=row_norm(M)
[n1 n2]=size(M);
result=M;
for i=1:n1
    if sum(M(i,:))==0
        result(i,:)=M(i,:);
    else
        result(i,:)=M(i,:)/sum(M(i,:));
    end
end