function invmat=invert_matrix_with_lu(L,U,P);

% function invmat=invert_matrix_with_lu(L,U);
% Calculate matrix inverse for given L and U from LU decomposition
% 
% (c) M. Huang, Ph.D.
% Based on Numerical Recipes in C

N=size(L,1);
A=L+U-eye(N,N); % pack L and U in the A matrix, minus the trival unit diag elements
index=zeros(N,1);
for i=1:N       % now compress the permutation matrix into a vector
    index(i)=find(P(i,:)==1);
end 
for k=1:N
    b=zeros(N,1);
    b(k)=1.0;
    ii=0;
    for i=1:N
        ip=i;
        sum=b(ip);
        b(ip)=b(i);
        if ii,
            for j=ii:(i-1)
                sum=sum-A(i,j)*b(j);
            end
        elseif sum,
            ii=i;
        end
        b(i)=sum;
    end    
    for i=N:(-1):1
        sum=b(i);
        for j=(i+1):N
            sum=sum-A(i,j)*b(j);
        end
        b(i)=sum/A(i,i);
    end
    invmat(:,k)=b; 
end
invmat(:,index)=invmat;
return

