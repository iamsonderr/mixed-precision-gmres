function x = BackUT( A,b )
%backward:   ?????????
% ??Ax=b
    [~,c] = size(A);
    A = A(1:c,:);
    x = zeros(c,1);
    x(end) = b(end)/A(c,c);
    for k = 2:c
        V = x(c-k+2:c);
        x(c-k+1) = (b(c-k+1)-A(c-k+1,c-k+2:end)*V)/A(c-k+1,c-k+1);
    end
end