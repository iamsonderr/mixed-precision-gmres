function x = BackwardUpperTriangular( A,b )

    [~,A_size] = size(A);
    x = zeros(A_size,1);
    x(end) = b(end)/A(A_size,A_size);
    for k = A_size:-1:1
        x(k) = (b(k) - A(k,k+1:A_size) * x(k+1:A_size))/A(k,k);
    end
end