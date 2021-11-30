function [T,bk] = Givens( H,b )
    [~,realm] = size(H);
    
    %Rotate Matrix need to recrate every iteration
    %R = eye(n,n);%Rotate Matrix
    for k = 1:realm
        R = eye(realm+1,realm+1);% eye return a unit matrix
        down = (H(k,k)^2+H(k+1,k)^2)^(1/2);
        s = H(k+1,k)/down;
        c = H(k,k)/down;
        R(k:k+1,k:k+1) = [c,s;-s,c];
        H = R*H;
        b = R*b;
    end
    T = H(1:realm,1:realm);
    bk = b(1:realm);
end