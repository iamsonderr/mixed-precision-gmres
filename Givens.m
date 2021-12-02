function [Rm_bar,gm_bar] = Givens( Hm_bar,beta_e1 )
    [~,real_m] = size(Hm_bar);
    
    %Rotate Matrix need to recrate every iteration
    %R = eye(n,n);%Rotate Matrix
    for k = 1:real_m
        R = eye(real_m+1,real_m+1);% eye return a unit matrix
        down = (Hm_bar(k,k)^2+Hm_bar(k+1,k)^2)^(1/2);
        s = Hm_bar(k+1,k)/down;
        c = Hm_bar(k,k)/down;
        R(k:k+1,k:k+1) = [c,s;-s,c];
        Hm_bar = R*Hm_bar;
        beta_e1 = R*beta_e1;
    end
    Rm_bar = Hm_bar(1:real_m,1:real_m);
    gm_bar = beta_e1(1:real_m);
end