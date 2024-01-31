clear all;

M = 5;
hm = [1:M];
hm = sort(hm,'descend');
R=2;
Pij = zeros(M,M);
Pij(1,1) = (2^R-1)/hm(1);
for m =2 : M
    bmi=[];
    for i = 1 : m
        bmi(i) = 1/(1+hm(m)*sum(Pij(i:m-1,i)));
    end

    A = []; % No other constraints
    b = [];
    Aeq = [];%
    beq = [];%zeros(M+1,1); 
    lb = [];
    ub = [];
    x0 =  ones(m,1); %initilization 
    options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
    x = fmincon(@(x) sum(x),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,m,bmi,hm,R),options);
    Pij(m,1:m) = x';

end

Pija = zeros(M,M);
Pija(1,1) = (2^R-1)/hm(1);
for m = 2 : M
    for i = 1 : m
        prod1 = 1;
        for p =1 : m
            prod1 = prod1 * hm(m)/(hm(m)*sum(Pija(p:m-1,p))+1);
        end

        Pija(m,i) = (exp(R)/prod1)^(1/m) - (hm(m)*sum(Pija(i:m-1,i)) +1 )/hm(m);


    end
end

Iij = zeros(M,M);
Iij(1,1) = Pija(1,1);
for m = 2 : M
    for i = 1 : m
        Iij(m,i) = sum(Pija(1:m,i));

    end
end
%test
% a = hm(2)/(hm(2)*Pij(1,1)+1);
% lambda = sqrt(a*hm(2)/exp(R));
% P21 = 1/lambda - 1/a
% P22 = 1/lambda - 1/hm(2)
interf = diag(hm)* Pij;

Eall = sum(sum(Pija))

for m = 1 : M
    Poma(m) = (exp(R)-1)/hm(m);
end
Ealloma = sum(Poma)

function [c,ceq] = mycons(x,m,bmi,hm,R)

sum1 = 0;
for i =1 : m
    sum1 = sum1 + log(1+bmi(i)*hm(m)*x(i));
end

c(1,1) = R-sum1;
c(2:m+1,1) = -x;
ceq = [];
end
