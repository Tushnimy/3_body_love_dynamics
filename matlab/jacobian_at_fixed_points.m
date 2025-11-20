function J = jacobian_at_fixed_point(x, params)
    syms x1 x2 x3 x4 x5 x6 x7 x8 real
    f1 = params(1)+(x3^2-x4^2)+params(4)*(x1+x4);
    f2 = 2*x3*x4+params(4)*(x2-x3);
    f3 = params(2)+(x1^2-x2^2)+params(5)*(x3+x2) +params(7)*x5;
    f4 = 2*x1*x2+params(5)*(x4-x1)+params(7)*x6;
    f5 = params(1)+(x7^2-x8^2)+params(4)*(x5+x8);
    f6 = 2*x7*x8+params(4)*(x6-x7);
    f7 = params(3)+(x5^2-x6^2)+params(6)*(x7+x6)+params(8)*x1;
    f8 = 2*x5*x6+params(6)*(x8-x5)+params(8)*x2;

    f = [f1; f2; f3; f4; f5; f6; f7; f8];

    J_sym = jacobian(f, [x1 x2 x3 x4 x5 x6 x7 x8]);
    J = double(subs(J_sym, [x1 x2 x3 x4 x5 x6 x7 x8], x));
end
