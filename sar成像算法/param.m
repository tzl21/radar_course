classdef param
%%%     
        %padzero
        % Na_zero
        % Nr_zero
        % Veq
        % Veq_m
        % f0
        % theta_rc
        % Fr
        % Kr
        % Br
        % Ba
        % Fa
%%%
    properties
        padzero
        Na_zero
        Nr_zero
        Veq
        Veq_m
        f0
        R0
        theta_rc
        Fr
        Kr
        Br
        Ba
        Fa
    end
    methods
        function obj = param (a,b,c,d,e,f,r,g,h,i,j,k,l)   %构造函数，函数类名一致，完成类中变量的初始化
            obj.padzero = a;
            obj.Na_zero = b;
            obj.Nr_zero = c;
            obj.Veq = d;
            obj.Veq_m = e;
            obj.f0 = f;
            obj.R0 = r;
            obj.theta_rc = g;
            obj.Fr = h;
            obj.Kr = i;
            obj.Br = j;
            obj.Ba = k;
            obj.Fa = l;
        end
    end
end