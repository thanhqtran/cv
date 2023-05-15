classdef Model < handle

    properties
        % =========================================
        % GLOBAL PARAMETERS SETTING
        % =========================================
        T = 30;
        alpha = 0.3;
        delta = 0.5; % research productivity shifter
        phi = 0.7; % past inventions' share
        lambda = 0.49; % researcher's spillover effects
        rho = 0.075; % time-cost of children
        e_bar = 0.012; % default education
        n_diff = (1.96 + 0.86) / (1.96 - 0.86); % difference in birth rates
        muL = 3; % skill conversion efficiency (low)
        muH = 5; % skill conversion efficiency (high)
        epsilon = 0.01; % good cost of child rearing

        % Experiment parameters
        theta;
        nu;
        n;

        % Initial conditions
        NL;
        NH;
        N;
        A;
        K;
        w_thes;

        % Simulation results
        time;
        x; %variety increase (triangle A(t-1))
        R;
        wL;
        wH;
        eL;
        nL;
        eH;
        nH;
        cL;
        cH;
        sL;
        sH;
        dL;
        dH;
        uL;
        uH;
        Y;
        Xalpha;
        L;
        H;
        HR;
        HY;
        PiL;
        PiH;
        xi;
        % SOME INDICATORS
        skill_prem;
        pop_ratio_H;
        pop_ratio_L;
        pop_ratio;
        labor_ratio;
        HY_ratio;
        HR_ratio;
        H_NH_ratio;
        gA;
        piA;
        k;
        kss;
        J;
        k_tilde;
        y;
        gy;
        gK;
        Pi_gap;
        pop_growth;
        sav_rate;
        edu_rate;
        fertility;
        total_pop;
        pop_growth2;
        % error checks 
        error_eq_KEY;
        error_A;
        error_N;
        %error_rc;

        % additional params property 
        beta;
        gamma;
        eta;
        LRTFR;
    end

    methods
        % Constructor
        function obj = Model(theta, nu, n, LRTFR)
            % Set experiment parameters
            obj.theta = theta;
            obj.nu = nu;
            obj.n = n;
            obj.LRTFR = LRTFR; %reproduction rate

            % Set up some parameters 
            obj.beta = 0.99 ^ (obj.T * 4); % consumption weight
            obj.eta = 1 - 1 / obj.n_diff; % returns to education
            obj.gamma = ((1 + obj.beta) * obj.rho * obj.LRTFR) / ((1 - obj.eta - obj.LRTFR * obj.rho)) %altruistic weight

            % Simulation results
            obj.time= [];
            obj.x = [];
            obj.error_eq_KEY = [];
            obj.R = [];
            obj.wL = [];
            obj.wH = [];
            obj.eL = [];
            obj.nL = [];
            obj.eH = [];
            obj.nH = [];
            obj.cL = [];
            obj.cH = [];
            obj.sL = [];
            obj.sH = [];
            obj.uL = [];
            obj.uH = [];
            obj.dL = [];
            obj.dH = [];
            obj.Y = [];
            obj.Xalpha = [];
            obj.L = [];
            obj.H = [];
            obj.HR = [];
            obj.HY = [];
            obj.PiL = [];
            obj.PiH = [];
            obj.xi = [];
            obj.piA = [];
            % SOME INDICATORS
            obj.skill_prem = [];
            obj.pop_ratio_H = [];
            obj.pop_ratio_L = [];
            obj.pop_ratio = [];
            obj.HY_ratio = [];
            obj.HR_ratio = [];
            obj.H_NH_ratio = [];
            obj.gA = [];
            obj.k = [];
            obj.J = [];
            obj.k_tilde = [];
            obj.y = [];
            obj.gy = [];
            obj.gK = [];
            obj.Pi_gap = [];
            obj.pop_growth = [];
            obj.pop_growth2 = [];
            obj.sav_rate = [];
            obj.edu_rate = [];
            obj.fertility = [];
            obj.total_pop = [];
            obj.labor_ratio = [];
            obj.kss = [];
            % error checks
            obj.error_A = [];
            obj.error_N = [];
            %obj.error_rc = [];
        end

        % Run Simulation
        function obj = run(obj)
            % at time zero (t=1)
            obj.time(1) = 0;
            obj.NL(1) = 2;
            obj.NH(1) = 1;
            obj.N(1) = obj.NL(1) + obj.NH(1);
            obj.A(1) = 1;
            obj.K(1) = 0.01;
            obj.y(1) = 0.12;
            obj.uL(1) = 0;
            obj.cL(1) = 0
            obj.uH(1) = 0;
            obj.cH(1) = 0; 
            % at the beginning of time 1 (t=2)
            obj.NL(2) = 2;
            obj.NH(2) = 1;
            obj.N(2) = obj.NL(2) + obj.NH(2);
            obj.K(2) = 0.01;
            obj.uL(2) = 0;
            obj.uH(2) = 0;
            obj.cL(2) = 0;
            obj.cH(2) = 0;
            % theshold level of wages
            obj.w_thes = (obj.e_bar - obj.eta * obj.epsilon) / (obj.eta * obj.rho);
            % loop
            for t = 2:obj.n
                syms x_var;
                assume(x_var > 0);
                obj.time(t) = t-1;
                obj.total_pop(t) = obj.N(t) + obj.N(t-1);

                % Solve for Variety Increase x
                obj.xi = (1 - obj.eta) * obj.rho * obj.gamma / (1 + obj.beta + obj.gamma);
                eq_KEY = obj.delta * (obj.A(t - 1) ^ obj.phi) * (obj.NH(t) * (1 - obj.xi / (obj.rho + ((obj.epsilon - obj.e_bar) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha))) * ((obj.A(t - 1) + x_var * (obj.alpha ^ (1 / (1 - obj.alpha)))) ^ obj.alpha)) / ((1 - obj.alpha) * (obj.K(t) ^ obj.alpha) * (obj.A(t - 1) ^ (obj.alpha * obj.phi)) * ((obj.A(t - 1) + x_var * (obj.alpha ^ (obj.alpha / (1 - obj.alpha)))) ^ (1 - obj.alpha))))) - ((obj.alpha ^ (1 / (1 - obj.alpha))) * (obj.A(t - 1) + x_var * (obj.alpha ^ (obj.alpha / (1 - obj.alpha))))) / (obj.delta * (obj.A(t - 1) ^ obj.phi))) - x_var == 0;
                obj.x(t) = vpasolve(eq_KEY, x_var);
                %% test if solution is true
                eq_KEY = obj.delta * (obj.A(t - 1) ^ obj.phi) * (obj.NH(t) * (1 - obj.xi / (obj.rho + ((obj.epsilon - obj.e_bar) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha))) * ((obj.A(t - 1) + obj.x(t) * (obj.alpha ^ (1 / (1 - obj.alpha)))) ^ obj.alpha)) / ((1 - obj.alpha) * (obj.K(t) ^ obj.alpha) * (obj.A(t - 1) ^ (obj.alpha * obj.phi)) * ((obj.A(t - 1) + obj.x(t) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha)))) ^ (1 - obj.alpha))))) - ((obj.alpha ^ (1 / (1 - obj.alpha))) * (obj.A(t - 1) + obj.x(t) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha))))) / (obj.delta * (obj.A(t - 1) ^ obj.phi))) - obj.x(t); %must be close to zero for the solution to be correct
                
                if obj.x(t) <= 0
                    obj.x(t) = 0;
                else
                    obj.x(t) = obj.x(t);
                end

                % SOLVE FOR A
                obj.A(t) = obj.A(t - 1) + obj.x(t);
                obj.HY(t) = (obj.alpha ^ (1 / (1 - obj.alpha))) * (obj.A(t - 1) + obj.x(t) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha)))) / (obj.delta * (obj.A(t - 1) ^ obj.phi));

                % DERIVE THE Factor Prices
                obj.R(t) = obj.alpha * (obj.K(t) / (obj.HY(t) * (obj.A(t - 1) + obj.x(t) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha)))))) ^ (obj.alpha - 1);
                obj.wL(t) = obj.theta * (obj.A(t - 1) ^ obj.nu);
                obj.wH(t) = (1 - obj.alpha) * ((obj.alpha / obj.R(t)) ^ ((obj.alpha) / (1 - obj.alpha))) * (obj.A(t - 1) + obj.x(t) * obj.alpha ^ ((obj.alpha) / (1 - obj.alpha)));

                % HOUSEHOLD EDU AND FERTI CHOICES
                %% Low-skilled parents
                if (obj.wL(t) <= obj.w_thes)
                    obj.eL(t) = 0;
                    obj.nL(t) = obj.gamma * obj.wL(t) / ((1 + obj.beta + obj.gamma) * (obj.rho * obj.wL(t) + obj.epsilon));
                else
                    obj.eL(t) = (obj.eta * (obj.rho * obj.wL(t) + obj.epsilon) - obj.e_bar) / (1 - obj.eta);
                    obj.nL(t) = (1 - obj.eta) * obj.gamma * obj.wL(t) / ((1 + obj.beta + obj.gamma) * (obj.rho * obj.wL(t) + obj.epsilon - obj.e_bar));
                end

                %% High-skilled parents
                if (obj.wH(t) <= obj.w_thes)
                    obj.eH(t) = 0;
                    obj.nH(t) = obj.gamma * obj.wL(t) / ((1 + obj.beta + obj.gamma) * (obj.rho * obj.wL(t) + obj.epsilon));
                else
                    obj.eH(t) = (obj.eta * (obj.rho * obj.wL(t) + obj.epsilon) - obj.e_bar) / (1 - obj.eta);
                    obj.nH(t) = (1 - obj.eta) * obj.gamma * obj.wL(t) / ((1 + obj.beta + obj.gamma) * (obj.rho * obj.wL(t) + obj.epsilon - obj.e_bar));
                end

                % HOUSEHOLD CONSUMPTION
                obj.cL(t) = obj.wL(t) / (1 + obj.beta + obj.gamma);
                obj.cH(t) = obj.wH(t) / (1 + obj.beta + obj.gamma);
                obj.sL(t) = obj.beta * obj.wL(t) / (1 + obj.beta + obj.gamma);
                obj.sH(t) = obj.beta * obj.wH(t) / (1 + obj.beta + obj.gamma);
                obj.dL(t) = obj.beta * obj.R(t) * obj.wL(t) / (1 + obj.beta + obj.gamma);
                obj.dH(t) = obj.beta * obj.R(t) * obj.wH(t) / (1 + obj.beta + obj.gamma);
                obj.PiL(t) = obj.muL * ((obj.eL(t) + obj.e_bar) / obj.A(t)) ^ obj.eta;
                obj.PiH(t) = obj.muH * ((obj.eH(t) + obj.e_bar) / obj.A(t)) ^ obj.eta;
                obj.uL(t-1) = log(obj.cL(t-1)) + obj.beta*(log(obj.dL(t))) + obj.gamma*(log(obj.nL(t)*obj.PiL(t)));
                obj.uH(t-1) = log(obj.cH(t-1)) + obj.beta*(log(obj.dH(t))) + obj.gamma*(log(obj.nH(t)*obj.PiH(t)));
                % LABOR EQUILIBRIUM
                obj.L(t) = (1 - obj.rho * obj.nL(t)) * obj.NL(t);
                obj.H(t) = (1 - obj.rho * obj.nH(t)) * obj.NH(t);

                % R&D SECTOR
                obj.HR(t) = obj.H(t) - obj.HY(t);
                obj.error_A = obj.delta * obj.A(t - 1) ^ obj.phi * obj.HR(t) - obj.x(t); %must be zero for the solution to be correct

                % SECTORAL PRODUCTION
                obj.Xalpha(t) = (obj.HY(t) ^ obj.alpha) * ((obj.alpha / obj.R(t)) ^ ((obj.alpha) / (1 - obj.alpha))) * (obj.A(t - 1) + obj.x(t) * obj.alpha ^ ((obj.alpha) / (1 - obj.alpha)));
                obj.Y(t) = obj.A(t - 1) * obj.L(t) + obj.Xalpha(t) * (obj.H(t) ^ (1 - obj.alpha));

                % SOME INDICATORS
                %% skill premium
                obj.skill_prem(t) = obj.wH(t) / obj.wL(t);
                obj.pop_ratio_H(t) = obj.NH(t) / (obj.NL(t) + obj.NH(t));
                obj.pop_ratio_L(t) = obj.NL(t) / (obj.NL(t) + obj.NH(t));
                obj.labor_ratio(t) = obj.L(t) / obj.H(t);
                obj.pop_ratio(t) = obj.NL(t) / obj.NH(t);
                obj.HY_ratio(t) = obj.HY(t) / obj.H(t);
                obj.HR_ratio(t) = obj.HR(t) / obj.H(t);
                obj.H_NH_ratio(t) = obj.H(t) / obj.NH(t);
                obj.gA(t) = obj.x(t) / obj.A(t);
                obj.k(t) = obj.K(t) / obj.N(t);
                obj.J(t) = obj.A(t - 1) + obj.x(t) * (obj.alpha ^ (obj.alpha / (1 - obj.alpha)));
                obj.k_tilde(t) = obj.K(t) / (obj.A(t) * obj.L(t) + obj.A(t) * obj.H(t));
                obj.y(t) = obj.Y(t) / obj.N(t);
                obj.gy(t) = obj.y(t) / obj.y(t - 1);
                obj.Pi_gap(t) = obj.PiH(t)/obj.PiL(t);
                obj.pop_growth(t) = obj.N(t) / obj.N(t - 1);
                obj.pop_growth2(t) = obj.N(t) / obj.N(t - 1) - 1;
                obj.gK(t) = obj.K(t) / obj.K(t - 1);
                obj.sav_rate(t) = (obj.sL(t) * obj.NL(t) + obj.sH(t) * obj.NH(t)) / (obj.wL(t) * obj.NL(t) + obj.wH(t) * obj.NH(t)); 
                obj.edu_rate(t) = (obj.eL(t) * obj.NL(t) + obj.eH(t) * obj.NH(t)) / (obj.wL(t) * obj.NL(t) + obj.wH(t) * obj.NH(t));
                obj.fertility(t) = obj.nL(t)/obj.NL(t) + obj.nH(t)/obj.NH(t);
                % UPDATE VALUES FOR THE NEXT PERIOD
                %% population
                obj.N(t + 1) = obj.NL(t) * obj.nL(t) + obj.NH(t) * obj.nH(t);
                obj.NL(t + 1) = obj.NL(t) * obj.nL(t) * (1 - obj.PiL(t)) + obj.NH(t) * obj.nH(t) * (1 - obj.PiH(t));
                obj.NH(t + 1) = obj.NL(t) * obj.nL(t) * obj.PiL(t) + obj.NH(t) * obj.nH(t) * obj.PiH(t);
                obj.error_N = obj.N(t + 1) - obj.NL(t + 1) - obj.NH(t + 1) %must be zero for the solution to be correct
                %% capital
                obj.K(t + 1) = obj.sL(t) * obj.NL(t) + obj.sH(t) * obj.NH(t);
                %innovator's profit 
                obj.piA(t) = obj.x(t) * (1-obj.alpha)* (obj.alpha^((1+obj.alpha)/(1-obj.alpha)))*(obj.R(t)^(-obj.alpha/(1-obj.alpha)))*(obj.HY(t));
                % SOME VERIFICATION CHECKS
                %% resource constraint
                %obj.error_rc(t) = obj.Y(t) - obj.wL(t) * obj.L(t) - obj.wH(t) * obj.H(t) - obj.dL(t) * obj.NL(t - 1) - obj.dH(t) * obj.NH(t - 1);
            end
        end
    end
end
