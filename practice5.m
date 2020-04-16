%%
% practice 5
% solve Lambert probrem
% �P�ʂ�km,s
% Reference:
% ���g�r�Y, "�~�b�V������͂ƋO���݌v�̊�b", pp.74-76.
clear;
%%
R1 = [  5000 10000 2100];
R2 = [-14600  2500 7000];
DELTA_T = 36000;
%����
NN = 0;

%�Œ�l
k       = [0,0,1];      %z�����P�ʃx�N�g��
MU_Ea_KM3S2 = 398600;   %�n���̏d�͒萔[km^3/s^2]

[V1, V2] = calc_Lam(R1, R2, DELTA_T, NN, k, MU_Ea_KM3S2);
%%
% Lambert problem
function [V1, V2] = calc_Lam(R1, R2, DELTA_T, NN, k, MU_Ea_KM3S2)
    % calc
    r1 = norm(R1);          %�o���n�_�ł̈ʒu�̃m����
    r2 = norm(R2);          %�����n�_�ł̈ʒu�̃m����

    %�J�ڊp�̌v�Z(���W�A���o��)
    if 0 <= dot(cross(R1, R2), k)
        theta = acos(dot(R1, R2) / (r1 * r2));
    else
        theta = 2 * pi - acos(dot(R1, R2) / (r1 * r2));
    end
    if theta < 0
        fprintf('ERROR: ��<-2��\n');
    end
    theta = theta + 2 * pi * NN;

    %�e�l�̌v�Z
    c       = sqrt(r1 ^ 2 + r2 ^ 2 - 2 * r1 * r2 * cos(theta));
    a_m     = (r1 + r2 + c) / 4;
    s       = 2 * a_m;
    T_m     = 2 * pi * sqrt(a_m ^ 3 / MU_Ea_KM3S2);
    beta_m  = 2 * asin(sqrt((s - c) / s));   % 0 <= beta <= pi

    %�O���`��̑I�� <<1:�ȉ~�C2:�������C3:�o�Ȑ�>>
    if (2 * pi * NN <= theta) && (theta <= pi * (2 * NN + 1))
        delta_tp1 = (1 / 3) * sqrt(2 / MU_Ea_KM3S2) * (s ^ 1.5 - (s - c) ^ 1.5);
        if DELTA_T > delta_tp1
            trajectory = 1;
            %�ȉ~�^�C�v�̏ꍇ�����̂��߂Ƀ�t�����߂�D
            delta_tm1 = NN * T_m + sqrt(a_m ^ 3 / MU_Ea_KM3S2) * (pi - (beta_m - sin(beta_m)));
        elseif DELTA_T == delta_tp1
            trajectory = 2;
        elseif DELTA_T < delta_tp1
            trajectory = 3;
        else
            fprintf('ERROR:�������O���v�Z��C0�`pi�̎��̑J�ڋO�����ނɎ��s\n');
        end
    elseif (pi * (2 * NN + 1) < theta) && (theta <= 2 * pi * (NN + 1))
        delta_tp2 = (1 / 3) * sqrt(2 / MU_Ea_KM3S2) * (s ^ 1.5 + (s - c) ^ 1.5);
        if DELTA_T > delta_tp2
            trajectory = 1;
            %�ȉ~�^�C�v�̏ꍇ�����̂��߂Ƀ�t�����߂�D
            delta_tm2 = NN * T_m + sqrt(a_m ^ 3 / MU_Ea_KM3S2) * (pi + (beta_m - sin(beta_m)));
        elseif DELTA_T == delta_tp2
            trajectory = 2;
        elseif DELTA_T < delta_tp2
            trajectory = 3;
        else
            fprintf('ERROR:�������O���v�Z��Cpi�`2pi�̎��̑J�ڋO�����ނɎ��s\n');
        end
    else
        fprintf('ERROR:�������O���̌v�Z�Ɏ��s\n');
    end
    %%
    if NN == 0
        [a, p] = Calc_Newton0(a_m, MU_Ea_KM3S2, s, c, trajectory, NN, theta, DELTA_T, delta_tm1, r1, r2);
        [V1, V2] = Calc_vel(MU_Ea_KM3S2, p, R1, R2, theta);
        fprintf('Transfer Angle[deg]: %.2f\n',theta/pi*180);
        fprintf('Velocity Dep.[km/s]: %.2f  %.2f  %.2f,\n',V1);
        fprintf('Velocity Arr.[km/s]: %.2f %.2f %.2f,\n',V2);
    else
        trajectory = 1;
        if T_m * NN > DELTA_T
            fprintf("�𖳂�\n");
        else
            delta_tm1 = NN * T_m + sqrt(a_m ^ 3 / MU_Ea_KM3S2) * (pi - (beta_m - sin(beta_m)));
            delta_tm2 = NN * T_m + sqrt(a_m ^ 3 / MU_Ea_KM3S2) * (pi + (beta_m - sin(beta_m)));
            if (2 * pi * NN < theta) && (theta < pi * (2 * NN + 1))
                if DELTA_T <= delta_tm1
                    calc_ellipse = 1;
                else
                    calc_ellipse = 2;
                end
            else
                if DELTA_T <= delta_tm2
                    calc_ellipse = 3;
                else
                    calc_ellipse = 4;
                end
            end
            %�j���[�g�����t�\���@
            n = 50;         %�J��Ԃ��v�Z�̏��
            if (calc_ellipse == 2)||(calc_ellipse == 4)
                m = 2;
            else 
                m = 1;
            end
            for j = 1 : m
                a = a_m * 1.1;    %�����l
                if (calc_ellipse == 2) && (j == 2)
                    calc_ellipse = 1;
                elseif (calc_ellipse == 4) && (j == 2)
                    calc_ellipse = 3;
                end
                for i = 1 : n
                    %�p�����[�^�̍Ē�`
                    T     = 2 * pi * sqrt(a ^ 3 / MU_Ea_KM3S2);
                    alpha = 2 * asin(sqrt(s / (2 * a)));
                    beta  = 2 * asin(sqrt((s - c) / (2 * a)));
                    gamma = 2 * asinh(sqrt(s / (2 * abs(a))));
                    delta = 2 * asinh(sqrt((s - c) / (2 * abs(a))));
                    switch(calc_ellipse)
                        case 1
                            delta_t2     = NN * T + sqrt(a ^ 3 / MU_Ea_KM3S2) * ((alpha - sin(alpha)) - (beta - sin(beta)));     %(4.3.8)
                            delta_t2_dot = 3 * delta_t2 / (2 * a) - 1 / sqrt(MU_Ea_KM3S2 * (a ^ 3)) * ((s ^ 2) / sin(alpha) - ((s - c) ^ 2) / sin(beta));
                        case 2
                            delta_t2     = (NN + 1) * T -sqrt(a ^ 3 / MU_Ea_KM3S2) * ((alpha - sin(alpha)) + (beta - sin(beta))); %(4.3.9)
                            delta_t2_dot = 3*delta_t2/(2*a)+1/sqrt(MU_Ea_KM3S2*(a^3))*((s^2)/sin(alpha)+((s-c)^2)/sin(beta));
                        case 3
                            delta_t2     = NN * T + sqrt(a ^ 3 / MU_Ea_KM3S2) * ((alpha - sin(alpha)) + (beta - sin(beta)));     %(4.3.10)
                            delta_t2_dot = 3 * delta_t2 / (2 * a) - 1 / sqrt(MU_Ea_KM3S2 * (a ^ 3)) * ((s ^ 2) / sin(alpha) + ((s - c) ^ 2) / sin(beta));
                        case 4
                            delta_t2     = (NN + 1) * T - sqrt(a ^ 3 / MU_Ea_KM3S2) * ((alpha - sin(alpha)) - (beta - sin(beta))); %(4.3.11)
                            delta_t2_dot = 3 * delta_t2 / (2 * a) + 1 / sqrt(MU_Ea_KM3S2 * (a ^ 3)) * ((s ^ 2) / sin(alpha) + ((s - c) ^ 2) / sin(beta));
                    end
                    if (abs(DELTA_T - delta_t2) > DELTA_T * (10 ^ (-5))) && (i < n)
                        a = a + (DELTA_T - delta_t2) / delta_t2_dot;
                    else
                        break
                    end
                end
                %�e�픻��
                if (abs(DELTA_T - delta_t2) > DELTA_T * (10 ^ (-5))) && (i == n)
                    fprintf('��������\n');
                    continue
                end
                if a < a_m
                    fprintf('a���͈͊O\n');
                    continue
                end
                if (alpha < 0) || (pi < alpha)
                    fprintf('alpha���͈͊O\n');
                    continue
                end
                if (beta < 0) || (pi < beta)
                    fprintf('beta���͈͊O\n');
                    continue
                end
                if (abs(DELTA_T - delta_t2) < DELTA_T * 10 ^ (-5)) && (a > 0)
                    %�O���v�f�̌v�Z
                    s = 2 * a * ((sin(alpha / 2)) ^ 2);
                    c = 2 * a * sin((alpha + beta) / 2) * sin((alpha - beta) / 2);
                    switch(trajectory)
                        case 1
                            switch(calc_ellipse)
                                case 1
                                    p = (4 * a * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((alpha + beta) / 2)) ^ 2); 
                                case 2
                                    p = (4 * a * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((alpha - beta) / 2)) ^ 2);
                                case 4
                                    p = (4 * a * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((alpha + beta) / 2)) ^ 2); 
                                case 3
                                    p = (4 * a * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((alpha - beta) / 2)) ^ 2);
                            end
                        case 3
                            switch(calc_hyperbola)
                                case 1
                                    p = (4 * abs(a) * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((gamma + delta) / 2)) ^ 2);
                                case 2
                                    p = (4 * abs(a) * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((gamma - delta) / 2)) ^ 2);
                            end
                    end
                    [V1, V2] = Calc_vel(MU_Ea_KM3S2, p, R1, R2, theta);
                    fprintf('Transfer Angle[deg]: %.2f\n',theta / pi * 180);
                    fprintf('Velocity Dep.[km/s]: %.2f  %.2f  %.2f,\n', V1);
                    fprintf('Velocity Arr.[km/s]: %.2f %.2f %.2f,\n', V2);
                end
            end
        end
    end
end
%%
% �j���[�g�����t�\���@
function [a, p] = Calc_Newton0(a_m, mu, s, c, trajectory, NN, theta, DELTA_T, delta_tm1, r1, r2)
    n = 50;         %�J��Ԃ��v�Z�̏��
    a = a_m * 1.1;    %�����l
    for i = 1 : n
        %�p�����[�^�̍Čv�Z
        T     = 2 * pi * sqrt(a ^ 3 / mu);
        alpha = 2 * asin(sqrt(s / (2 * a)));              % 0<= alpha <=pi
        beta  = 2 * asin(sqrt((s - c) / (2 * a)));          % 0<= beta  <=pi
        gamma = 2 * asinh(sqrt(s / (2 * abs(a))));
        delta = 2 * asinh(sqrt((s - c) / (2 * abs(a))));
        switch(trajectory)
            case 1
                if (2 * pi * NN <= theta) && (theta <= pi * (2 * NN + 1))
                    if DELTA_T < delta_tm1
                        delta_t2     = NN * T + sqrt(a ^ 3 / mu) * ((alpha - sin(alpha)) - (beta - sin(beta)));     %(4.3.8)
                        delta_t2_dot = 3 * delta_t2 / (2 * a) - 1 / sqrt(mu * (a ^ 3)) * ((s ^ 2) / sin(alpha) - ((s - c) ^ 2) / sin(beta));
                        calc_ellipse = 1;
                    elseif DELTA_T > delta_tm1
                        delta_t2     = (NN + 1) * T - sqrt(a ^ 3 / mu) * ((alpha - sin(alpha)) + (beta - sin(beta))); %(4.3.9)
                        delta_t2_dot = 3 * delta_t2 / (2 * a) +...
                                       1 / sqrt(mu * (a ^ 3)) * ((s ^ 2) / sin(alpha) + ((s - c) ^ 2) / sin(beta));
                        calc_ellipse = 2;
                    else
                        fprintf('ERROR:case1��delta_T1��delta_tm1�̒l�̔�r�Ɏ��s(�������Ȃ�)(i=%d)\n',i);
                        break
                    end
                elseif (pi * (2 * NN + 1) < theta) && (theta <= 2 * pi * (NN + 1))
                    if DELTA_T < delta_tm2
                        delta_t2     = NN * T + sqrt(a ^ 3 / mu) * ((alpha - sin(alpha)) + (beta - sin(beta)));     %(4.3.10)
                        delta_t2_dot = 3 * delta_t2 / (2 * a) - 1 / sqrt(mu * (a ^ 3)) * ((s ^ 2) / sin(alpha) + ((s - c) ^ 2) / sin(beta));
                        calc_ellipse = 2;
                    elseif DELTA_T > delta_tm2
                        delta_t2 = (NN + 1) * T - sqrt(a ^ 3 / mu) * ((alpha - sin(alpha)) - (beta - sin(beta))); %(4.3.11)
                        delta_t2_dot = 3 * delta_t2 / (2 * a) + 1 / sqrt(mu * (a ^ 3)) * ((s ^ 2) / sin(alpha) + ((s - c) ^ 2) / sin(beta));
                        calc_ellipse = 1;
                    else
                        fprintf('ERROR:case1��delta_T1��delta_tm2�̒l�̔�r�Ɏ��s(�������Ȃ�)(i=%d)\n',i);
                        break
                    end
                else
                    fprintf('ERROR:case1�ŏꍇ�������s(i=%d)\n',i);
                    break
                end

            case 2
                fprintf('�������O���ƂȂ�(i=%d)\n',i);
                break

            case 3
                if (0 <= theta) && (theta <= pi)
                    delta_t2       = sqrt((abs(a)) ^ 3 / mu) * ((sinh(gamma) - gamma) - (sinh(delta) - delta));     %(4.3.12)
                    delta_t2_dot   = 3 * delta_t2 / (2 * abs(a)) - 1 / sqrt(mu * (abs(a) ^ 3)) * (s ^ 2 / sinh(gamma) - ((s - c) ^ 2) / sinh(delta));
                    calc_hyperbola = 1;
                elseif (pi < theta) && (theta <= 2*pi)
                    delta_t2       = sqrt((abs(a)) ^ 3 / mu) * ((sinh(gamma) - gamma) + (sinh(delta) - delta));     %(4.3.13)
                    delta_t2_dot   = 3 * delta_t2 / (2 * abs(a)) - 1 / sqrt(mu * (abs(a) ^ 3)) * (s ^ 2 / sinh(gamma) + ((s - c) ^ 2) / sinh(delta));
                    calc_hyperbola = 2;
                else
                    fprintf('ERROR:case3�ŏꍇ�������s(i=%d)\n',i);
                    break
                end
        end
        if (abs(DELTA_T - delta_t2) > DELTA_T * (10 ^ (-5))) && (i < n)
            a = a + (DELTA_T - delta_t2) / delta_t2_dot;
        elseif (abs(DELTA_T - delta_t2) > DELTA_T * (10 ^ (-5))) && (i == n)
            fprintf('��������(i=%d)\n',i);
            break
        else
            break
        end
    end

    %���Z�o��̏���
    if abs(DELTA_T - delta_t2) < DELTA_T * 10 ^ (-5)
        %�O���v�f
        s = 2 * a * ((sin(alpha / 2)) ^ 2);
        c = 2 * a * sin((alpha + beta) / 2) * sin((alpha - beta) / 2);
        switch(trajectory)
            case 1
                switch(calc_ellipse)
                    case 1
                        p = (4 * a * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((alpha + beta) / 2)) ^ 2); 
                    case 2
                        p = (4 * a * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((alpha - beta) / 2)) ^ 2);
                end
            case 3
                switch(calc_hyperbola)
                    case 1
                        p = (4 * abs(a) * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((gamma + delta) / 2)) ^ 2);
                    case 2
                        p = (4 * abs(a) * (s - r1) * (s - r2)) / (c ^ 2) * ((sin((gamma - delta) / 2)) ^ 2);
                end
        end
    end
end
%%
% calculate velocity
function [V1, V2] = Calc_vel(mu, p, R1, R2, theta)
    r1 = norm(R1);
    r2 = norm(R2);
    %���O�����W���̌W���̌v�Z
    g_bar = sqrt(mu * p) / (r1 * r2 * sin(theta));
    f     = 1 - r2 / p * (1 - cos(theta));
    g_dot = 1 - r1 / p * (1 - cos(theta));
    %���x�x�N�g���̌v�Z
    V1 = g_bar * (R2 - f * R1);
    V2 = g_bar * (g_dot * R2-R1);
end