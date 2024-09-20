close all;
clear;
set(0, 'DefaultAxesFontSize', 11, 'DefaultLineLineWidth',2);
load('mos2_fit_data.mat');
% FileName = 'MoS2_1um_a_freq_eps1_eps2.txt';
% data = load(FileName);
wavelength = linspace(0.38, 0.68, 100);
freq= 1./wavelength;  %% input wavelength
% eps1 = data(:,2);
% eps2 = data(:,3);
% 
% eps_complex = eps1 + 1j .* eps2;

eps_opt = abs(x0(1)) ...
    + abs(x0(2)).^2.*abs(x0(3))./(abs(x0(2))^2 - freq.^2 - 1i*freq*x0(4)) +...
    + abs(x0(5)).^2.*abs(x0(6))./(abs(x0(5))^2 - freq.^2 - 1i*freq*x0(7))...
    + abs(x0(8)).^2.*abs(x0(9))./(abs(x0(8))^2 - freq.^2 - 1i*freq*x0(10))...
    + abs(x0(11)).^2.*abs(x0(12))./(abs(x0(11))^2 - freq.^2 - 1i*freq*x0(13));

n_fit = sqrt(eps_opt);
n_glass = ones(size(eps_opt))*1.4;
n_air = ones(size(eps_opt));
n_mat = vertcat(n_air,n_fit, n_glass);
% n_mat = vertcat(n_air,n_fit, n_air);

n_layer = 3;
theta = pi()/180*[0 20, 40, 60, 80];
% theta = pi()./180*[0];



r_s = zeros(length(theta), length(n_fit), n_layer-1);
t_s = zeros(length(theta), length(n_fit), n_layer-1);
r_p = zeros(length(theta), length(n_fit), n_layer-1);
t_p = zeros(length(theta), length(n_fit), n_layer-1);
M_s_11 = zeros(length(theta), length(n_fit), n_layer-1);
M_s_12 = zeros(length(theta), length(n_fit), n_layer-1);
M_s_21 = zeros(length(theta), length(n_fit), n_layer-1);
M_s_22 = zeros(length(theta), length(n_fit), n_layer-1);

M_p_11 = zeros(length(theta), length(n_fit), n_layer-1);
M_p_12 = zeros(length(theta), length(n_fit), n_layer-1);
M_p_21 = zeros(length(theta), length(n_fit), n_layer-1);
M_p_22 = zeros(length(theta), length(n_fit), n_layer-1);
P_11 = zeros(length(theta), length(n_fit), n_layer-2);
P_22 = zeros(length(theta), length(n_fit), n_layer-2);

thick = 0.61/1000;

k0 = 2*pi./wavelength;
for i = 1:length(theta)
    k_x = k0.*sin(theta(i)).*n_mat(1,:);
    for j = 1:n_layer-1
        k_j = k0.*n_mat(j,:);
        k_y_j = sqrt(k_j.^2 - k_x.^2);
        k_jp1 = k0.*n_mat(j+1,:);
        k_y_jp1 = sqrt(k_jp1.^2 - k_x.^2);
        r_s(i,:,j) = (k_y_jp1-k_y_j)./(k_y_jp1+k_y_j);
        t_s(i,:,j) = (2.*k_y_j)./(k_y_jp1+k_y_j);
        r_p(i,:,j) = -((k_y_jp1./k_jp1).*k_j - (k_y_j./k_j).*k_jp1)./((k_y_jp1./k_jp1).*k_j + (k_y_j./k_j).*k_jp1);
        t_p(i,:,j) = (2.*k_y_j)./((k_y_jp1./k_jp1).*k_j + (k_y_j./k_j).*k_jp1);

        M_s_11(i,:,j) = 1./t_s(i,:,j);
        M_s_12(i,:,j) = r_s(i,:,j)./t_s(i,:,j);
        M_s_21(i,:,j) = r_s(i,:,j)./t_s(i,:,j);
        M_s_22(i,:,j) = 1./t_s(i,:,j);

        M_p_11(i,:,j) = 1./t_p(i,:,j);
        M_p_12(i,:,j) = r_p(i,:,j)./t_p(i,:,j);
        M_p_21(i,:,j) = r_p(i,:,j)./t_p(i,:,j);
        M_p_22(i,:,j) = 1./t_p(i,:,j);

        if j~=1
            P_11(i,:,j-1) = exp(-1i.*k_y_j.*thick(j-1));
            P_22(i,:,j-1) = exp(+1i.*k_y_j.*thick(j-1));  
        end
    end
end


for j = 1:n_layer-1
    if j==1
        Q_p_11 = M_p_11(:,:,1);
        Q_p_12 = M_p_12(:,:,1);
        Q_p_21 = M_p_21(:,:,1);
        Q_p_22 = M_p_22(:,:,1);
        Q_s_11 = M_s_11(:,:,1);
        Q_s_12 = M_s_12(:,:,1);
        Q_s_21 = M_s_21(:,:,1);
        Q_s_22 = M_s_22(:,:,1);
    else
        % disp(j)
        
        Q_p_11 = Q_p_11.*P_11(:,:,j-1);
        Q_p_12 = Q_p_12.*P_22(:,:,j-1);
        Q_p_21 = Q_p_21.*P_11(:,:,j-1);
        Q_p_22 = Q_p_22.*P_22(:,:,j-1);
        Q_s_11 = Q_s_11.*P_11(:,:,j-1);
        Q_s_12 = Q_s_12.*P_22(:,:,j-1);
        Q_s_21 = Q_s_21.*P_11(:,:,j-1);
        Q_s_22 = Q_s_22.*P_22(:,:,j-1);

        
        Q_p_11 = Q_p_11.*M_p_11(:,:,j)+Q_p_12.*M_p_21(:,:,j);
        Q_p_12 = Q_p_11.*M_p_12(:,:,j)+Q_p_12.*M_p_22(:,:,j);
        Q_p_21 = Q_p_21.*M_p_11(:,:,j)+Q_p_22.*M_p_21(:,:,j);
        Q_p_22 = Q_p_21.*M_p_12(:,:,j)+Q_p_22.*M_p_22(:,:,j);
        Q_s_11 = Q_s_11.*M_s_11(:,:,j)+Q_s_12.*M_s_21(:,:,j);
        Q_s_12 = Q_s_11.*M_s_12(:,:,j)+Q_s_12.*M_s_22(:,:,j);
        Q_s_21 = Q_s_21.*M_s_11(:,:,j)+Q_s_22.*M_s_21(:,:,j);
        Q_s_22 = Q_s_21.*M_s_12(:,:,j)+Q_s_22.*M_s_22(:,:,j);
    end
end


T_s = abs((1./Q_s_11).^2);
R_s = abs((Q_s_21./Q_s_11).^2);


T_p = abs((1./Q_p_11).^2);
R_p = abs((Q_p_21./Q_p_11).^2);

% plot(wavelength, T_p(2,:));
%%
close all;
figure;  
plot(1./freq, 1-n_glass.*T_p(1,:)-R_p(1,:));
figure;

plot(1./freq, n_glass.*T_p(1,:));
figure;
plot(1./freq, R_p(1,:));




%%

close all;
for ang = round(theta(2:end)*180/pi)
    T_ref = load(strcat('MoS2_trans_Ez_', string(ang), '.txt'));
    abs_ref = load(strcat('MoS2_abs_Ez_', string(ang), '.txt'));
    R_ref = load(strcat('MoS2_refl_Ez_', string(ang), '.txt'));

    figure;
    plot(1./freq, T_s(round(ang/20)+1,:));
    hold on;
    plot(T_ref(:,1), T_ref(:,2),LineStyle='none', Marker='+');
    hold off;

    figure;
    plot(1./freq, 1- R_s(round(ang/20)+1,:)-T_s(round(ang/20)+1,:));
    hold on;
    plot(T_ref(:,1), abs_ref(:,2),LineStyle='none', Marker='+');
    hold off;
    % 
    % figure;
    % plot(1./freq, R_s(round(ang/20)+1,:));
    % hold on;
    % plot(T_ref(:,1), R_ref(:,2),LineStyle='none', Marker='+');
    % hold off;
    % 

end

distFig('Rows', 2, 'Columns',4);
%%
close all;
for ang = round(theta(2:end)*180/pi)
    T_ref = load(strcat('MoS2_trans_Hz_', string(ang), '.txt'));
    abs_ref = load(strcat('MoS2_abs_Hz_', string(ang), '.txt'));
    refl_ref = load(strcat('MoS2_refl_Hz_', string(ang), '.txt'));

    figure;
    plot(1./freq, T_p(round(ang/20)+1,:));
    hold on;
    plot(T_ref(:,1), T_ref(:,2),LineStyle='none', Marker='+');
    hold off;

    figure;
    plot(1./freq, 1- R_p(round(ang/20)+1,:)-T_p(round(ang/20)+1,:));
    hold on;
    plot(T_ref(:,1), abs_ref(:,2),LineStyle='none', Marker='+');
    hold off;
    % figure;
    % plot(1./freq, R_p(round(ang/20)+1,:));
    % hold on;
    % plot(T_ref(:,1), refl_ref(:,2),LineStyle='none', Marker='+');
    % hold off;

end

distFig('Rows', 2, 'Columns',4);


% 
%         r_s(i,j,:) = (k_y_jp1-k_y_j)/(k_y_jp1+k_y_j);
%         t_s(i,j,:) = (2*k_y_j)/(k_y_jp1+k_y_j);
%         r_p(i,j,:) = 
% 
% 
% phi = 2*pi*0.61./1000.*n_fit.*freq;
% 
% t_const = 1./( (n_fit+1).^2./(4.*n_fit).*exp(-1i.*phi) - (n_fit-1).^2./(4.*n_fit).*exp(+1i.*phi) );
% r_const = ( (n_fit.^2-1)./(4.*n_fit).*(-exp(-1i.*phi) + exp(+1i.*phi))  )./( (n_fit+1).^2./(4.*n_fit).*exp(-1i.*phi) - (n_fit-1).^2./(4.*n_fit).*exp(+1i.*phi) );
% 
% abs_ref = load('1D_MoS2_absorb_arr.txt');
% refl_ref = load("1D_MoS2_reflect_arr.txt");
% trans_ref = load("1D_MoS2_trans_arr.txt");
% 
% plot(1./freq, abs(r_const).^2)
% hold on;
% plot(abs_ref(:,1), refl_ref(:,2),LineStyle='--')
% hold off;
% 
% figure;
% plot(1./freq, abs(t_const).^2)
% hold on;
% plot(abs_ref(:,1), trans_ref(:,2),LineStyle='--')
% hold off;
% 
% 
% figure;
% plot(1./freq, 1- abs(t_const).^2- abs(r_const).^2);
% hold on;
% plot(abs_ref(:,1), abs_ref(:,2),LineStyle='--')
% hold off;
% 

