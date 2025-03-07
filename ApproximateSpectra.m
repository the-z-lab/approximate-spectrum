function [] = ApproximateSpectra()

% Compare the approximate extinction spectrum from a single electrostatics
% measurement to the exact spectrum. Only works with eps_inf = 1 for now.

% Specify variable parameters
N = 1000; % number of particles
eta = [0.01, 0.05, 0.10:0.10:0.50]'; % colloid volume fraction
gamma = 0.10; % damping (relative to k_p)
eps_inf = 1; % high-freq dielectric constant
path = sprintf('/Users/zacharysherman/Documents/Scattering/HS/MAT_Files/N_%.0f/', N); % path to data

% Initialize
N_eta = length(eta);
omega_edges = (0.01:0.01:0.81)' - 0.005;
addpath('/Users/zacharysherman/Documents/Scattering/Stratified_Media')

% Set up figures
close all
h_Cr = figure; hold on
h_ext = figure; hold on
h_cdf = figure; hold on
h_n_re = figure; hold on
h_n_im = figure; hold on
h_R = figure; hold on
h_T = figure; hold on
colors = parula(length(eta));

% Loop through volume fractions
for i = 1:N_eta
    
    % Load dipole data
    filename = sprintf('%seta%.2f_gamma%.2f_eps%d.mat', path, eta(i), gamma, eps_inf);
    data = load(filename, 'k', 'C', 'p');
    omega = data.k;
    C = data.C; % dims: freq., component-polarization
    p = data.p; % dims: particle, component-polarization, freq.
    
    % Extract only dipole components in the field direction and combine polarizations
    p = squeeze([p(:,1,:); p(:,5,:); p(:,9,:)]); % dims: particle, frequency
    
    % Exact capacitance
    C = (C(:,1)+C(:,5)+C(:,9))/3;
    C_re = real(C);
    C_im = imag(C);
    
    % Exact extinction
    ext = 3/(4*pi)*omega.*C_im;
    ext_i = 3/(4*pi)*omega'.*imag(p);
    
    % Exact resonance frequencies
    [~,ind] = max(ext_i, [], 2);
    omega_res = omega(ind);
    
    % Use the lowest freq as the electrostatic dipole
    p_0 = real(p(:,1));
    
    % Approximate the resonant frequencies
    omega_res_approx = sqrt(4*pi./(3*p_0));
    %omega_res_approx = omega_res;

    % Approximate particle spectra
    h = 9/(gamma*(2+eps_inf)^2);
    W = gamma;
    ext_i_approx = h./(1+(omega_res_approx/W).^2.*(omega'./omega_res_approx - omega_res_approx./omega').^2);
    C_re_i_approx = p_0.*(1 - (omega'./omega_res_approx).^2)./((omega'./omega_res_approx).^4 + ((W./omega_res_approx).^2 - 2).*(omega'./omega_res_approx).^2 + 1);
    C_im_i_approx = 4*pi/3*ext_i_approx./omega';
    
    % Approximate average spectrum
    ext_approx = squeeze(mean(ext_i_approx)).';
    C_re_approx = squeeze(mean(C_re_i_approx)).';
    C_im_approx = squeeze(mean(C_im_i_approx)).';
    C_approx = C_re_approx + 1i*C_im_approx;
    
    % Bin the resonant frequencies
    cdf = histcounts(omega_res, omega_edges, 'Normalization', 'cdf');
    cdf_approx = histcounts(omega_res_approx, omega_edges, 'Normalization', 'cdf');
    
    % Effective permittivity and refractive index
    eps_eff = 1 + 3*eta(i)/(4*pi)*C;
    eps_eff_approx = 1 + 3*eta(i)/(4*pi)*C_approx;
    n_eff = sqrt(eps_eff);
    n_eff_approx = sqrt(eps_eff_approx);
    
    % Transmittance and reflectance: exact
    h = [-inf, 0.1, inf];
    n_in = repmat([ones(length(omega),1), n_eff, ones(length(omega),1)], 1, 1, 3); % assume isotropic; replicate N_omega-by-3 (layers) into N_omega-by-3-by-3 (nx, ny, nz)
    [E_r_s, E_t_s, E_r_p, E_t_p] = Stratified_Media_Anisotropic(omega, 0, n_in, h);
    
    T_s = abs(E_t_s).^2;
    T_p = abs(E_t_p).^2;
    R_s = abs(E_r_s).^2;
    R_p = abs(E_r_p).^2;
    
    % Transmittance and reflectance: approx
    n_in = repmat([ones(length(omega),1), n_eff_approx, ones(length(omega),1)], 1, 1, 3); % assume isotropic; replicate N_omega-by-3 (layers) into N_omega-by-3-by-3 (nx, ny, nz)
    [E_r_s, E_t_s, E_r_p, E_t_p] = Stratified_Media_Anisotropic(omega, 0, n_in, h);
    
    T_s_approx = abs(E_t_s).^2;
    T_p_approx = abs(E_t_p).^2;
    R_s_approx = abs(E_r_s).^2;
    R_p_approx = abs(E_r_p).^2;
    
    % Plot the exact and approximate C' spectra
    figure(h_Cr)
    plot(omega, C_re, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, C_re_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate extinction spectra
    figure(h_ext)
    plot(omega, ext, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, ext_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate resonant distributions
    figure(h_cdf)
    plot(omega, 1-cdf, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, 1-cdf_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate refractive index
    figure(h_n_re)
    plot(omega, real(n_eff), 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, real(n_eff_approx), '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate extinction coefficient
    figure(h_n_im)
    plot(omega, imag(n_eff), 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, imag(n_eff_approx), '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate reflectance
    figure(h_R)
    plot(omega, R_p, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, R_p_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate transmittance
    figure(h_T)
    plot(omega, T_p, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, T_p_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
end

% Format real capacitance plot
figure(h_Cr)
xlabel('\omega')
ylabel('C''')
xlim([0.2, 0.8])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

% Format extinction plot
figure(h_ext)
xlabel('\omega')
ylabel('\sigma')
xlim([0.2, 0.8])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

% Format cdf plot
figure(h_cdf)
xlabel('\omega')
ylabel('1-cdf')
xlim([0.3, 0.7])
ylim([0, 1])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

% Format refractive index plot
figure(h_n_re)
xlabel('\omega')
ylabel('n''')
xlim([0.2, 0.8])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

% Format extinction coeff plot
figure(h_n_im)
xlabel('\omega')
ylabel('n''''')
xlim([0.2, 0.8])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

% Format reflectance plot
figure(h_R)
xlabel('\omega')
ylabel('reflectance')
xlim([0.2, 0.8])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

% Format transmittance plot
figure(h_T)
xlabel('\omega')
ylabel('transmittance')
xlim([0.2, 0.8])
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1],'XDir','reverse')

end