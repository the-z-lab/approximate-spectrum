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

% Set up figures
close all
h_ext = figure; hold on
h_cdf = figure; hold on
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
    
    % Exact extinction
    ext = 3/(4*pi)*omega.*imag(C(:,1)+C(:,5)+C(:,9))/3;
    ext_i = 3/(4*pi)*omega'.*imag(p);
    
    % Exact resonance frequencies
    [~,ind] = max(ext_i, [], 2);
    omega_res = omega(ind);
    
    % Use the lowest freq as the electrostatic dipole
    p_0 = real(p(:,1));
    
    % Approximate the resonant frequencies
    omega_res_approx = sqrt(4*pi./(3*p_0));

    % Approximate particle spectra
    h = 9/(gamma*(2+eps_inf)^2);
    W = gamma;
    ext_i_approx = h./(1+(omega_res_approx/W).^2.*(omega'./omega_res_approx - omega_res_approx./omega').^2);
    
    % Approximate average spectrum
    ext_approx = squeeze(mean(ext_i_approx));
    
    % Bin the resonant frequencies
    cdf = histcounts(omega_res, omega_edges, 'Normalization', 'cdf');
    cdf_approx = histcounts(omega_res_approx, omega_edges, 'Normalization', 'cdf');
    
    % Plot the exact and approximate extinction spectra
    figure(h_ext)
    plot(omega, ext, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, ext_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
    % Plot the exact and approximate resonant distributions
    figure(h_cdf)
    plot(omega, 1-cdf, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, 1-cdf_approx, '--', 'Color', colors(i,:), 'LineWidth', 1)
    
end

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

end