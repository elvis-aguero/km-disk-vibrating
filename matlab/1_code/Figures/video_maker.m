%clc
clear vidObj
close all;
cd ..
addpath(fullfile(pwd, "simulation_code" ));

p = uigetdir();
%p = fullfile(pwd, "..",  "D50Quant100\rho1000sigma7220nu98muair0\RhoS1000SigmaS7220\R0350mm\ImpDefCornerAng180U39\N=20tol=5.00e-05");
cd(p);

global errored
errored = ~isfile('z.mat');
try
    load('ProblemConditions.mat'); NN = N;
catch
    load('U0.mat');
    load('Fr.mat');
    disp("Couldn't find Problem Conditions");
end
load_vars('vz.mat'); Vo = abs(vz(1));


try
    load_vars('etas.mat');
    etaMatPer = etas;
catch
    files = dir(fullfile(pwd, "etaMatPer*.mat"));
    N = length(files);
    etaAux = [];
    for i = 1:N
        load_vars(files(i).name);
        etaAux = [etaAux, etaMatPer];
    end
    etaMatPer = etaAux;
end
load_vars('z.mat')
load_vars('etaOri.mat')
load_vars('tvec.mat')

% Loading oscillating conditions. CHANGED
resFile = dir("simulationResults*.mat");
if ~isempty(resFile)
    load(resFile(1).name, 'PROBLEM_CONSTANTS');
    Gamma = PROBLEM_CONSTANTS.bath_forcing_amplitude;
    Fr = PROBLEM_CONSTANTS.froude;
    w = PROBLEM_CONSTANTS.force_frequency;
else
    % Fallback to existing or default values
    Gamma = 0; 
    try load('Fr.mat','Fr'); catch; Fr = 1; end
    w = 1; 
end

zb = Gamma/(Fr*w^2)*cos(w*tvec); % Elevation of the pool in lab frame. CHANGED
zbplot = zb;

load_vars('oscillation_amplitudes.mat');

h = figure();
h.Position = [100 100 1260 400];
for ii = floor(linspace(1, size(etaMatPer,2), 300))

    %% Drop video
    subplot(1, 3, 2);
    thetaplot = linspace(0, pi, 100);
    zsTop = zs_from_spherical(thetaplot, oscillation_amplitudes(:, ii));
    xsTop = r_from_spherical(thetaplot, oscillation_amplitudes(:, ii)); 
    plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop] + z(ii) + zb(ii),'k','Linewidth',2); % CHANGED: added zb to show lab frame
    hold on

    plot([fliplr(-1*r),r(2:end)],(zbplot(ii)+[0;flipud(etaMatPer(:,ii));etaMatPer(2:end,ii);0]),...
        'color',[.4 .4 .4],'LineWidth',2)
    hold on
    %axis equal
    grid on
    set(gca,'xlim',[-2.5 2.5],'ylim',[-2 3],'Xtick',-5:5,'FontName','Times','FontSize',14);
    xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',14)
    ylabel('$z/R_o$','interpreter','Latex','FontName','Times',...
        'FontSize',20,'rotation',90)
%     t = (ii-1)/360;
%     to = floor(t);
%     t1 = floor(10*(t-to));
%     t2 = round(100*(t-to-t1/10));
%     if t2 == 10
%         t2=0;
%         t1=t1+1;
%     end
%     if t1 == 10
%         t1 = 0;
%         to = to+1;
%     end
    hold off
    title(sprintf("$ t/t_\\sigma =\\ $ %3.2f", tvec(ii)),'FontSize',18,...
            'interpreter','latex','FontName','Times')   
        
    %% Pressure field distribution
    
    f = zeta_generator(pressure_amplitudes(:, ii));
    pfield = f(thetaplot) - sum(pressure_amplitudes(:, ii));
    pmean = mean(pfield(1:50));
    
    %position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
    %    1/3 * subplotHeight + (1-subplotHeight)/2, subplotWidths/subplots, subplotHeight/3];
    %subplot('Position', position);
    subplot(1, 3, 1);
    plot(thetaplot*180/pi, pfield-pmean);
    set(gca,'xlim',[0 180], 'ylim', [-0.5, 2], 'Xtick', [45, 90, 135], ...
        'Ytick', [0 1 2], 'FontName','Times','FontSize',14);
    grid on
    %axis equal
    %xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',18)
    xlabel('$\theta$','interpreter','Latex','FontName','Times','FontSize',14); 
    ylabel('$p(\theta)$','interpreter','Latex','FontName','Times',...
        'FontSize',20,'rotation',90)
    
    
    a = gca;
    a.XRuler.TickLabelGapOffset = -4;

    %% Pressure amplitudes
    
    %position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
    %    (1-subplotHeight)/2-0.04, subplotWidths/subplots, subplotHeight/3];
    subplot(1, 3, 3);
    bar(0:NN, [-sum(pressure_amplitudes(:, ii)); pressure_amplitudes(:, ii)]);
    set(gca, 'ylim', [-0.7, 0.7], 'Xtick', [0, floor(NN/2)], 'Ytick', [-.5 0 .5], ...
        'FontName','Times','FontSize',14);
    grid on
    xlabel('$\ell$','interpreter','Latex','FontName','Times','FontSize',14); 
    ylabel('$B_{\ell}$','interpreter','Latex','FontName','Times',...
        'FontSize',20,'rotation',90)

    a = gca;
    a.XRuler.TickLabelGapOffset = -4;
        
    drawnow
    %g = gcf;
    %g.WindowState = 'maximized';
    currFrame = getframe(h);
    writeVideo(vidObj,currFrame);
    hold off
end
close(vidObj);


function load_vars(str)
    global errored
    
    if errored == true
        str = "errored_" + str; 
    end
    
    vars = load(str);
    fn = fieldnames(vars);
    for ii = 1:length(fn)
        assignin('caller', fn{ii}, vars.(fn{ii}));
    end
    
end



