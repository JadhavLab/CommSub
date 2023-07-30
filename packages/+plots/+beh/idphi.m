function plot_idphi(behavior)
    scatter3(behavior.time, behavior.lindist, idphi(behavior.X, behavior.Y, 8)./8, 'CData', behavior.time);
    xlabel('time');
    ylabel('lindist');
    zlabel('idphi\nVicariousTrialAndError');
    alpha(0.1);
    colorbar;

