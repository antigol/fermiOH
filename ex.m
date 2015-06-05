r = linspace(0, 10, 1000);

ef = 25; % énergie de Fermi en unité de hbar omega et en ignorant le 3/2

num = fermiOH(ef, r); % somme sur les états d'énergie inférieur à ef
ana = approx(ef, r);  % approximation de densité locale

plot(r, num, r, ana);

% Conclusion personnelle :
% la différence de hauteur est due à l'approximation de comptage du nombre d'états occupés
% les vagues sont due au caractère non quantique de l'approximation
