function [By_mean, By_del] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, n_per_group, verbose)
[Bx,By,Bz] = MT.Field2D(FOV_1,FOV_2,FOV_3,surface,n_per_group);

% Only keep the circular region
for idx_1=1:length(FOV_1)
    for idx_2=1:length(FOV_2)
        r = sqrt(FOV_1(idx_1)^2 + FOV_2(idx_2)^2);
        if r > half_FOV
            Bx(idx_2,idx_1) = NaN;
            By(idx_2,idx_1) = NaN;
            Bz(idx_2,idx_1) = NaN;
        end
    end
end

Bx = Bx*1e3;
By = By*1e3;
Bz = Bz*1e3;

%     font_size = 18;
%     figure;
%     pcolor(FOV_1,FOV_2,By);
%     axis square;
%     shading flat;
%     cb = colorbar;%('Ticks',[-21,-18,-15,-12,-9,-6]);
%     set(get(cb,'Title'),'string','mT');
%     cb.FontSize = font_size;
%     colormap jet;
%     xlabel('x/mm','FontSize',font_size);
%     if surface == 'y'
%         ylabel('z/mm','FontSize',font_size);
%         ttl = sprintf('By at y = %.1f mm',FOV_3);
%     else
%         ylabel('y/mm','FontSize',font_size);
%         ttl = sprintf('By at z = %.1f mm',FOV_3);
%     end
%     title(ttl,'FontSize',font_size);
%     %caxis([-21,-6]);
%     set(gcf,'color','w');
%     ax = gca;
%     ax.FontSize = font_size;

% Homogeneity & linearity
By_mean = mean(By(~isnan(By)));
By_del = range(By(~isnan(By)));
if verbose
    fprintf('Bandwidth: %.3f%%\n', By_del/By_mean*100);
    fprintf('Mean: %.1f mT\n', By_mean);
end
%     weight = Halbach_MT.Weight();
%     fprintf('Weight: %.1f kg\n', weight);
%     force = max(forces);

end