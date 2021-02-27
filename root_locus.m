function [rc,departure_angles,arrival_angles,break_away_in,asymptotes,borderline_stability] = root_locus(sys)
% Plotting function for all elements of the root locus
[r,dep_an,arr_an,break_points,asymptotes,borderline_stability] = root_locus_calculation(sys);

legend_plots = [];
legend_descriptions = {};
% % Plotting
% Plotting via rlocus
figure;
rlocus(sys)
% Plotting Departure Angles
hold on
axis equal


radius = 0.2;
for i=1:length(dep_an)
    angle_pr = dep_an(i)/360;
    angfrac = 2*pi* angle_pr;
    angv = linspace(0, angfrac);
    x = radius*cos(angv);
    y = radius*sin(angv);
    x0 = real(r(i,1));
    y0 = imag(r(i,1));
    plot(x+x0, y+y0, '-r', 'LineWidth',1);
    
    plot([0 0; 1.1 1.1]*x(1)+x0, [0 0; 1.1 1.1]*y(1)+y0, [0 0; 1.1 1.1]*x(end)+x0, [0 0; 1.1 1.1]*y(end)+y0, '-r', 'LineWidth',1);
    
    
    text(x(length(x)/2)+x0, y(length(y)/2)+y0, sprintf('%.1f\\circ', angfrac*180/pi));   
end

% Plotting Arrival Angles
index = 0;
for i=1:length(arr_an)
    while(1/r(i+index,end)==0)        
        index = index + 1;
    end
    
    angle_pr = arr_an(i)/360;
    angfrac = 2*pi* angle_pr;
    angv = linspace(0, angfrac);
    x = radius*cos(angv);
    y = radius*sin(angv);
    x0 = real(r(i+index,end));
    y0 = imag(r(i+index,end));
    plot(x+x0, y+y0, '-m', 'LineWidth',1);
    plot([0 0; 1.1 1.1]*x(1)+x0, [0 0; 1.1 1.1]*y(1)+y0, [0 0; 1.1 1.1]*x(end)+x0, [0 0; 1.1 1.1]*y(end)+y0, '-m', 'LineWidth',1);
    text(x(length(x)/2)+x0, y(length(y)/2)+y0, sprintf('%.1f\\circ', angfrac*180/pi));   
end


% Plotting Breakaway/in points
if not(isempty(break_points))
    p_breakin = zeros(1,length(break_points));
    break_desc = {};
    for i=1:length(break_points)
        p_breakin(i) = plot(real(break_points(i)),imag(break_points(i)),'ks');
        x_temp = real(break_points(i));
        y_temp = imag(break_points(i));
        if x_temp<0
            sign_x = '-';
        else
            sign_x = '+';
        end
        if y_temp<0
            sign_y = '-';
        else
            sign_y = '+';
        end
        break_desc = [break_desc  strcat('Breakpoint:',sign_x,num2str(round(abs(x_temp),4)),sign_y,num2str(round(abs(y_temp),4)),'i')];
    end
    legend_plots = [legend_plots p_breakin];
    legend_descriptions = [legend_descriptions [break_desc]];
end

% Plotting asymptotes
syms x;
lamda = cell2mat(asymptotes(3));
asym_point = cell2mat(asymptotes(1));
if not(isempty(lamda))
    for i=1:length(lamda)
        if lamda(i) ~= 0
            if 1/lamda(i) == 0
                xx = [asym_point,asym_point];
                yy = [-100,100];
                plot(xx,yy,'k-.');
                p_asymp = plot(asym_point,0,'k-.');
                desc_asymp = strcat('Asymptote at x=',num2str(round(asym_point,4)));
            else
                yi = lamda(i)*x - lamda(i)*asym_point;
                fplot(yi,'k-.');
                p_asymp = plot(asym_point,0,'k-.');
                desc_asymp = strcat('Asymptote at x=',num2str(round(asym_point,4)));
            end
        else
            p_asymp = plot(asym_point,0,'k.');
            desc_asymp = 'Asymptote: Real axis.';
        end
    end
    legend_plots = [legend_plots p_asymp];
    legend_descriptions = [legend_descriptions [desc_asymp]];
end


% Plotting intersection points with y'y

desc_inter = [];
k_critical = cell2mat(borderline_stability(2));
roots_inter = cell2mat(borderline_stability(1));
if not(isempty(k_critical))    
    p_inter = zeros(1,length(k_critical));
    for i=1:length(k_critical)
        x_temp = real(roots_inter(i));
        y_temp = imag(roots_inter(i));
        if x_temp<0
            sign_x = '-';
        else
            sign_x = '+';
        end
        if y_temp<0
            sign_y = '-';
        else
            sign_y = '+';
        end
        desc_inter = [desc_inter; strcat(sign_x,num2str(round(abs(x_temp),4)),sign_y,num2str(abs(y_temp)),'i, K=',num2str(k_critical(i)))];
        p_inter(i) = plot(x_temp,y_temp,'k*');
    end
    legend_plots = [legend_plots p_inter(1)];
    legend_descriptions = [legend_descriptions [desc_inter]];
end

legend(legend_plots,legend_descriptions);
hold off;

rc = r;
departure_angles = dep_an;
arrival_angles = arr_an;
break_away_in = break_points;
end