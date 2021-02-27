function  [rc,departure_angles,arrival_angles,break_away_in,asymptotes,borderline_stability] = root_locus_calculation(sys)
% Root locus calculation
[rc,k] = rlocus(sys);
% Departure and arrival angles calculation
len = length(k);
num_tragectories = length(rc(:,1));
P = pole(sys);
Z = zero(sys);
n = length(P);
m = length(Z);
departure_angles = [];
arrival_angles = [];
for i=1:num_tragectories
   y1 = imag(rc(i,1));
   y2 = imag(rc(i,2));
   x1 = real(rc(i,1));
   x2 = real(rc(i,2));
   theta = atan2(y2-y1,x2-x1)*180/pi;
   departure_angles = [departure_angles theta];
   if 1/imag(rc(i,len))~=0 && 1/real(rc(i,len))~=0
       y1 = imag(rc(i,len-1));
       y2 = imag(rc(i,len));
       x1 = real(rc(i,len-1));
       x2 = real(rc(i,len));    
       theta = 180 + atan2(y2-y1,x2-x1)*180/pi;
       arrival_angles = [arrival_angles theta];
   end   
end

% Breakaway/in points calculation
syms s;
[Num,Den] = tfdata(sys);
sys_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
sys_derivative = diff(sys_syms);
possible_points = double(solve(sys_derivative==0));
break_away_in = [];
for i=1:length(possible_points)
    temp = double(subs(sys_syms,s,possible_points(i)));
    if temp<0&&imag(temp)==0
        break_away_in = [break_away_in possible_points(i)];
    end
end

% Asymptotes
asymp_coord = - (sum(-P)-sum(-Z))/(n-m);
asymp_angles = [];
for p=0:n-m-1
temp_theta = (2*p+1)*180/(n-m);
asymp_angles = [asymp_angles temp_theta];  
end
asymp_lamda = tand(asymp_angles);
asymptotes = {asymp_coord asymp_angles asymp_lamda};

% Intersection with y'y - Borderline stability
inter_pts = [];
inter_K = [];
for i=1:num_tragectories
    curr_trajectory = real(rc(i,2:end));
    condition = 0.1;
    while(length(find(abs(curr_trajectory)<condition))>1)
        condition = condition*0.1;
    end
    temp_index = find(abs(curr_trajectory)<condition)+1;
    inter_pts = [inter_pts rc(i,temp_index)];
    inter_K = [inter_K k(temp_index)];
end
borderline_stability = {inter_pts,inter_K};
end
