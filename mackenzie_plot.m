function [dislist,disaxlist] = mackenzie_plot(n,input_str)

%% Plot mackenzie distribution for n random rotations 
% Set input_str to 'on' to show disorientation histogram (this is where
% mackenzie distribution is shown), as well as axes plot in SST


if nargin == 1
    input_str = 'off';
end 

q_cubic_sym_flip = ... 
[1 0 0 0; ...
 sqrt(2)/2 0 0 sqrt(2)/2; ... 
 sqrt(2)/2 sqrt(2)/2 0 0; ... 
 0 0 0 1; ... 
 0.5 0.5 -0.5 0.5; ... 
 sqrt(2)/2 0 0 -sqrt(2)/2; ... 
 0 0 -sqrt(2)/2 sqrt(2)/2; ... 
 0.5 0.5 0.5 -0.5; ... 
 0 1 0 0;
 0 0 sqrt(2)/2 sqrt(2)/2; ... 
 0 sqrt(2)/2 0 sqrt(2)/2; ... 
 0.5 0.5 -0.5 -0.5; ... 
 0 -sqrt(2)/2 0 sqrt(2)/2; ... 
 sqrt(2)/2 -sqrt(2)/2 0 0; ...
 0 0 1 0; ... 
 0.5 -0.5 -0.5 -0.5; ... 
 0 -sqrt(2)/2 sqrt(2)/2 0; ... 
 0.5 -0.5 0.5 0.5; ... 
 sqrt(2)/2 0 -sqrt(2)/2 0; ... 
 0.5 -0.5 0.5 -0.5; ... 
 sqrt(2)/2 0 sqrt(2)/2 0; ... 
 0.5 -0.5 -0.5 0.5; ... 
 0 sqrt(2)/2 sqrt(2)/2 0; ... 
 0.5 0.5 0.5 0.5];

%putting in format [qx qy qz q0]
q_cubic_sym = zeros(size(q_cubic_sym_flip));
q_cubic_sym(:,1:3) = q_cubic_sym_flip(:,2:4);
q_cubic_sym(:,4) = q_cubic_sym_flip(:,1);

om_cubic_sym_list = qu2omlist(q_cubic_sym); %9x24 matrix
om_cubic_sym = om_reconstruct(om_cubic_sym_list); %cell with 24 3x3 matrices

dislist = zeros(1,n);
disaxlist = zeros(n,3);

for j = 1:n
    axis_un1 = 2*rand(1,3)-1;
    axis_un2 = 2*rand(1,3)-1;
    axis_n1 = axis_un1/sqrt(sum(axis_un1.^2));
    axis_n2 = axis_un2/sqrt(sum(axis_un2.^2));
    angle1 = rand(1)*2*pi;
    angle2 = rand(1)*2*pi;

    % Random Axis Angle Pair
    aa_1 = [axis_n1 angle1];
    aa_2 = [axis_n2 angle2];
    qu_1 = ax2qu(aa_1);
    qu_2 = ax2qu(aa_2);

    % Converted to orientation matrices
%     O1 = vrrotvec2mat(aa_1);
%     O2 = vrrotvec2mat(aa_2);
    

    % Define misorientation g as rotation required to bring ori 1 into ori 2: O2 = misorientation*O1

    %g_unsymmetrized = O2*O1'; %inverse of OM = transpose of OM
    % structure to store symmetric variants of orientation
%     O1_sym_variants = cell(1,24);

    %structure to store symmetric variants of misorientation in different
    %representations
%     g_sym_variants_omcell = cell(1,24);
%     g_sym_variants_omlist = zeros(24,9);
    % g_sym_variants_aa = zeros(24,4);
   % g_sym_variants_qu = zeros(24,4);
    curr_max = 0;
    for i = 1:24
%         O1_sym_var = om_cubic_sym{1,i}*O1;
%         O1_sym_variants{1,i} = O1_sym_var;

        %calculate misorientations
%         g_sym_var_om = O2*O1_sym_var';
%         g_sym_variants_omcell{1,i} = g_sym_var_om;
%         g_sym_variants_omlist(i,:) = reshape(g_sym_var_om',[1 9]);
        
        %with quaternions
        q1_sym_var = qmult(q_cubic_sym(i,:),qu_1);
        qg_sym_var = qmult((q1_sym_var),qinv(qu_2));
        qg_4 = abs(qg_sym_var(4)); %we seek to maximize fourth component
        
        if qg_4 > curr_max
            curr_max = qg_4;
            ax_choose = qg_sym_var(1:3);
        end  
    
 %       g_sym_variants_qu(i,:) = qg_sym_var;
    end

    %postom variants control for conversion math (just using orientation
    %matrices)
%     postom_variants_qulist = om2qulist(g_sym_variants_omlist);
%     postom_variants_alist = qu2alist(postom_variants_qulist);

    disorientation_angle = rad2deg(2*acos(curr_max));
    sortedax = sort(abs(ax_choose));
% 
%     [disorientation_angle, dis_index] = min(rad2deg(2*acos((abs(g_sym_variants_qu(:,4))))));%min(postom_variants_alist(:,4));
    dislist(j) = disorientation_angle;
%     
    %angle placed into the SST convenient for my plotting routine
%     sortedax = sort(abs(g_sym_variants_qu(dis_index,1:3)));
    disaxlist(j,:) = sortedax(:,[2 1 3]);
end

if ~strcmp(input_str,'off')
    figure
    histogram(dislist(dislist ~= 0),100);
    title(['Disorientation Statistics for Randomly Sampled Cubic Orientations, sample size = ',...
        num2str(n)])
    xlabel('Disorientation Angle (°)')
    ylabel('Count')

    [th_dis,r_dis] = stereo(disaxlist);
    figure
    pax = polaraxes;
    fz_edge_p = fzedge_find_m(linspace(0,pi/4));
    polarplot(pax,fz_edge_p(:,1),fz_edge_p(:,2),'k')
    hold on;
    polarscatter(th_dis,r_dis,'filled')
    hold off
    pax.ThetaLim = [0 45];
    title(['Axis Statistics for Randomly Sampled Cubic Orientations, sample size = ',...
        num2str(n)])
end

end

function g = axis_insert(g_un,s)
%Insert axis into FZ

if nargin == 1
    s = 'off';
end 

g_n = abs(g_un)./sqrt(sum(abs(g_un.^2),2)); %normalize rows, take absolute value
g1 = zeros(length(g_n(:,1)),3);

for i = 1:length(g_n(:,1))
    g1(i,:) = (sort((g_n(i,:))));
end

g = g1(:,[2 1 3]); %this is the permutation that works to put axes in FZ!

if strcmp(s,'on')
    [th_g,r_g] = stereo(g);
    figure
    pax = polaraxes;
    fz_edge_p = fzedge_find_m(linspace(0,pi/4));
    polarplot(pax,fz_edge_p(:,1),fz_edge_p(:,2),'k')
    hold on;
    polarscatter(th_g,r_g,'filled')
    hold off
    pax.ThetaLim = [0 45];
    pax.RLim = [0 fzedge_find(pi/4)];
end


end
% take similar elements of an omlist and only considers them once. 
function [omout, szout, th_out, r_out] = omcluster(supertest, input_str)
% supertest = super_cu{185,6};
% qutest = supertest(:,1:4);

if nargin == 1
    input_str = 'off';
end
if length(supertest(1,:)) == 14
    omtest = supertest(:,5:13); 
    sztest = supertest(:,14);
elseif length(supertest(1,:)) == 10
    omtest = supertest(:,1:9);
    sztest = supertest(:,10);
else 
    omtest = supertest;
    sztest = 50000;
end

n = length(omtest(:,1));
omnew = zeros(size(omtest));

%sort each basis vector 
for i = 1:n
    omline = abs(omtest(i,:));
    newomline = [sort(omline(1:3)) sort(omline(4:6)) sort(omline(7:9))];
    omnew(i,:) = newomline;
end 

%for each orientation line, compute difference between other sorted orientation lines
%choose some cutoff below which we call an orientation "similar enough" to
%consider the same. Output unique OM's only to new list.
omold = omnew;
szold = sztest;
indold = 1:length(sztest);

omuniq = zeros(size(omnew));
szuniq = zeros(size(sztest));
induniq = zeros(size(sztest));

k = 0;
while ~isempty(szold)
    
    k = k+1;
    % take peak with most atoms
    [~,ind_max] = max(szold);
    
    %find similar peaks
    cutoff = 0.15;
    omdiff = sum(abs(omold - omold(ind_max,:)),2);
    ind_to_combine = find(omdiff < cutoff);
    ind_to_keep = setdiff(1:length(szold),ind_to_combine);
    %adding sizes of peaks that will be collapsed
    sz_sum = sum(szold(ind_to_combine));

    %output unique orientation and size sum, and index in original model
    omuniq(k,:) = omold(ind_max,:);
    szuniq(k,:) = sz_sum;
    induniq(k,:) = indold(ind_max);
    

    %delete current equivalent orientations from structure
    omold(ind_to_combine,:) = [];
    szold(ind_to_combine) = [];
    indold(ind_to_combine) = [];
end 

induniq(szuniq == 0,:) = [];
omuniq(szuniq == 0,:) = [];
szuniq(szuniq == 0) = [];

omuniq(:,end+1) = szuniq;
omuniq(:,end+1) = induniq;

omuniq = abs(sortrows(-omuniq,10));

% omout = omuniq(:,1:9);
szout = omuniq(:,10);
indout = omuniq(:,11);
omout = omtest(indout,:);

[th_out, r_out] = triple_stereo(omout);

if ~strcmp(input_str,'off')
ratio_sz = szout/sum(szout);
    for i = 1:length(omout(:,1))
        disp(['Peak ',num2str(i),' : ', num2str(round(ratio_sz(i)*100,2)), '% of atoms in layer'])
        disp(reshape(omout(i,:),[3 3])')
    end
end

end

function omlist = qu2omlist(qq)
%quaternion matrix to orientation matrix cell



qbar = qq(:,4).*qq(:,4)-(qq(:,1).*qq(:,1)+qq(:,2).*qq(:,2)+qq(:,3).*qq(:,3));

q = cell(3);

q{1,1} = qbar + 2.0*qq(:,1).*qq(:,1);
q{2,2} = qbar + 2.0*qq(:,2).*qq(:,2);
q{3,3} = qbar + 2.0*qq(:,3).*qq(:,3);

q{1,2} = 2.0*(qq(:,1).*qq(:,2)-qq(:,4).*qq(:,3));
q{2,3} = 2.0*(qq(:,2).*qq(:,3)-qq(:,4).*qq(:,1));
q{3,1} = 2.0*(qq(:,3).*qq(:,1)-qq(:,4).*qq(:,2));
q{2,1} = 2.0*(qq(:,2).*qq(:,1)+qq(:,4).*qq(:,3));
q{3,2} = 2.0*(qq(:,3).*qq(:,2)+qq(:,4).*qq(:,1));
q{1,3} = 2.0*(qq(:,1).*qq(:,3)+qq(:,4).*qq(:,2));

omlist = [q{1,:} q{2,:} q{3,:}];

thr = 1e-7;
for i = 1:length(omlist(:,1))
    omi = omlist(i,:);
    omi(abs(omi) < thr) = 0.0;
    omlist(i,:) = omi;
end 


end

function omcell = om_reconstruct(omlist)
%takes n x 9 list, reconstructs rows as 9 linear indices of matrix 

dim = size(omlist);
if dim(1) == 1 && dim(2) == 9
    omcell = zeros(3);
    omcell(1,:) = omlist(1:3);
    omcell(2,:) = omlist(4:6);
    omcell(3,:) = omlist(7:9);
else
    n = length(omlist(:,1));
    omcell = cell(1,n);

    for i = 1:n
        omtemp = zeros(3);
        omtemp(1,:) = omlist(i,1:3);
        omtemp(2,:) = omlist(i,4:6);
        omtemp(3,:) = omlist(i,7:9);
        omcell{1,i} = omtemp;
    end 

end

end 

function [th,r] = triple_stereo(test_omp1)

if length(test_omp1(1,:)) >= 9
    o_x1 = test_omp1(:,1:3);
    o_y1 = test_omp1(:,4:6);
    o_z1 = test_omp1(:,7:9);

    th = zeros(size(o_x1));
    r = zeros(size(o_x1));

    [thx,rx] = stereo(axis_insert(o_x1));
    [thy,ry] = stereo(axis_insert(o_y1));
    [thz,rz] = stereo(axis_insert(o_z1));

    th(:,1) = thx;
    th(:,2) = thy;
    th(:,3) = thz;

    r(:,1) = rx;
    r(:,2) = ry;
    r(:,3) = rz;
    
elseif length(test_omp1(1,:)) < 9 %quaternion
    o_x1 = test_omp1(:,1:3);
    [th,r] = stereo(axis_insert(o_x1));

end

end

function [Theta,R] = stereo(aa)
%Stereographic projection of rotation axes
% n x 3 matrix, columns are unit vectors [ax ay az]
n = length(aa(:,1));
% r = zeros(1,n);
R = zeros(1,n);
Theta = zeros(1,n);
%Phi = zeros(1,n);

for i = 1:n
    a = aa(i,1:3);
    r = norm(a);
    if a(1) < 0
        Theta(i) = atan(a(2)/a(1)) + pi;
    else 
        Theta(i) = atan(a(2)/a(1));
    end 
    phi = acos(a(3)/r);
    R(i) = sin(pi-phi)/(1-cos(pi-phi));
end 

R(isnan(R)) = 0;
Theta(isnan(Theta)) = 0;
% 
% figure
% pax = polaraxes;
% polarscatter(Theta,R)
% title('stereographic projection of rotation axes')
% pax.ThetaLim = [0 90];
% pax.RLim = [0 2];
% 
% figure
% polarhistogram(Theta)
% 
% figure
% histogram(R)
% figure
% histogram2(Theta,R)
% xlabel('Theta')
% ylabel('R')

end

function fz_edge = fzedge_find_m(th_list)
%Finds edge of standard triangle for specified quadruplet of r / theta
%pairs. 
%

fz_edge = zeros(length(th_list),2);
fz_shape = stereomat([0 0 1; 1 1 1; 1 0 1; -1 0 -1]);

r_new = (fz_shape(4,2) + fz_shape(3,2))/2.;
r_old = fz_shape(2,2);

theta1_span = fz_shape(2,1);
theta2_span = pi-acos((r_old^2 - r_new^2 - (r_new-fz_shape(3,2))^2)/(2*r_new*(r_new-fz_shape(3,2)))); %law of cosines, draw picture
%value above is insensitive to large angles 

th_ratio = theta1_span/theta2_span; 
% now we express the new curve in terms of polar coordinates of the first
fz_edge_pre = th_list/th_ratio;
fz_edge_x = r_new*cos(fz_edge_pre)-(r_new-fz_shape(3,2));%(r_new-r_old*cos(theta1_span));
fz_edge_y = r_new*sin(fz_edge_pre);
fz_edge_r = sqrt(fz_edge_x.^2 + fz_edge_y.^2);
fz_edge_th = atan(fz_edge_y ./ fz_edge_x);

fz_edge(:,1) = fz_edge_th;
fz_edge(:,2) = fz_edge_r;

end

function thr = stereomat(aa)
%Stereographic projection of rotation axes, theta column and r column
% n x 3 matrix, columns are unit vectors [ax ay az]
% assumes only positive axes values
n = length(aa(:,1));
thr = zeros(n,2);

ax = aa(:,1);
ay = aa(:,2);
az = aa(:,3);
r = sqrt(ax.^2 + ay.^2 + az.^2);
thr(:,1) = atan(ay./ax); %theta
phi = acos(az./r);
thr(:,2) = sin(pi*ones(length(phi),1)-phi)./(ones(length(phi),1)-cos(pi*ones(length(phi),1)-phi)); %R

% set [100] projection to origin
thr(isnan(thr)) = 0;

% 
% Theta = thr(:,1);
% R = thr(:,2);
% 
% figure
% pax = polaraxes;
% polarscatter(Theta,R)
% title('stereographic projection of rotation axes')
% pax.ThetaLim = [0 90];
% pax.RLim = [0 2];
% 
% figure
% polarhistogram(Theta)
% 
% figure
% histogram(R)
% figure
% histogram2(Theta,R)
% xlabel('Theta')
% ylabel('R')

end
function fz_edge_r = fzedge_find(th_list)
%Finds edge of standard triangle for specified quadruplet of r / theta
%pairs. 
%


fz_shape = stereomat([0 0 1; 1 1 1; 1 0 1; -1 0 -1]);

r_new = (fz_shape(4,2) + fz_shape(3,2))/2.;
r_old = fz_shape(2,2);

theta1_span = fz_shape(2,1);
theta2_span = pi-acos((r_old^2 - r_new^2 - (r_new-fz_shape(3,2))^2)/(2*r_new*(r_new-fz_shape(3,2)))); %law of cosines, draw picture
%value above is insensitive to large angles 

th_ratio = theta1_span/theta2_span; 
% now we express the new curve in terms of polar coordinates of the first
fz_edge_pre = th_list/th_ratio;
fz_edge_x = r_new*cos(fz_edge_pre)-(r_new-fz_shape(3,2));%(r_new-r_old*cos(theta1_span));
fz_edge_y = r_new*sin(fz_edge_pre);
fz_edge_r = sqrt(fz_edge_x.^2 + fz_edge_y.^2);
%fz_edge_th = atan(fz_edge_y ./ fz_edge_x);


end

% from axis-angle pair to quaternions

function q = ax2qu(ax)

thr = 1e-10;
if (abs(ax(4)-0.0)<thr)
   q = [ 1.0, 0.0, 0.0, 0.0 ];
else
   c = cos(ax(4)*0.5);
   s = sin(ax(4)*0.5);
   q = [ c, ax(1)*s, ax(2)*s, ax(3)*s ];
end 

% set values very close to 0 as 0
if (abs(q(1))-0)<thr
    q(1)=0;
elseif (abs(q(2))-0)<thr
    q(2)=0;
elseif (abs(q(3))-0)<thr
    q(3)=0;
elseif (abs(q(4))-0)<thr
    q(4)=0;
end
end

function q = qmult(q1,q2)
%quaternion multiplication
%quaternion input and output in form (qx, qy, qz, w)

q1v = q1(:,1:3);
q1w = q1(:,4);
q2v = q2(:,1:3);
q2w = q2(:,4);

qw = q1w.*q2w - dot(q1v,q2v,2);
qv = bsxfun(@times,q1w,q2v) + bsxfun(@times,q2w,q1v) + cross(q1v,q2v,2);

q = zeros(size(q1));
q(:,1:3) = qv;
q(:,4) = (qw); %this was changed May 2017

end

function q = qinv(q1)
%inverse of unit quaternion
%input and output: scalar element last 

q = zeros(size(q1));
q(:,1:3) = -q1(:,1:3);
q(:,4) = q1(:,4);

end




