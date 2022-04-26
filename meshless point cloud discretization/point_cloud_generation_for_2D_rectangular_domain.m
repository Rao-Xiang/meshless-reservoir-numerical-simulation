function [ point_cloud_information ] = point_cloud_generation_for_2D_rectangular_domain(  )
% this function is used for point cloud generation of 2D_rectangular_domain
% [fracture_lines, fracture_length] = fracture_lines_solver_1();
%% inputs
load('fracture_lines.mat')
x_size = 800;
y_size = 610;
char_length = 10; % characteristic length for matrix nodes
char_fracture = 10; % characteristic length for fracture nodes
matching_flag = 1; % 0 denotes non-matching point cloud, 1 denotes matching point cloud
input_judge_overlapping_distance_matr_frac = 5;
%% default parameters
char_deltx = char_length;
char_delty = char_length;
char_welloc = char_length;
if matching_flag == 0
    judge_overlapping_distance_matr_frac = -100;% it is the default
else
    judge_overlapping_distance_matr_frac = input_judge_overlapping_distance_matr_frac; % if matching, need to input
end
test_fract_intersections = 1e-4;% It is the default
nf = size(fracture_lines,1)/2;
fracture_length = zeros(nf,1);
for i = 1:nf
    fracture_length(i,1) = norm(fracture_lines(2*i,:)-fracture_lines(2*i-1,:));
end
%% figure 1: plot the reservoir domain
figure('color','w')
for i = 1:nf
    fracture_2_node_1 = fracture_lines(2*i-1,:);
    fracture_2_node_2 = fracture_lines(2*i,:);
    plot([fracture_2_node_1(1,1);fracture_2_node_2(1,1)],[fracture_2_node_1(1,2);fracture_2_node_2(1,2)],'k');
    hold on;
end
x_framework = [0;x_size;x_size;0];
y_framework = [0;0;y_size;y_size];
line(x_framework,y_framework,'color','k');
hold on;
% % % vertical_wells = [
% % %         55, 455;
% % % %         735, 90;
% % % %         95, 515;
% % %         755, 55;
% % %         ];
% % % scatter(vertical_wells(:,1),vertical_wells(:,2),20,'b');%,'filled'
axis equal
xlabel('x, m');
ylabel('y, m');
xlim([0,x_size]);
ylim([0,y_size]);
title('Computational reservoir domain');
hold off;
%% STEP #1: compute fracture intersections
interfractures = [];
inter_ksais = [];
numbering_fracture_intersections = [];
coordinates_fracture_intersections = [];
flag_fracture_intersections = 0;
interfrac_cross = zeros(nf);
for i = 1:nf-1
    fracture_1_node_1 = fracture_lines(2*i-1,:);
    fracture_1_node_2 = fracture_lines(2*i,:);
    vector_fracture_1 = fracture_1_node_2 - fracture_1_node_1;
    for j = i+1:nf
        fracture_2_node_1 = fracture_lines(2*j-1,:);
        fracture_2_node_2 = fracture_lines(2*j,:);
        vector_fracture_2 = fracture_2_node_2 - fracture_2_node_1;
        coef_matrix = [vector_fracture_1(1),-vector_fracture_2(1);
            vector_fracture_1(2),-vector_fracture_2(2);];
        coef_vector = [fracture_2_node_1(1)-fracture_1_node_1(1);
            fracture_2_node_1(2)-fracture_1_node_1(2);];
        if abs(det(coef_matrix))<1e-12 % the two fractures are parallel
        else
            ksai = coef_matrix\coef_vector;
            if ksai(1)>=0 && ksai(1)<=1 && ksai(2)>=0 && ksai(2)<=1
                the_coordinates_1 = fracture_1_node_1+ ksai(1)*(fracture_1_node_2-fracture_1_node_1);
                the_coordinates_2 = fracture_2_node_1+ ksai(2)*(fracture_2_node_2-fracture_2_node_1);
                % Check that the point is indeed an intersection point to improve the robustness of the program
                if norm(the_coordinates_2-the_coordinates_1)>=test_fract_intersections
                    error('This may be an error because the difference between the coordinates of the two points is relatively large and cannot be considered as two coincident points.');
                end
                interfrac_cross(i,j) = 1; interfrac_cross(j,i) = 1;
                flag_fracture_intersections = flag_fracture_intersections +1;
                interfractures = [interfractures; i, j];
                inter_ksais = [inter_ksais; ksai(1), ksai(2)];
                flag_fracture_intersections = flag_fracture_intersections +1;
                numbering_fracture_intersections = [numbering_fracture_intersections; flag_fracture_intersections];
                coordinates_fracture_intersections = [coordinates_fracture_intersections; the_coordinates_1];
            end
        end
        
    end
end
%% STEP #2: generation of fracture nodes
fractures_ksais = cell(nf,1);
%
for i = 1:nf
    fractures_ksais{i,1}=[0;1];
end
%
number_of_frac_intersect = size(interfractures,1);
for i = 1:number_of_frac_intersect
    fractures_ksais{interfractures(i,1),1} = [fractures_ksais{interfractures(i,1),1};inter_ksais(i,1)];
    fractures_ksais{interfractures(i,2),1} = [fractures_ksais{interfractures(i,2),1};inter_ksais(i,2)];
end
for i = 1:nf
    sub_fracture_ksais = fractures_ksais{i,1};
    sorted_sub_fracture_ksais = unique(sub_fracture_ksais);
    fractures_ksais{i,1} = sorted_sub_fracture_ksais;
end
% fracture nodes denoted by ksai
for i = 1:nf
    sub_fracture_ksais = fractures_ksais{i,1};
    number_sub_fracture_ksais = size(sub_fracture_ksais,1);
    temp_sub_fracture_ksais = [];
    for j = 1:(number_sub_fracture_ksais-1)
        ksai_1 = sub_fracture_ksais(j);
        ksai_2 = sub_fracture_ksais(j+1);
        if round(fracture_length(i)*(ksai_2-ksai_1)/char_fracture)==0
            number_fracture_segments = 1;
        else
            number_fracture_segments = round(fracture_length(i)*(ksai_2-ksai_1)/char_fracture);
        end
        real_char_fracture_legnth = fracture_length(i)*(ksai_2-ksai_1)/number_fracture_segments;
        sub_sub_fracture_ksais = (linspace(ksai_1,ksai_2,number_fracture_segments+1))';
        temp_sub_fracture_ksais = [temp_sub_fracture_ksais; sub_sub_fracture_ksais];
    end
    real_sub_fracture_ksais = unique(temp_sub_fracture_ksais);
    fractures_ksais{i,1} = real_sub_fracture_ksais;
end
% fracture-node coordinates
fracture_nodes = cell(nf,1);
fracture_nodes_numbering = [];
fracture_nodes_coordinates = [];
node_flag = 0;
for i = 1:nf
    fracture_i_node_1 = fracture_lines(2*i-1,:);
    fracture_i_node_2 = fracture_lines(2*i,:);
    sorted_sub_fracture_ksais = fractures_ksais{i,1};
    sub_fracture_nodes = [];
    for j = 1:size(sorted_sub_fracture_ksais,1)
        node_flag = node_flag +1;
        sub_fracture_nodes = [sub_fracture_nodes; node_flag];
        fracture_nodes_numbering = [fracture_nodes_numbering; node_flag];
        fracture_nodes_coordinates = [fracture_nodes_coordinates; fracture_i_node_1+ sorted_sub_fracture_ksais(j)*(fracture_i_node_2-fracture_i_node_1)];
    end
    fracture_nodes{i,1} = sub_fracture_nodes;
end
fracture_nodes_information = [fracture_nodes_numbering, fracture_nodes_coordinates];
% find the points on the two fractures corresponding to the intersection of
% the two fractures
overlapping_pairs_interfractures = [];
need_changed_flag = zeros(size(fracture_nodes_numbering,1),1);
for i= 1:size(interfractures,1)
    the_one_fracture = interfractures(i,1);
    the_other_one_fracture = interfractures(i,2);
    the_coordinates_intersection = coordinates_fracture_intersections(i,:);
    the_one_sub_fracture_nodes = fracture_nodes{the_one_fracture,1};
    the_other_one_sub_fracture_nodes = fracture_nodes{the_other_one_fracture,1};
    
    [aaa,bbb] = min((fracture_nodes_coordinates(the_one_sub_fracture_nodes,1)-the_coordinates_intersection(1)).^2+(fracture_nodes_coordinates(the_one_sub_fracture_nodes,2)-the_coordinates_intersection(2)).^2);
    [ccc,ddd] = min((fracture_nodes_coordinates(the_other_one_sub_fracture_nodes,1)-the_coordinates_intersection(1)).^2+(fracture_nodes_coordinates(the_other_one_sub_fracture_nodes,2)-the_coordinates_intersection(2)).^2);
    if aaa(1)>test_fract_intersections || ccc(1)>test_fract_intersections
        error('There is a problem with the code because there must be a fracture point that overlaps with this inter-fracture intersection')
    end
    overlapping_pairs_interfractures = [overlapping_pairs_interfractures; the_one_sub_fracture_nodes(bbb(1),1),the_other_one_sub_fracture_nodes(ddd(1),1)];
    need_changed_flag(the_other_one_sub_fracture_nodes(ddd(1),1),1) = 1;
end

% Delete duplicate fracture nodes
sorted_fracture_nodes_numbering = [];
for i = 1:size(fracture_nodes_numbering,1)
    if ismember(i,overlapping_pairs_interfractures(:,2))
    else
        sorted_fracture_nodes_numbering = [sorted_fracture_nodes_numbering; i];
    end
end
%% STEP #3: generation of matrix nodes
number_x_matrix_nodes = round(x_size/char_deltx);
number_y_matrix_nodes = round(y_size/char_delty);
real_deltx = x_size/number_x_matrix_nodes;
real_delty = y_size/number_y_matrix_nodes;
number_x_matrix_nodes = number_x_matrix_nodes +1;
number_y_matrix_nodes = number_y_matrix_nodes +1;
matrix_nodes_numbering = size(fracture_nodes_numbering,1)+(1:(number_x_matrix_nodes*number_y_matrix_nodes))';
matrix_nodes_coordinates = zeros(number_x_matrix_nodes*number_y_matrix_nodes,2);
for i = 1:(number_x_matrix_nodes)
    for j = 1:(number_y_matrix_nodes)
        matrix_nodes_coordinates(i+(j-1)*number_x_matrix_nodes,:) = [(i-1)*real_deltx, (j-1)*real_delty];
    end
end
matrix_nodes_information = [matrix_nodes_numbering, matrix_nodes_coordinates];

nodes_information = [fracture_nodes_information;matrix_nodes_information];

% in the matching point cloud, to determine which matrix nodes can be replaced by a crack node because
% they are too close to the crack node, since the crack node will accordingly be present in the matrix system as a matrix node
% in the non-matching point cloud, although it is no need to do so, it is not necessary to set
% the code separately by adding a judgement, but simply by making judge_overlapping_distance_matr_frac < 0 in the code below,
% i.e. making judge_overlapping_distance_matr_frac the default of -100 as set.
overlapping_matrix_fracture_also_matrix_pairs = [];
need_changed_flag_matrix_nodes = zeros(size(matrix_nodes_numbering,1),1);
for i = 1:size(matrix_nodes_numbering,1)
    if need_changed_flag_matrix_nodes(i,1) ~=1
        the_matrix_node = matrix_nodes_numbering(i);
        for j = 1:size(sorted_fracture_nodes_numbering,1)
            the_fracture_also_matrix_node = sorted_fracture_nodes_numbering(j);
            if norm(matrix_nodes_coordinates(i,:) - fracture_nodes_coordinates(the_fracture_also_matrix_node,:)) <= judge_overlapping_distance_matr_frac
                need_changed_flag_matrix_nodes(i) = 1;
                overlapping_matrix_fracture_also_matrix_pairs = [overlapping_matrix_fracture_also_matrix_pairs; the_fracture_also_matrix_node, the_matrix_node];
                break;
            end
        end
    end
end

sorted_matrix_nodes_numbering = [];
if size(overlapping_matrix_fracture_also_matrix_pairs,1)>=1
    for i = 1:size(matrix_nodes_numbering,1)
        if ismember(matrix_nodes_numbering(i),overlapping_matrix_fracture_also_matrix_pairs(:,2))
        else
            sorted_matrix_nodes_numbering = [sorted_matrix_nodes_numbering; matrix_nodes_numbering(i)];
        end
    end
elseif judge_overlapping_distance_matr_frac <=0
    for i = 1:size(matrix_nodes_numbering,1)
        sorted_matrix_nodes_numbering = [sorted_matrix_nodes_numbering; matrix_nodes_numbering(i)];
    end
end

sorted_nodes_numbering = [sorted_fracture_nodes_numbering; sorted_matrix_nodes_numbering];

real_nodes_nubering = (1:size(sorted_nodes_numbering,1))';
real_nodes_coordinates = zeros(size(sorted_nodes_numbering,1),2);
for i = 1:size(real_nodes_nubering,1)
    real_nodes_coordinates(i,:) = nodes_information(sorted_nodes_numbering(i),2:3);
end
%% STEP 4
%% figure 2
for i = 1:size(sorted_fracture_nodes_numbering,1)
    pre_fracture_nodes_coordinates(i,:) = nodes_information(sorted_fracture_nodes_numbering(i),2:3);
end
for i = 1:size(sorted_matrix_nodes_numbering,1)
    pre_matrix_nodes_coordinates(i,:) = nodes_information(sorted_matrix_nodes_numbering(i),2:3);
end

if matching_flag ~=0 
    Type_I_matrix_nodes = pre_fracture_nodes_coordinates;
    Type_II_matrix_nodes = pre_matrix_nodes_coordinates;
    matrix_nodes_coordinates = real_nodes_coordinates;
    fracture_nodes_coordinates = pre_fracture_nodes_coordinates;
    point_cloud_information.Type_I_matrix_nodes = Type_I_matrix_nodes;
    point_cloud_information.Type_II_matrix_nodes = Type_II_matrix_nodes;
    point_cloud_information.fracture_nodes_coordinates = fracture_nodes_coordinates;
    point_cloud_information.matrix_nodes_coordinates = matrix_nodes_coordinates;
else 
    matrix_nodes_coordinates = pre_matrix_nodes_coordinates;
    fracture_nodes_coordinates = pre_fracture_nodes_coordinates;
    point_cloud_information.fracture_nodes_coordinates = fracture_nodes_coordinates;
    point_cloud_information.matrix_nodes_coordinates = matrix_nodes_coordinates;
end
relation_between_real_and_sorted_nodes_numbering = [real_nodes_nubering,sorted_nodes_numbering];

% revise numbering information of fracture nodes because there are duplicate fracture nodes deleted 
original_nodes_numbering = [fracture_nodes_numbering; matrix_nodes_numbering];
relation_between_sorted_and_original_nodes_numbering = [original_nodes_numbering, original_nodes_numbering];
total_overlapping_pairs = [overlapping_pairs_interfractures; overlapping_matrix_fracture_also_matrix_pairs];
for i = 1:size(original_nodes_numbering,1)
    if ismember(original_nodes_numbering(i),total_overlapping_pairs(:,2))
        a = find(total_overlapping_pairs(:,2)==original_nodes_numbering(i));
        if size(a,1)>1
            error('代码有错误，某节点多次消失')
        else
            relation_between_sorted_and_original_nodes_numbering(i,1) = total_overlapping_pairs(a,1);
        end
    end
end

for i = 1:nf
    sub_fracture_nodes = fracture_nodes{i,1};
    for j = 1:size(sub_fracture_nodes,1)
        a = find(relation_between_sorted_and_original_nodes_numbering(:,2)==sub_fracture_nodes(j));
        if size(a,1)>1
            error('There is an error in the code. A node disappears many times')
        else
            the_sorted_node = relation_between_sorted_and_original_nodes_numbering(a,1);
            b = find(relation_between_real_and_sorted_nodes_numbering(:,2)==the_sorted_node);
            if size(b,1)>1
                error('There is an error in the code. Please check it')
            else
                the_real_node = relation_between_real_and_sorted_nodes_numbering(b,1);
                sub_fracture_nodes(j) = the_real_node;
            end
        end
    end
    fracture_nodes{i,1} = sub_fracture_nodes;
end
point_cloud_information.fracture_nodes = fracture_nodes;
%%  figure 2: plot all nodes
figure('color','w')
scatter(matrix_nodes_coordinates(:,1),matrix_nodes_coordinates(:,2),5,'filled');
hold on;
scatter(fracture_nodes_coordinates(:,1),fracture_nodes_coordinates(:,2),5,'filled','k'); % 裂缝节点
hold on;
x_framework = [0;x_size;x_size;0];
y_framework = [0;0;y_size;y_size];
line(x_framework,y_framework,'color','k');
hold off;
xlabel('x, m');
ylabel('y, m');
axis equal;
xlim([0,x_size]);
ylim([0,y_size]);
title('All nodes: matrix nodes and fracture nodes')
hold off
%% figure 3: plot matrix nodes
figure('color','w')
if matching_flag ~=0
    scatter(Type_II_matrix_nodes(:,1),Type_II_matrix_nodes(:,2),5,'filled');
    hold on;
    scatter(Type_I_matrix_nodes(:,1),Type_I_matrix_nodes(:,2),5,'filled','r');
    hold on
    x_framework = [0;x_size;x_size;0];
    y_framework = [0;0;y_size;y_size];
    line(x_framework,y_framework,'color','k');
    hold off;
    xlabel('x, m');
    ylabel('y, m');
    axis equal;
    xlim([0,x_size]);
    ylim([0,y_size]);
    title('Matrix nodes')
else % 代表非匹配性点云
    scatter(pre_matrix_nodes_coordinates(:,1),pre_matrix_nodes_coordinates(:,2),5,'filled');
    hold on;
    x_framework = [0;x_size;x_size;0];
    y_framework = [0;0;y_size;y_size];
    line(x_framework,y_framework,'color','k');
    hold off;
    xlabel('x, m');
    ylabel('y, m');
    axis equal;
    xlim([0,x_size]);
    ylim([0,y_size]);
    title('Matrix nodes')
end
%% figure 4
figure('color','w')
% scatter(real_nodes_coordinates(:,1),real_nodes_coordinates(:,2),8,'filled');
scatter(pre_fracture_nodes_coordinates(:,1),pre_fracture_nodes_coordinates(:,2),5,'filled','k');
hold on
for i = 1:nf
    sub_fracture_nodes = fracture_nodes{i,1};
    %         num_sub_fracture_nodes = size(sub_fracture_nodes,1);
    plot(real_nodes_coordinates(sub_fracture_nodes,1),real_nodes_coordinates(sub_fracture_nodes,2),'k');
    hold on;
end
x_framework = [0;x_size;x_size;0];
y_framework = [0;0;y_size;y_size];
line(x_framework,y_framework,'color','k');
hold off;
xlabel('x, m');
ylabel('y, m');
axis equal
xlim([0,x_size]);
ylim([0,y_size]);
title('Fracture nodes and fractures')
%% non-matching well-node treatment
% fractured horizontal well
fractured_horizontal_wells = [10,310;790,310;];
number_fractured_horizontal_wells = size(fractured_horizontal_wells,1)/2;
welloc = cell(number_fractured_horizontal_wells,1);
welloc_coordinates = cell(number_fractured_horizontal_wells,1);
figure('color','w')
for i = 1:number_fractured_horizontal_wells
    node_1_coord = fractured_horizontal_wells(2*i-1,:);
    node_2_coord = fractured_horizontal_wells(2*i,:);
    welloc_i = [];
    welloc_coordinates_i = [];
    for j = 1:nf
        sub_fracture_nodes = fracture_nodes{j,1};
        dis_fracture_nodes_to_well = zeros(size(sub_fracture_nodes,1),1);
        for k = 1:size(sub_fracture_nodes,1)
            the_fracture_node = sub_fracture_nodes(k);
            the_fracture_node_coord = real_nodes_coordinates(the_fracture_node,:);
            vector_1 = the_fracture_node_coord-node_1_coord;
            vector_2 = node_2_coord-node_1_coord;
            pre = dot(vector_1,vector_2)/norm(vector_2);
            dis_fracture_nodes_to_well(k,1) = sqrt(norm(vector_1).^2-pre.^2);
        end
        [aa, bb] = min(dis_fracture_nodes_to_well);
        if aa<=char_welloc
            if ~ismember(sub_fracture_nodes(bb(1),1),welloc_i)
                welloc_i = [welloc_i; sub_fracture_nodes(bb(1),1)];
                welloc_coordinates_i = [welloc_coordinates_i; real_nodes_coordinates(sub_fracture_nodes(bb(1),1),:)];
            end
        end
    end
    if size(welloc_i,1)==0
        error('No perforation points were obtained for this fractured horizontal well');
    else
        welloc{i,1} = welloc_i;
        welloc_coordinates{i,1} = welloc_coordinates_i;
        scatter(real_nodes_coordinates(welloc_i,1),real_nodes_coordinates(welloc_i,2));hold on;
    end
end
x_framework = [0;x_size;x_size;0];
y_framework = [0;0;y_size;y_size];
line(x_framework,y_framework,'color','k');
hold off;
xlabel('x, m');
ylabel('y, m');
axis equal
xlim([0,x_size]);
ylim([0,y_size]);
title('Perforating points of fractured horizontal wells');
% vertical well
vertical_wells = [
    55, 455;
    735, 90;
    95, 515;
    755, 145;        ];
number_vertical_wells = size(vertical_wells,1);
node_numbering_vertical_wells = zeros(number_vertical_wells,1);
for i = 1:number_vertical_wells
    the_vertical_well_loc = vertical_wells(i,:);
    distance = (matrix_nodes_coordinates(:,1)-the_vertical_well_loc(1)).^2 + (matrix_nodes_coordinates(:,2)-the_vertical_well_loc(2)).^2;
    [aa,bb] = min(distance);
    node_numbering_vertical_wells(i) = bb(1);
end
point_cloud_information.welloc = welloc;
point_cloud_information.welloc_coordinates = welloc_coordinates;
point_cloud_information.fracture_nodes = fracture_nodes;
point_cloud_information.node_numbering_vertical_wells = node_numbering_vertical_wells;
end
%

%