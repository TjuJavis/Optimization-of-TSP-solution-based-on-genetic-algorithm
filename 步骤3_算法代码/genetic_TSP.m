function genetic_TSP
    clc; close all; clear;
    warning off;
    % 设定参数

    % 选取数据集编号(即城市数量)。可供选择的数据集有三：101、30、16.其中30是自编数据，其余来自于TSPLIB数据集
    DataSet = 101; 
    % 获取城市数量位置、数据集最优解，并设置样本数及迭代次数
    [opt, times, s, N, city_coordinate] = tsp(DataSet); 

    c = 0.3 * s; %设置精英个数
    pc = 0.9; % 交叉概率
    pm = 0.2; % 变异概率

    % 也可自行重设样本数和迭代次数及精英个数
%     s = 100;
%     times = 5000; % 最大迭代次数
%     c = 30; % 精英个数

    time = 0; % 实际迭代次数 
    pop = zeros(s, N + 1); % 初始种群+适应度
    pop_fit_aver = []; % 总适应度
    min_dis = []; % 最短距离
    pop_min = []; % 最短距离的个体
    
    % 初始化种群
    for i = 1:s
        pop(i, 1:N) = randperm(N);
    end
    
    % 画出城市位置分布图
    clf
    plot(city_coordinate(:, 1), city_coordinate(:, 2), 'ro'); % 画城市位置散点图
    for i = 1:N
        test_t = num2str(i);
        text(city_coordinate(i, 1), city_coordinate(i, 2), test_t); % 标号
    end
    grid on;
    title('城市位置分布图');
    xlabel('x');
    ylabel('y');
    
    % 计算城市间距离
    city_distance = CityDistance(city_coordinate, N);
    
    % 计算初始种群的适应度
    [individual_fit, sum, min1, min_index] = GroupFit(city_distance, N, pop, s);
    pop_fit_aver = [pop_fit_aver; sum];
    min_dis = [min_dis; min1];
    pop(:, N + 1) = individual_fit;
    pop_min = [pop_min; pop(min_index, :)];
    pop = ChooseParents(pop, N, s, c); % 选择父代
    
    % 开始遗传算法迭代并计时
    disp('正在使用遗传算法迭代求解TSP中……');
    tic; % 开始计时
    for i = 1:times
        time = time + 1;
        E_new_new = zeros(s, N + 1); % 子代
        for j = 1:s/2
            a = rand(1);
            b = rand(1);
            if a <= pc % 交叉
                [E_new_new(j, 1:N), E_new_new(j + s/2, 1:N)] = CrossVariation(pop(j, 1:N), pop(j + s/2, 1:N), N);
            else
                E_new_new(j, 1:N) = pop(j, 1:N);
                E_new_new(j + s/2, 1:N) = pop(j + s/2, 1:N);
            end
            if b <= pm % 变异
                E_new_new(j, 1:N) = Mutation(E_new_new(j, 1:N), N);
                E_new_new(j + s/2, 1:N) = Mutation(E_new_new(j + s/2, 1:N), N);
            end
        end
        [individual_fit, sum, min1, min_index] = GroupFit(city_distance, N, E_new_new, s);
        pop_fit_aver = [pop_fit_aver; sum];
        min_dis = [min_dis; min1];
        E_new_new(:, N + 1) = individual_fit;
        pop_min = [pop_min; E_new_new(min_index, :)];
        pop = ChooseParents(E_new_new, N, s, c);
    end

    % 结束计时并输出经过的时间
    elapsed_time = toc;
    disp(['遗传求解TSP迭代消耗时间：', num2str(elapsed_time), ' 秒']);

    % 输出结果
    [a, min_index] = min(min_dis);
    disp(['最短路径（km） = ', num2str(a)]);
    dr = (a - opt) / opt * 100; % 计算偏差量
    disp(['与最优解的偏差量 = ', num2str(dr), '%']);
    disp('最优路径如图所示');
    time1 = 1:time + 1;
    
    % 画出迭代过程中的最小值散点图和总适应度折线图
    disp('其余相关图像一并展示');
    figure
    plot(time1, min_dis, 'k.');
    grid on;
    title('每代最小值散点图');
    xlabel('迭代次数');
    ylabel('最短距离');
    
    figure
    plot(time1, pop_fit_aver);
    grid on;
    title('总适应度折线图');
    xlabel('迭代次数');
    ylabel('每代总适应度');
    
    % 画出最优路径图
    figure
    DrawPath(city_coordinate, pop_min, min_index, N)
    grid on;
    title('最优路径图');
    xlabel('x');
    ylabel('y');
end


% 计算城市间距离矩阵
function [city_distance] = CityDistance(city_coordinate, N)
    city_distance = zeros(N, N);
    for i = 1:N
        for j = 1:N
            % 计算城市间的欧几里得距离
            city_distance(i, j) = ((city_coordinate(i, 1) - city_coordinate(j, 1))^2 + ...
                (city_coordinate(i, 2) - city_coordinate(j, 2))^2)^0.5;
        end
    end
end

% 计算种群的适应度
function [individual_fit, num, min_distance, a] = GroupFit(city_distance, N, pop, s)
    individual_distance = zeros(s, 1);
    for j = 1:s
        sum_distance = 0;
        for i = 1:N-1
            % 计算每个染色体的路径长度
            sum_distance = sum_distance + city_distance(pop(j, i), pop(j, i+1));
        end
        sum_distance = sum_distance + city_distance(pop(j, N), pop(j, 1));
        individual_distance(j, 1) = sum_distance;
    end
    [min_distance, a] = min(individual_distance);
    individual_fit = 1 ./ individual_distance; % 适应度是路径长度的倒数
    num = 0;
    for i = 1:s
      num = num + individual_fit(i, 1); % 总适应度
    end
end

% % 父代精英选择
% function [pop_ok] = ChooseParents(pop, N, s, c)
%     % 按照染色体适应度排序
%     pop = sortrows(pop, N + 1);
%     
%     % 选择精英个体
%     for i = 1:c
%         pop(i, :) = pop(s + 1 - i, :);
%     end
%     
%     % 随机排列种群
%     randIndex = randperm(size(pop, 1));
%     pop = pop(randIndex, :);
%     
%     % 返回选择后的种群
%     pop_ok = pop;
% end


% 轮盘赌精英选择
% 输入参数：
%   - pop: 种群
%   - N: 染色体长度
%   - s: 种群大小
%   - c: 精英数量
% 输出参数：
%   - pop_ok: 选择后的种群
function [pop_ok] = ChooseParents(pop, N, s, c)
    % 轮盘赌选择
    fit = pop(:, N + 1);
    fit_sum = sum(fit);
    fit_prob = fit / fit_sum; % 计算每个个体的选择概率
    cum_prob = cumsum(fit_prob); % 计算累积概率
    r = rand(s, 1); % 生成s个随机数
    [~, idx] = histc(r, [0; cum_prob]); % 根据随机数选择个体
    pop = pop(idx, :);
    
    % 精英选择策略
    pop = sortrows(pop, N + 1); % 再次排序
    elite = pop(s + 1 - c:s, :); % 选择最好的c个个体
    pop(1:c, :) = elite; % 将精英个体放回种群
    
    % 随机打乱种群
    randIndex = randperm(size(pop, 1));
    pop = pop(randIndex, :);
    
    pop_ok = pop;
end

% 染色体循环交叉
% 输入参数：
%   pop1: 父代个体1的染色体
%   pop2: 父代个体2的染色体
%   N:    染色体长度
% 输出参数：
%   a:    子代个体1的染色体
%   b:    子代个体2的染色体
function [a, b] = CrossVariation(pop1, pop2, N)
    % 初始化子代
    a = zeros(1, N);
    b = zeros(1, N);
    
    % 随机选择一个交叉点
    crosspoint = randi([1, N]);
    
    % 将交叉点的基因复制到子代
    a(crosspoint) = pop1(crosspoint);
    b(crosspoint) = pop2(crosspoint);
    
    % 循环交叉
    while 1
        % 在父体B中找到与子代A交叉点相同的基因的位置
        index_in_pop2 = find(pop2 == a(crosspoint));
        
        % 如果该位置在子代A中已经有基因，则停止循环
        if a(index_in_pop2) ~= 0
            break;
        end
        
        % 将父体A在该位置的基因复制到子代A
        a(index_in_pop2) = pop1(index_in_pop2);
        
        % 更新交叉点
        crosspoint = index_in_pop2;
    end
    
    % 将父体B中剩余的基因复制到子代A的剩余位置
    a(a == 0) = pop2(a == 0);
    
    % 对于子代B，重复上述步骤
    crosspoint = find(b == pop1(crosspoint));
    while 1
        index_in_pop1 = find(pop1 == b(crosspoint));
        if b(index_in_pop1) ~= 0
            break;
        end
        b(index_in_pop1) = pop2(index_in_pop1);
        crosspoint = index_in_pop1;
    end
    b(b == 0) = pop1(b == 0);
end

% 染色体逆转变异
% 输入参数：
%   - pop0: 待变异的染色体
%   - N: 染色体长度
% 输出参数：
%   - a: 变异后的染色体
function [a] = Mutation(pop0, N)
    % 生成两个随机的交叉点
    crosspoint = randi([1, N], 1, 2);
    crosspoint = sort(crosspoint); % 确保交叉点按顺序排列
    
    % 从染色体中取出交叉段
    sub = pop0(crosspoint(1):crosspoint(2));
    
    % 反转交叉段
    sub = flip(sub);
    
    % 将反转后的片段放回染色体中
    pop0(crosspoint(1):crosspoint(2)) = sub;
    
    % 返回变异后的染色体
    a = pop0;
end

% 画路径图
function DrawPath(city_coordinate, E_new_new, min_index, N)
    k = E_new_new(min_index, 1:N);
    plot(city_coordinate(:, 1), city_coordinate(:, 2), 'bo'); % 画城市坐标散点图
    hold on;
    for i = 1:N-1
        % 画出最优路径
        plot([city_coordinate(k(i), 1), city_coordinate(k(i+1), 1)], [city_coordinate(k(i), 2), city_coordinate(k(i+1), 2)], 'r', 'LineWidth', 2);
        test_t = num2str(i);
        text(city_coordinate(k(i), 1), city_coordinate(k(i), 2), test_t);
        hold on;
    end
    test_t = [num2str(N)];
    text(city_coordinate(k(N), 1), city_coordinate(k(N), 2), test_t);
end

% 生成城市数量位置数据及初始化样本数迭代次数
function [opt, times, s, n_citys, city_position] = tsp(n)
    % 根据城市数量选择文件路径和样本数
    if n == 101
        filename = 'TSPDataSet\eil101.tsp\eil101.tsp';
        opt = 629; % 读取数据集最优解距离
        s = 200; % 设置样本数
        times = 4000; % 设置迭代次数
    elseif n == 16
        filename = 'TSPDataSet\ulysses16.tsp\ulysses16.tsp';
        opt = 74;
        s = 100;
        times = 500;
    elseif n == 30
        % 加载城市位置坐标
        city_position = readmatrix('TSPDataSet\SelfDataSet\city30.xlsx');  
        n_citys = 30; % 初始化城市数量
        opt = 420;
        s = 100;
        times = 1000;
        return
    end
    
    % 打开文件并读取数据
    fid = fopen(filename,'rt');
    location = [];
    A = [1 2];
    tline = fgetl(fid);
    while ischar(tline)
        if(strcmp(tline,'NODE_COORD_SECTION'))
            while ~isempty(A)
                A = fscanf(fid,'%f',[3,1]);
                if isempty(A)
                    break;
                end
                location = [location; A(2:3)'];
            end
        end
        tline = fgetl(fid); 
        if strcmp(tline,'EOF')
            break;
        end
    end
    
    % 获取城市数量和位置
    [m, ~] = size(location);
    n_citys = m;
    city_position = location;
    fclose(fid);
    
    % 其他城市数量的代码...
end
