clc
clear
addpath(genpath(pwd))
qlist =[302 502 1002;
        303 503 1003;
        304 504 1004;
        305 505 1005];   % simulation datasets
result_local = 'simu';          % 
dataset = "302";         % Select a dataset to be evaluated.

% data information
num_gene = str2double(dataset{1}(1:strlength(dataset)-1));
num_cell = str2double(dataset{1}(strlength(dataset)));
fprintf("\n *********************   scNAE Result  *************************** \n");
fprintf('             Dataset #%s contains %d genes and %d cells.             \n',dataset,num_gene,num_cell)
final_optp = ones(num_cell);
result_type = cell(num_cell,1);

%% run main

for celltype = 1:num_cell
    result = [];
    load(sprintf('./simu_data/simu_data_%s.mat',dataset))
  %% data preprocessing
    dat.orig = datan{celltype,1};
    dat.fit = datan{celltype,2};
    adj(adj>1)=1;
    optadj = adj;
    number = 1000;
   %% rewire
    tic;
    [adjlist,bug] = rewire(adj,number,1);    
    if bug == 1
       fprintf('The %s dataset is deleted because the network structure is too simple.\n', num2str(q));
    end

   %% run scNAE
    if bug ~=1
    newadjlist = cell(number,1);
    newobj=zeros(number,1);
    optpvalue = 1;
    LP  = 1 ;
    for i =-6:2
        c.i = i;
        parfor k = 1:number
            [newgra,obj,par] = refer_scNAE(dat,sparse(adjlist{k}),LP,c);
            newadjlist{k,1}=newgra;
            newobj(k)=obj(end);
        end
    % calculate p-value
        pvalue = length(find(newobj(1:end)<=newobj(1)))/length(newobj);
        if optpvalue > pvalue
            optpvalue = pvalue;
            save(sprintf('./result/%s/resobj_simu_%s_%d.mat',result_local,dataset,celltype))   
            final_optp(celltype) = optpvalue;
        end
       result(i+7) = pvalue;
    end 
    t = toc;
    [a1,b1]=min(result);
    fprintf('Cell type %s, p-value is %s, log10(mu) = %s. \n', num2str(celltype),num2str(a1),num2str(b1));
    result_type{celltype} = result;
    end
    A1 = adjlist{1};   % 原始
    A2 = newadjlist{find(newobj == min(newobj))};
    A21 = adjlist{find(newobj == min(newobj))};  % 最优
    A22 = A2.*A21;
    fprintf('最优网络为 ID %d, 变动边数为 %d \n',find(newobj == min(newobj)),sum(sum(abs((A1~=0)-(A21~=0)))))
end