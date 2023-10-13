clc
clear
addpath(genpath(pwd))
dizhi = 'kegg'; % kegg wiki reactome
num_pathway = length(dir(sprintf('./data/%s',dizhi)))-2;
num_cell = 8;
final_optp = ones(num_cell,num_pathway);
result_type = cell(num_cell,1);
exist1 = 1;
for celltype = 1:num_cell
    result = [];
    for q = 1:184
         if exist(sprintf('./data/%s/data_%s_new.mat',dizhi,num2str(q)))
            load(sprintf('./data/%s/data_%s_new.mat',dizhi,num2str(q)))
            exist1 = 1;
         elseif exist(sprintf('./data/%s/data_%s.mat',dizhi,num2str(q)))
            load(sprintf('./data/%s/data_%s.mat',dizhi,num2str(q)))
            exist1 = 1;
         else
            exist1 = 0;
         end
          %% 数据处理
         if size(adj,1)<4 && exist1 == 1
                fprintf('第 %s 个数据集因基因节点太少而删除\n', num2str(q));
         elseif size(adj,1)>1000 && exist1 == 1
                fprintf('第 %s 个数据集因基因节点太多而跳过\n', num2str(q));
         else
            dat.orig = datan{celltype,1};
            dat.fit = datan{celltype,2};
            optadj = adj;
            number = 1000;
           %% rewire
            tic;
            [adjlist,bug] = rewire(adj,number);    
            if bug == 1
               fprintf('第 %s 个数据集因网络结构单一而删除\n', num2str(q));
            else
            newadjlist = cell(number,1);
            newobj=zeros(number,1);
            optpvalue = 1;
            for i =-6:2
                c.i = i;
                for k = 1:1
                    [newgra,obj,par] = refer_NAE(dat,sparse(adjlist{k}),1,c);
                    newadjlist{k,1}=newgra;
                    newobj(k)=obj(end);
                end
              %% p_value
                pvalue = length(find(newobj(1:end)<=newobj(1)))/length(newobj);
                if optpvalue > pvalue
                    optpvalue = pvalue;
                    save(sprintf('./result/%s/resobj_IAV_%s_%d.mat',dizhi,num2str(q),celltype))   
                    final_optp(celltype,q) = optpvalue;
                end
              %  [~,minID] = min(newobj);
                minL(q) = min(newobj);     
               result(q,i+7) = pvalue;
            end
            t = toc;
            [a1,b1]=min(result(q,:));
            fprintf('第 %s 个数据集,第 %s 个细胞类型,p值为 %s，参数为 i = %s. \n', num2str(q), num2str(celltype),num2str(a1),num2str(b1));
            result_type{celltype,1} = result;
            end
        end
    end
end
%save (sprintf('./result/%s/finaldata_2.mat',dizhi),"final_optp","result_type")




