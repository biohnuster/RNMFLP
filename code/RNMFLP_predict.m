function [scoreMatrix_2] = RNMFLP_predict(k,beita,gama,iterate)
    interactions_ori=importdata('../data/association.xlsx');
    disSim_ori=importdata('../data/disease semantic similarity.txt');
    [rows,cols]=size(interactions_ori);
    % 计算disease,cicRNA高斯核相似性
    [cgk,dgk] = gkl(rows,cols,interactions_ori);
    % 计算circRNA的功能相似性
    circSim_ori = circRNASS( interactions_ori, 0.5*(disSim_ori+dgk));
    % 整合 disease,circRNA的相似性矩阵
    circSimi = zeros(rows,rows);
    disSimi = zeros(cols,cols);
    for i=1:rows
        for j=1:rows
            if circSim_ori(i,j)~=0
                circSimi(i,j)=circSim_ori(i,j);
            else
                circSimi(i,j)=cgk(i,j);
            end
        end
    end
    for i=1:cols
        for j=1:cols
            if disSim_ori(i,j)~=0
                disSimi(i,j)=disSim_ori(i,j);
            else
                disSimi(i,j)=dgk(i,j);
            end
        end
    end
    
    % 模型预测
    [scoreMatrix_2] = improve_NMFLP(interactions_ori,circSimi,disSimi,beita,gama,k,iterate);
    save scoreMatrix_2 scoreMatrix_2;
end
    
    
    
