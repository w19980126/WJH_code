function Iout = matmove(Iin,displc)
% 本函数对目标矩阵进行平移
% Iin是输入的矩阵
% displc是一个一维向量，第一个元素指定行变动，第二个元素
% 指定列变动，比如（-5,3）表示向上移动5行，向右移动3列
% 正值表示向行/列更大的方向移动，负值仿此
% Iout是变换之后的矩阵
    Iout = rowmove(Iin,displc(1));
    Iout = colmove(Iout,displc(2));
end

function Iout = rowmove(Iin,r)
    E = eye(size(Iin,1));
    trnsmat = eye(size(Iin,1));
    Iout = Iin;
    r  = round(r);
    if r>0
        trnsmat(:,1:size(Iin,1)-r) = E(:,r+1:size(Iin,1));
        trnsmat(:,size(Iin,1)-r+1:size(Iin,1)) = E(:,1:r);
    elseif r<0
        r = abs(r);
        trnsmat(1:size(Iin,1)-r,:) = E(r+1:size(Iin,1),:);
        trnsmat(size(Iin,1)-r+1:size(Iin,1),:) = E(1:r,:);
    end
    Iout = trnsmat*Iin;
end

function Iout = colmove(Iin,r)
    E = eye(size(Iin,2));
    trnsmat = eye(size(Iin,2));
    Iout = Iin;
    r  = round(r);
    if r<0
        r = abs(r);
        trnsmat(:,1:size(Iin,1)-r) = E(:,r+1:size(Iin,1));
        trnsmat(:,size(Iin,1)-r+1:size(Iin,1)) = E(:,1:r);
    elseif r>0
        trnsmat(1:size(Iin,1)-r,:) = E(r+1:size(Iin,1),:);
        trnsmat(size(Iin,1)-r+1:size(Iin,1),:) = E(1:r,:);
    end
    Iout = Iin*trnsmat;
end
