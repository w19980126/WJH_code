function m = maskGen(I)
%% This function is used to generate masks of size(I)
    prompt = 'circle[c] or  rectangle[r]:\n';
    modle = input(prompt);
    [m,n] = size(I);
    if modle == 'c'
        m = circle_mask(m,n);
    elseif modle == 'r'
        m = retangle_mask(m,n);
    else
        m = [];
    end
    figure
    imagesc(abs(I).*m);
end

function mask = circle_mask(m,n)
    prompt = 'Please enter radius, row and column:\n';
    D = input(prompt);
    mask = zeros(m,n);
    for ii = 1:m
        for jj = 1:n
            r = sqrt((ii - D(2))^2 + (jj - D(3))^2);
            if r <= D(1)
                mask(ii,jj) = 1;
            end
        end
    end
end

function mask = retangle_mask(m,n)
    prompt = 'Please enter side length, row and column:\n';
    D = input(prompt);
    mask = zeros(m,n);
    temp = 0.5*D(1);
    for ii = 1:m
        for jj = 1:n
            if abs(ii - D(2))<=temp && abs(jj - D(2))<=temp
                mask(ii,jj) = 1;
            end
        end
    end
end


