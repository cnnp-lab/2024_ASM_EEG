function [p,h] = gmb_ShadowPlot_V2 (X,Y,col,per)

% Plots the mean value of X with a shadow of +- STD
% act indicates the points in Y that will be accentuated

L = size(X,3);

for i=1:L
    x = X(:,:,i);
    aux = nanmedian(x,1);
    [b,~] = size(x);
    if b>1
        SEM = nanstd(x,1)/sqrt(size(x,1));
        ts = tinv([per 1-per],size(x,1)-1);
        CI = SEM.*ts';
    else
        CI = zeros(2,size(x,2));
    end
    xconf{i} = [Y, flip(Y)];
    yconf{i} = [aux+CI(1,:), flip(aux+CI(2,:))]';
end

hold on

for i=1:L
    p(i) = fill(xconf{i},yconf{i},'b');
    p(i).FaceAlpha = 0.3;
    p(i).EdgeColor = 'none';
    p(i).FaceColor = col{i};
end
for i=1:L
    h(i) = plot(Y,nanmedian(X(:,:,i),1),'color',col{i},'LineStyle','-','Marker','none');
end

xlim(sort([Y(1),Y(end)]))
end
