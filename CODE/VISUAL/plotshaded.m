function varout = plotshaded(x,y,fstr, alph)
% x: x coordinates
% y: either just one y vector, or 2xN or 3xN matrix of y-data
% fstr: format ('r' or 'b--' etc)
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r', 0.7);

if size(y,1)>size(y,2)
    y=y';
end

if size(y,1)==1 % just plot one line
    plot(x,y,fstr);
end

if size(y,1)==2 %plot shaded area
    px=[x,fliplr(x)]; % make closed patch
    py=[y(1,:); fliplr(y(2,:))]';
    varout =  patch(px,py,1,'FaceColor',fstr,'EdgeColor','none', 'FaceAlpha', alph);
end

if size(y,1)==3 % also draw mean
    px=[x,fliplr(x)];
    py=[y(1,:), fliplr(y(3,:))];
    varout = patch(px,py,1,'FaceColor',fstr,'EdgeColor','none', 'FaceAlpha', alph);
    plot(x,y(2,:),fstr);
end

% alpha(alph); % make patch transparent
end
