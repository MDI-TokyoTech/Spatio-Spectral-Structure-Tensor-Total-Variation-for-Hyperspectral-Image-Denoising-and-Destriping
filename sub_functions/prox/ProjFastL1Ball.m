function[X] = ProjFastL1Ball(X,alpha)

% The trick used here is that we compute the threshold like above, and if it is negative, 
% then y is inside the ball, so there is nothing to do and the threshold is set to zero.
x = X(:);
x = max(abs(x)-max(max((cumsum(sort(abs(x),1,'descend'),1)-alpha)./(1:size(x,1))'),0),0).*sign(x);
X = reshape(x,size(X));