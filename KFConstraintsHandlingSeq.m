% [x,a]=KFConstraintsHandling(1,1,[1;1;1],[0;-3],1)
function [xTrunc,PTrunc,Success]=KFConstraintsHandlingSeq(xTrunc,PTrunc,D,d,x1,x2,X,Z1)
try
Success=1;
k=1;
try
    [Utrunc, Wtrunc, Vtrunc] = svd(PTrunc);
catch
    sprintf('Variance has NAN values')
    Success=0;
    return;
end
Ttrunc = Utrunc;
TTT = Ttrunc * Ttrunc';
if (norm(eye(size(TTT)) - TTT) > 1e-8)
    disp('Error - Ttrunc is not orthogonal.');
    Success=0;
    return;
end
if (norm(Utrunc*Wtrunc*Utrunc' - PTrunc) > 1e-8)
    disp('Error - SVD failed for pdf trunction');
    Success=0;
    return;
end
% figure(1);
% contour(x1,x2,Z1);
% hold on
% plot(x1,d(1).*ones(size(x1, 2),1));
% plot(x1,d(2).*ones(size(x1, 2),1));
% hold off
% figure(2);
% Z1 = mvnpdf(X,xTrunc',Wtrunc);
% Z1 = reshape(Z1,length(x2),length(x1));
% contour(x1,x2,Z1);
% hold on
% plot(x1,d(1).*ones(size(x1, 2),1));
% plot(x1,d(2).*ones(size(x1, 2),1));
% hold off
Amgs = sqrt(Wtrunc) * Ttrunc' * D(k,:)'; % n x 1, where n = number of states
[Wmgs, S] = MGS(Amgs);
S = S * sqrt(D(k,:) * PTrunc * D(k,:)') / Wmgs;
cTrunc = (d(k) - D(k,:) * xTrunc) / sqrt(D(k,:) * PTrunc * D(k,:)');
dTrunc = (d(k+1) - D(k+1,:) * xTrunc) / sqrt(D(k+1,:) * PTrunc * D(k+1,:)');
alpha = sqrt(2/pi) / (erf(dTrunc/sqrt(2)) - erf(cTrunc/sqrt(2)));
mu = alpha * (exp(-cTrunc^2/2) - exp(-dTrunc^2/2));
sigma2 = alpha * (exp(-cTrunc^2/2) * (cTrunc - 2 * mu) - exp(-dTrunc^2/2) * (dTrunc - 2 * mu)) + mu^2 + 1;

figure(1);
x3=min(-3,dTrunc):0.01:max(3,cTrunc);
Z2 = normpdf(x3,mu,sigma2);
x4=0:0.1:max(Z2);
plot(x3,Z2);
hold on
plot(x3,normpdf(x3,0,1));
plot(cTrunc.*ones(length(x4),1),x4);
plot(dTrunc.*ones(length(x4),1),x4);
legend('Truncated PDF', 'Standard Normal PDF','c_2','c_1')
hold off
zTrunc = zeros(size(xTrunc));
CovZ = eye(length(zTrunc));
if ~isnan(mu) %&& ~isinf(mu)
    zTrunc(1) = mu;
    CovZ(1,1) = sigma2;
else
    sprintf('Constraints failed x=%d,xbar=%d',xTrunc(1),xTrunc(2))
    Success=0;
    return;
end
xTrunc = Ttrunc * sqrt(Wtrunc) * S' * zTrunc + xTrunc;
PTrunc = Ttrunc * sqrt(Wtrunc) * S' * CovZ * S * sqrt(Wtrunc) * Ttrunc';
catch
    sprintf('Calculation Failure!')
    Success=0;
end
end