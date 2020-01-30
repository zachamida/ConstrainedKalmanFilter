Mu_1 = app.mu_xSlider.Value;
Mu_2 = app.mu_ySlider.Value;
Sigma_1=app.sigma_xSlider.Value;
Sigma_2=app.sigma_ySlider.Value;
Corr = app.CorrSlider.Value;
LB=app.LowerBoundyEditField.Value;
UB=app.UpperBoundyEditField.Value;
mu = [Mu_1 Mu_2];
Cov = [Sigma_1^2 Corr*Sigma_1*Sigma_2; Corr*Sigma_1*Sigma_2 Sigma_2^2];
x1 = -3*Sigma_1:Sigma_1/10:3*Sigma_1;
x2 = -3*Sigma_2:Sigma_2/10:3*Sigma_2;
D=[0 1;0 1];
d=[UB LB];
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];
Z = mvnpdf(X,mu,Cov);
Z = reshape(Z,length(x2),length(x1));

surf(app.UIAxes,x1,x2,Z,'FaceAlpha',0.9,'FaceColor','interp')

hold(app.UIAxes,'on');

[xTrunc,PTrunc,~]=KFConstraintsHandlingSeq(mu',Cov,D,d,x1,x2,X,Z);
try
    Z_trunc = mvnpdf(X,xTrunc',PTrunc);
    Z_trunc = reshape(Z_trunc,length(x2),length(x1));
    axis(app.UIAxes,[-3*Sigma_1  3*Sigma_1    -3*Sigma_2  3*Sigma_2    0  max(Z_trunc(:))])

    surf(app.UIAxes,x1,x2,Z_trunc,'FaceAlpha',0.2);
    hold(app.UIAxes,'off');
catch
    print('Truncation Failed!')
end
