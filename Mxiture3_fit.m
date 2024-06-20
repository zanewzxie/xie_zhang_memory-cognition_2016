function Model = Mxiture3_fit(ROCdata,ifplot)

lb = [0 0 0  -Inf -Inf -Inf -Inf -Inf ]; % Lower bound for parameters
ub = [1 4 1 Inf Inf Inf Inf Inf];
startvalue = [0.3 1 0.1 0.1 0.5 1 1.5 2];
options = optimoptions('patternsearch','MaxIterations',1000000,'MaxFunctionEvaluations',1000000,'MeshTolerance',1.0000e-08);

[Model.para,FVAL,EXITFLAG] = patternsearch({@Mxiture3_modelpred_ROC,ROCdata},startvalue,[],[],[],[],lb,ub,options);
[Model.NegLL output] = Mxiture3_modelpred_ROC(Model.para, ROCdata);
[Model.ROC] = Mxiture3_modelpred_ROC_curve(Model.para, ROCdata);
Model.AIC = 2*Model.NegLL+2*(size(Model.para,2)-5);
Model.BIC = 2*Model.NegLL+(size(Model.para,2)-5)*log(size(ROCdata.nhi,2)); %2*negll+k*log(n); n is the data points

npercon = (ROCdata.nhi(1)+ROCdata.nfa(1)+ROCdata.nmi(1)+ROCdata.ncr(1))/2;
if ifplot
%     figure;
    plot(output.pfa,output.phi,'o');
    hold on;plot(ROCdata.nfa/npercon,ROCdata.nhi/npercon,'*');
    hold on;plot(smooth(Model.ROC.pfa),smooth(Model.ROC.phi),'-');
    axis([0 1 0 1]);
end

end


function [NLL output] = Mxiture3_modelpred_ROC(pars, data)


Po = pars(1);    % Po
d = pars(2);    % d
Pn = pars(3);    % Pn
crit(1) = pars(4);
crit(2) = pars(5);
crit(3) = pars(6);
crit(4) = pars(7);
crit(5) = pars(8);

Nvec = data.N;
nhi = data.nhi;
nmi = data.nmi;
nfa = data.nfa;
ncr = data.ncr;

Phi = NaN(1,length(crit));
Pfa = NaN(1,length(crit));


for ict=1:length(crit)
    phi(1,ict) = (1-Po)*normcdf(crit(ict),0,1) + (Po * normcdf(crit(ict),-d,1));
    pfa(1,ict) =  (1-Pn)*normcdf(crit(ict),0,1);
end


output.phi= phi;
output.pfa= pfa;

LL =  sum(nhi .* log(phi)) + sum(nmi .* log(1-phi))  ...
    + sum(nfa .* log(pfa)) + sum(ncr .* log(1-pfa)) ;
NLL = -LL; % negative log likelihood of the parameters
end



function [output] = Mxiture3_modelpred_ROC_curve(pars, data)


Po = pars(1);    % Po
d = pars(2);    % d
Pn = pars(3);    % Pn
crit = linspace(-100, 100, 1000);

Nvec = data.N;
nhi = data.nhi;
nmi = data.nmi;
nfa = data.nfa;
ncr = data.ncr;

Phi = NaN(1,length(crit));
Pfa = NaN(1,length(crit));


for ict=1:length(crit)
    phi(1,ict) = (1-Po)*normcdf(crit(ict),0,1) + (Po * normcdf(crit(ict),-d,1));
    pfa(1,ict) =  (1-Pn)*normcdf(crit(ict),0,1);
end

output.phi= phi;
output.pfa= pfa;
end


