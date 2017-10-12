            
nn=41;

hyp.cov = [0; 0];sn = 0.1;hyp.lik = log(sn);hyp.mean=[]; % Clean memory. !IMPORTANT!
            hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X(1:nn,:), Z(1:nn,i));
            exp(hyp.lik)
            nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:nn,:), Z(1:nn,i))
            
            [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:nn,:), Z(1:nn,i), X_star);

    %         %MCMC
    %         hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i));
    %         exp(hyp.lik)
    %         nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i))
    %         [m(:,i) s(:,i)] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);

            mup=m(:,i)+2*sqrt(s(:,i));
            mdown=m(:,i)-2*sqrt(s(:,i));
            figure()
            plot(m(:,i),'-');
            hold on
            plot(mup,'--');
            plot(mdown,'--');
            hold off