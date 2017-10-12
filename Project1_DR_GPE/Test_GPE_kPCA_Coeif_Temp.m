       


        scale=10000000;
        
        hyp.cov = [zeros(Dim_X+1,1);0];
        hyp.lik = log(sn);
        hyp.mean=[];
    
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i).*scale);
        exp(hyp.lik)
        nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i).*scale)        
        [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i).*scale, X_star);
        
        %wipe memory
%         hyp.cov = [zeros(Dim_X+1,1);-Inf];

        
        
        
       dim=i;         
        mu=m(:,dim);
        sigma=s(:,dim);
              
        [mu,index]=sort(mu);
        sigma=sigma(index);
%         figure(i)
        figure
        errorbar(mu,sigma,'rx')