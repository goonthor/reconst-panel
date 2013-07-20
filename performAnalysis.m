function [guesses fs modelOut modelfactor] = performAnalysis(pnts, vals, newpnts, weights, modelIn, options)
disp('Starting analysis...')
tic
fs = options.fsize;
modelfactor = options.modelfactor;
modelOut = modelIn;
if options.flag==1
    if options.isWin
        [guesses fs] = svr_window(pnts,vals,newpnts,weights,options);
    else
        [guesses fs modelOut modelfactor] = svr(pnts,vals,newpnts,weights,modelIn,options);
    end
else
    guesses = runTriScat(pnts, vals, newpnts);
end
t=toc;
disp(['Analysis completed in ' num2str(t) ' seconds.']);

