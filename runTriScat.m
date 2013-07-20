function guesses = runTriScat(pnts, vals, newpnts)
if size(pnts,2) == 2
        switch size(vals,2)
            case 1
                F = TriScatteredInterp(pnts(:,1), pnts(:,2), vals);
                guesses = F(newpnts(:,1),newpnts(:,2));
            case 2
                F1 = TriScatteredInterp(pnts(:,1), pnts(:,2), vals(:,1));
                F2 = TriScatteredInterp(pnts(:,1), pnts(:,2), vals(:,2));
                guesses(:,1) = F1(newpnts(:,1),newpnts(:,2));
                guesses(:,2) = F2(newpnts(:,1),newpnts(:,2));
            otherwise
                disp('Bad data format.');
        end
elseif size(pnts,2) == 3
        switch size(vals,2)
            case 1
                F = TriScatteredInterp(pnts(:,1), pnts(:,2), pnts(:,3), vals);
                guesses = F(newpnts(:,1),newpnts(:,2),newpnts(:,3));
            case 3
                F1 = TriScatteredInterp(pnts(:,1), pnts(:,2), pnts(:,3), vals(:,1));
                F2 = TriScatteredInterp(pnts(:,1), pnts(:,2), pnts(:,3), vals(:,2));
                F3 = TriScatteredInterp(pnts(:,1), pnts(:,2), pnts(:,3), vals(:,3));
                guesses(:,1) = F1(newpnts(:,1),newpnts(:,2),newpnts(:,3));
                guesses(:,2) = F2(newpnts(:,1),newpnts(:,2),newpnts(:,3));
                guesses(:,3) = F3(newpnts(:,1),newpnts(:,2),newpnts(:,3));
            otherwise
                disp('Bad data format.');
        end
else
    disp('Cannot run TriScat in 4D.')
end