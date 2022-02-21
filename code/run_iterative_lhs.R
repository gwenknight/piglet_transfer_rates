### Iterative lhs run 

### INITIAL CONDITON 1
Initial.Values = 

#run LHS on this + / - 100%
  m <- lhs_build_run(Initial.Values, 1, "scnplay", nsamples = 5)
  
Find max1 = max likelihood
Run LHS on this + / - 100% 
Find max2 = max likelihood

If max 2 more than 10% different from max1 then run +/- 100% again
Etc until reach peak
Check different save_places

If not then run +/-50%, +/- 33%, +/-10%

Combine all of above

Take top 1000 likelihoods and look at parameter sets
