% Run a parameter sweep over similarity/weight ratio for consanguineus
% clade implementation of vocabulary generator
w_p = 1;
w_s = 0.4:0.05:1;
params = combvec(w2,w_s);
parfor i = 1:length(params)
    vocab_generator_consang(w_p,w_s(i));
end
