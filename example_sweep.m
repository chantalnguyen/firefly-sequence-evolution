%% Example parameter sweep over similarity/predation weight ratio
num_species = 5;
w_p = 1; % predation weight -- usually set to one, while varying similarity weight
w_s = .2:.2:3; % similarity weight
mprob = 0.1; % mutation probability
max_iter = 1000; % number of iterations per chunk
num_runs = 15; % number of chunks total, for num_runs*max_iter total iterations
seq_length = 50; % sequence length
dstr = datestr(now,30); % current time
num_flies = 10; % number of flies
parfor i = 1:length(w_s)
    a = w_s(i);
    rng('shuffle');
    % First we initialize the vocab generator by running vocab_generator_init
    [sequences,final_seq,costs,seq_mat,init_seq,first_sequences] = vocab_generator_init(w_p,a,num_species,num_flies,mprob,max_iter,seq_length);
    fname = ['ns' num2str(num_species) '_wp' num2str(w_p) '_ws' num2str(a) '_' dstr num2str(randi(30)) '.mat'];
    svstf(fname,sequences,final_seq,costs,seq_mat,init_seq,first_sequences);
    for j = 2:num_runs
        % Then we feed the output sequences into vocab_generator_continue for
        % more chunks in order to save along the way
        [seqs,final_seq,cost,seq_mat] = vocab_generator_continue(seq_mat,w_p,a,num_species,num_flies,mprob,max_iter,seq_length);
        sequences = cat(3,sequences,seqs);
        costs = cat(2,costs,cost);
        svstf(fname,sequences,final_seq,costs,seq_mat);
    end
    
end


function svstf(filename,sequences,final_sequence,costs,sequence_matrix,init_seq,first_sequences)
if nargin > 6
    save(filename,'sequences','final_sequence','costs','sequence_matrix','init_seq','first_sequences');
else
    save(filename,'sequences','final_sequence','costs','sequence_matrix');
end
end
