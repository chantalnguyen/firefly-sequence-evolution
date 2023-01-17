function fin_seq = vocab_generator_consang(w_p, w_s, num_flies)
% runs vocab generator along consangineus clade: P. ignitus, indictus, aquilonius, macdermotti,
% greeni, consangineus
% initialized with macdermotti and ignitus as fixed sequences, with 1
% free-to-evolve sequence
% then, fix the 3 and add a new sequence that can evolve
% then repeat until we achieve 6 species
rng('shuffle')
num_species = 1; % num_species is always 1 - we evolve 1 species at a time
%%
fixed_seq(1,:) = [1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]; % macdermotti
fixed_seq(2,:) = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]; % ignitus
seq_length = length(fixed_seq);
sequence_matrix = ones(num_flies,seq_length);
dt = datestr(now,30);

filename = ['consang_first2plus1_wp' num2str(w_p) '_ws' num2str(w_s) '_' dt];
fin_seq = vocab_gen_tree(sequence_matrix,fixed_seq,num_flies,num_species,w_p,w_s,filename);
fixed_seq = vertcat(fixed_seq,fin_seq); % add the output sequence to the list of fixed sequences
sequence_matrix = repmat(fin_seq,num_flies,1);
filename = ['consang_first2plus2_wp' num2str(w_p) '_ws' num2str(w_s) '_' dt];
fin_seq = vocab_gen_tree(sequence_matrix,fixed_seq,num_flies,num_species,w_p,w_s,filename);

fixed_seq = vertcat(fixed_seq,fin_seq);
sequence_matrix = repmat(fin_seq,num_flies,1);
filename = ['consang_first2plus3_wp' num2str(w_p) '_ws' num2str(w_s) '_' dt];
fin_seq = vocab_gen_tree(sequence_matrix,fixed_seq,num_flies,num_species,w_p,w_s,filename);

fixed_seq = vertcat(fixed_seq,fin_seq);
sequence_matrix = repmat(fin_seq,num_flies,1);
filename = ['consang_first2plus4_wp' num2str(w_p) '_ws' num2str(w_s) '_' dt];
fin_seq = vocab_gen_tree(sequence_matrix,fixed_seq,num_flies,num_species,w_p,w_s,filename,1);


end


function fin_seq = vocab_gen_tree(sequence_matrix,fixed_seq,num_flies,num_species,w_p,w_s,filename,save_flag)
if nargin < 9
    save_flag = 0;
end
% 2 known sequences and one of zeros
seq_length = length(fixed_seq);

scores = cell(num_flies,1);

max_iter = 12000;
mutate_prob = 0.1;
sequences = zeros(num_flies,seq_length,max_iter);

fitness_matrix = zeros(num_flies,max_iter,num_species);
fitness_orig = zeros(num_flies,1,num_species);

combos = combnk(1:num_flies*num_species,2);
for t = 1:max_iter
    for j = 1:length(combos) % iterate over every pair of flies
        combo = combos(j,:);
        fly1 = sequence_matrix(mod(combo(1)-1,num_flies)+1,:,ceil(combo(1)/num_flies));
        fly2 = sequence_matrix(mod(combo(2)-1,num_flies)+1,:,ceil(combo(2)/num_flies));
        
        if ceil(combo(1)/num_flies) == ceil(combo(2)/num_flies) % if flies are same species, compare scores
            if t > 1
                score1 = w_s*(mean(scores{combo(1)})) + calc_pred_risk(fly1,w_p);
                score2 = w_s*(mean(scores{combo(2)})) + calc_pred_risk(fly2,w_p);
                if score1 > score2
                    sequence_matrix(mod(combo(1)-1,num_flies)+1,:,ceil(combo(1)/num_flies)) = fly2;
                    if rand < mutate_prob
                        if sum(fly1) == 0 || sum(fly1) == length(fly1)
                            fly1 = flip_bit(fly1);
                        else
                            r = randi(2);
                            if r == 1
                                fly1 = transpose_bits(fly1);
                            else
                                fly1 = flip_bit(fly1);
                            end
                        end
                        sequence_matrix(mod(combo(1)-1,num_flies)+1,:,ceil(combo(1)/num_flies)) = fly1;
                        scores{combo(1)} = [];
                    end
                elseif score2 > score1
                    sequence_matrix(mod(combo(2)-1,num_flies)+1,:,ceil(combo(2)/num_flies)) = fly1;
                    if rand < mutate_prob
                        if sum(fly2) == 0 || sum(fly2) == length(fly2)
                            fly2 = flip_bit(fly2);
                        else
                            r = randi(2);
                            if r == 1
                                fly2 = transpose_bits(fly2);
                            else
                                fly2 = flip_bit(fly2);
                            end
                        end
                        sequence_matrix(mod(combo(2)-1,num_flies)+1,:,ceil(combo(2)/num_flies)) = fly2;
                        scores{combo(2)} = [];
                    end
                else % if scores are equal
                    if rand < mutate_prob % mutate one or the other
                        if randi(2) == 1
                            if sum(fly1) == 0 || sum(fly1) == length(fly1)
                                fly1 = flip_bit(fly1);
                            else
                                r = randi(2);
                                if r == 1
                                    fly1 = transpose_bits(fly1);
                                else
                                    fly1 = flip_bit(fly1);
                                end
                            end
                            sequence_matrix(mod(combo(1)-1,num_flies)+1,:,ceil(combo(1)/num_flies)) = fly1;
                            scores{combo(1)} = [];
                        else
                            if sum(fly2) == 0 || sum(fly2) == length(fly2)
                                fly2 = flip_bit(fly2);
                            else
                                r = randi(2);
                                if r == 1
                                    fly2 = transpose_bits(fly2);
                                else
                                    fly2 = flip_bit(fly2);
                                end
                            end
                            sequence_matrix(mod(combo(2)-1,num_flies)+1,:,ceil(combo(2)/num_flies)) = fly2;
                            scores{combo(2)} = [];
                            
                        end
                    end
                end
            end
        else % if flies are not same species, calculate similarity score
            distance = calc_similarity(fly1,fly2);
            score1 = scores{combo(1)};
            score2 = scores{combo(2)};
            score1 = [score1; distance];
            score2 = [score2; distance];
            scores{combo(1)} = score1;
            scores{combo(2)} = score2;
        end
    end
    for j = 1:size(fixed_seq,1) % calc similarity with the fixed sequence
        for k = 1:num_flies
            for m = 1:num_flies
                distance = calc_similarity(sequence_matrix(k,:),fixed_seq(j,:));
                score1 = scores{k};
                score1 = [score1; distance];
                scores{k} = score1;
            end
        end
    end
    if t == 1
        for i = 1:num_species
            for j = 1:num_flies
                fitness = calc_pred_risk(sequence_matrix(j,:,i),w_p);
                score = mean(scores{(i-1)*num_flies+j});
                fitness_orig(j,t,i) = w_s*score + fitness;
            end
        end
    end
    
    
    for i = 1:num_species
        for j = 1:num_flies
            fitness = calc_pred_risk(sequence_matrix(j,:,i),w_p);
            score = mean(scores{(i-1)*num_flies+j});
            fitness_matrix(j,t,i) = w_s*score + fitness;
        end
    end
    sequences(:,:,t) = sequence_matrix;
    
    
   
end

%%
fin_seq = mode(sequence_matrix(:,:,1));
if sum(fin_seq)> 0
    longest=max(accumarray(nonzeros((cumsum(~fin_seq)+1).*fin_seq),1));
    seq1_string = num2str(fin_seq);
    seq2_string = num2str(ones(1,longest));
    k = strfind(seq1_string,seq2_string);
    fin_seq = circshift(fin_seq,seq_length-floor(k/3));
end

if save_flag == 1
    save([filename '.mat']);
end
end



function longest = calc_similarity(seq1, seq2)
% calculate length of longest common sequence in two wrapped sequences
s1 = horzcat(seq1,seq1(1:end-1));
s2 = horzcat(seq2,seq2(1:end-1));
str_len = zeros(length(s2)-1,1);
for i = 1:length(s2) % shift sequence by one bit
    s2 = circshift(s2,1);
    sdiff = ~(abs(s1 - s2));
    if sum(sdiff) > 0 % if there are some shared bits
        str_len(i) = max(accumarray(nonzeros((cumsum(~sdiff)+1).*sdiff),1));
    else
        str_len(i) = 0;
    end
end
longest = max(str_len)/length(s1);


end

function pred_risk = calc_pred_risk(sequence,w1)
% fitness is defined as total time light is on (predation risk)
% plus coeff * times the flash switches on (energy cost)
prisk = sum(sequence)/length(sequence);
energy = numel(find(diff(sequence)==-1))+1;

w2 = 0; % although here i just hard-set energy cost to 0. can adjust later
pred_model = @(pred_risk,energy) (w1)*pred_risk + (w2)*energy;
pred_risk = pred_model(prisk,energy);

end

function sequence = transpose_bits(original_sequence)
% switch two adjacent bits in sequence
ind = randi(length(original_sequence)-1,1);
sequence = original_sequence;
sequence(ind+1) = original_sequence(ind);
sequence(ind) = original_sequence(ind+1);
end

function sequence = flip_bit(original_sequence)
% flip a bit in sequence
ind = randi(length(original_sequence),1);
sequence = original_sequence;
temp = sequence(ind);
temp(:) = ~temp; % changes 0 to 1 and 1 to 0
sequence(ind) = temp;
end


