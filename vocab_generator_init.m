function [sequences,final_sequence,costs,sequence_matrix,init_seq,first_sequences] = vocab_generator_init(w_p,w_s,num_species,num_flies,mutate_prob,max_iter,seq_length)
%% Initializes the naming game from random initial conditions
% Inputs:
% w_p: predation weight
% w_s: similarity weight
% num_species: # of species
% num_flies: # of fly agents per species
% mutate_prob: probability of mutation
% max_iter: max # of iterations
% seq_length: length of sequence

% initial conditions: all sequences are random
sequence_matrix = randi(2,num_flies,seq_length,num_species)-1;
init_seq = sequence_matrix;
costs = zeros(num_flies,max_iter,num_species); % store the costs for each sequence

% keep track of similarity scores
scores = NaN*zeros(num_flies*num_species,1000);

% save to memory every 100 steps
timesaves = 100:100:max_iter;
sequences = zeros(num_species,seq_length,length(timesaves));

first_sequences = zeros(num_species,seq_length,100);
for i = 1:num_species
    first_sequences(i,:,1) = mode(sequence_matrix(:,:,i));
end
%% Naming game loop
% tic;
for t = 1:max_iter
    
    % compute similarity scores between flies of different species
    for i = 1:num_species
        for k = 1:num_flies
            for j = i+1:num_species
                for l = 1:num_flies
                    fly1 = sequence_matrix(k,:,i);
                    fly2 = sequence_matrix(l,:,j);
                    distance = calc_similarity(fly1,fly2);
                    t1 = scores((i-1)*num_flies+k,:);
                    t2 = scores((j-1)*num_flies+l,:);
                    t1(find(isnan(t1),1,'first')) = distance;
                    t2(find(isnan(t2),1,'first')) = distance;
                    scores((i-1)*num_flies+k,:) = t1;
                    scores((j-1)*num_flies+l,:) = t2;
                end
            end
        end
    end
    
    
    % compare scores between pairs of flies of the same species
    for i = 1:num_species
        for j = 1:num_flies
            for k = j+1:num_flies
                fly1 = sequence_matrix(j,:,i); 
                fly2 = sequence_matrix(k,:,i);
                
                if t > 1
                    % compute score for each fly
                    score1 = (w_s/(w_p+w_s))*(nanmean(scores((i-1)*num_flies+j,:))) + (1/(w_p+w_s))*calc_pred_risk(fly1,w_p);
                    score2 = (w_s/(w_p+w_s))*(nanmean(scores((i-1)*num_flies+k,:))) + (1/(w_p+w_s))*calc_pred_risk(fly2,w_p);
                    if score1 > score2 % if fly1 has higher score, fly1 takes on fly2's sequence 
                        sequence_matrix(j,:,i) = fly2; 
                        fly1 = fly2;
                        if rand < mutate_prob % mutate: either flip or transpose bits
                            if sum(fly1) == 0 || sum(fly1) == length(fly1) % if all bits are the same, flip a bit
                                fly1 = flip_bit(fly1); 
                            else % otherwise, decide whether to flip or transpose
                                r = randi(2); 
                                if r == 1
                                    fly1 = transpose_bits(fly1);
                                else
                                    fly1 = flip_bit(fly1);
                                end
                            end
                            sequence_matrix(j,:,i) = fly1;
                        end
                        scores((i-1)*num_flies+j,:) = NaN*scores((i-1)*num_flies+j,:);
                    elseif score2 > score1  % if fly2 has higher score, fly2 takes on fly1's sequence 
                        sequence_matrix(k,:,i) = fly1;
                        fly2 = fly1;
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
                            sequence_matrix(k,:,i) = fly2;
                        end
                        scores((i-1)*num_flies+k,:) = NaN*scores((i-1)*num_flies+k,:);
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
                                sequence_matrix(j,:,i) = fly1;
                                scores((i-1)*num_flies+j,:) = NaN*scores((i-1)*num_flies+j,:);
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
                                sequence_matrix(k,:,i) = fly2;
                                scores((i-1)*num_flies+k,:) = NaN*scores((i-1)*num_flies+k,:);
                            end
                        end
                    end
                end
                
                
            end
        end
    end
    
    % save the mode of the sequences every 100 steps
    if t < 100
        for i = 1:num_species
            first_sequences(i,:,t+1) = mode(sequence_matrix(:,:,i));
        end
    elseif mod(t,100)==0
        for i = 1:num_species
            curr_seq = mode(sequence_matrix(:,:,i));
            sequences(i,:,t/100) = curr_seq;
        end
    end
    % compute and store the cost for each sequence
    for i = 1:num_species
        for j = 1:num_flies
            score = nanmean(scores((i-1)*num_flies+j,:));
            costs(j,t,i) = (w_s/(w_p+w_s))*score + (1/(w_p+w_s))*calc_pred_risk(sequence_matrix(j,:,i),w_p);
        end
    end
    
end

%%
final_sequence = zeros(num_species,seq_length);

for i = 1:num_species
    sequence = mode(sequence_matrix(:,:,i));
    final_sequence(i,:) = sequence;
end

end

function longest = calc_similarity(seq1, seq2)
% calculate length of longest common sequence in two wrapped sequences
s1 = horzcat(seq1,seq1(1:end-1));
s2 = horzcat(seq2,seq2(1:end-1));
str_len = zeros(length(s2),1);
for i = 1:length(s2) % shift sequence by one bit
    s2 = circshift(s2,1);
    sdiff = ~(abs(s1 - s2));
    if sum(sdiff) > 0 % if there are some shared bits
        A  = nonzeros((cumsum(~sdiff)+1).*sdiff);
        edges = .5:(max(A)+.5);
        str_len(i) = max(histcounts(A, edges));
        
        %         str_len(i) = max(accumarray(nonzeros((cumsum(~sdiff)+1).*sdiff),1));
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
