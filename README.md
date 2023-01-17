# firefly-sequence-evolution

Code for simulating the co-evolution of firefly sequences, subject to a cost function representing the tradeoff between minimizing mutual similarity and individual predation risk.

vocab_generator_init.m initializes the vocabulary generator model from random initial conditions.

vocab_generator_continue.m runs the vocabulary generator model from an arbitrary starting point.

example_sweep.m runs a parameter sweep over the similarity/predation weight ratios for a given number of species.



For reverse-engineering the cost function,

vocab_generator_consang.m runs the vocabulary generator along the consanguineus clade, initialized to 2 earliest species.

tree_sweep.m runs a parameter sweep over similarity/predation weight ratio for the consanguineus clade.
