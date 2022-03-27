# QL_Pacb
Information about the implementation of a Q-Learning Pacb mechanism for access control in cellular networks

The MATLAB file results_100_simulations.zip contains the first 100 results shown in the paper "Reinforcement Learning-Based ACB in LTE-A Networks for Handling Massive M2M and H2H Communications" L. Tello-Oquendo, D. Pacheco-Paramo, V. Pla, J. Martinez Bauset. IEEE International Conference on Communications ICC 2018.

The Q-Table defines the Q values according to a given state-action combination. The rows represent all the combinations of states, and the columns all the combinations of actions.

State s={E[Npt],CVnpt, Dnpt, Pacb}
Action a={1,..16}: define the value of Pacb


The other files are used to generate the results of the paper. 

The original LTE simulator with access control was implemented by Luis Tello https://github.com/lptelloq/LTE-A_RACHprocedure 

My contribution is an implementation of Q-Learning to adapt Pacb using the LTE RACH simulator. 

The main file is called LTEA_M_H_ACB_QL8.m This file can be used to train a system that learns the policy that optimizes the access probability.




