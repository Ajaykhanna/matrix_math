% 
% SAMPLE INPUT
% This is the input used to execute Test C
%

pin(1).x = [0.09 1.17 1.26];
pin(2).x = [0.09 1.17 1.26];
pin(3).x = [0.09 1.17 1.26];
pin(4).x = [0.09 1.17 1.26];
pin(5).x = [0.09 1.17 1.26];
pin(6).x = [21.42];

pin(1).mat = {'mod' 'm43' 'mod'};
pin(2).mat = {'mod' 'm70' 'mod'};
pin(3).mat = {'mod' 'm87' 'mod'};
pin(4).mat = {'mod' 'uo2' 'mod'};
pin(5).mat = {'mod' 'gt'  'mod'};
pin(6).mat = {'mod'};

N = 714;
dx = 0.09;

pin_map = [4 4 5 4 4 5 4 4 5 4 4 5 4 4 5 4 4,...
		   1 2 5 3 3 5 3 3 5 3 3 5 3 3 5 2 1,...
		   6];