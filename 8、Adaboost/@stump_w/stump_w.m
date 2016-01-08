function stump = stump_w

stump.threshold = 0;
stump.signum = 1;
stump.t_dim = 1;
stump=class(stump, 'stump_w') ;
%tr=class(tr, 'threshold_w', learner(idim, odim), learner_w) ;