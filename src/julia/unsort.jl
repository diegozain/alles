function unsort(u)
# diego domenzain
# fall 2020 @ CSM
# 
# get index permutation such that
# sorted u(permutation) = unsorted u:
# 
# u,i_us = unsort(u);
# 
# ------------------------------------------------------------------------------
# sort
i_us = sortperm(u);
u = sort(u);
# get indicies for reversing sort
i_us = sortperm(i_us);
# now sorted u(i_us) == unsorted u
return u,i_us
end