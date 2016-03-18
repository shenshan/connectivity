%{
connectivity.Slice (manual) # slice xiaolong patched
->connectivity.Animal
slice_id                 : int              # id of slice
-----
slice_name               : varchar(30)      # name of the slice in xiaolong's system
patch_date               : date             # date of experiment
region="V1"              : enum('V1','V2','PFC','Temporal cortex','barrel cortex','subcortical area','other')  # brain region
ref1_x                   : double           # reference point to calculate slice angel
ref1_y                   : double           # y of reference point 1
ref2_x                   : double           # x of reference point 2
ref2_y                   : double           # y of reference point 2
notes=""                 : varchar(4096)    # anything you want to say

%}

classdef Slice < dj.Relvar
end