%{
connectivity.Animal (manual) # animals xiaolong patched for connectivity studies
species            : enum('Monkey', 'Mouse') # species
animal_id          : int                     # animal id internal to database
-----
name=null          : varchar(20)             # name of a monkey, null for a mouse
date_of_birth=null : date                    # animal's date of birth
sex="unknown"      : enum('M','F','unknown') # sex of the animal
line="Unknown"     : enum('Unknown','SST-Cre','PV-Cre','VIP-Cre','WFS1-Cre','Wfs1-Ai9','Viaat-Ai9','PV-Ai9','SST-Ai9','VIP-Ai9','PV-ChR2-tdTomato','SST-ChR2-tdTomato','C57/BK6(WT)','Others') #  mouse line, monkey belongs to others
animal_notes=""    : varchar(4096)           # anything you want to say
animal_ts=CURRENT_TIMESTAMP : timestamp      # automatic
%}

classdef Animal < dj.Relvar   
    
end