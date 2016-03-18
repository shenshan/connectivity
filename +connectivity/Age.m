%{
connectivity.Age (computed) # compute the age of the animal
-> connectivity.Animal
-----
age       : int      # age of the animal
%}

classdef Age < dj.Relvar & dj.AutoPopulate
	
    properties
        popRel = connectivity.Animal & connectivity.Slice
    end
    
    methods(Access=protected)
		function makeTuples(self, key)
            
            birth_date = fetch1(connectivity.Animal & key, 'date_of_birth');
            slices = fetch(connectivity.Slice & key);
            
            patch_date = fetch1(connectivity.Slice & slices(1), 'patch_date');
            
            key.age = datenum(patch_date) - datenum(birth_date);
			self.insert(key)
		end
	end

end