module SyntheticCancerGenome
  input :somatic_mutations, :array, "Somatic mutations", nil, :required => true
  input :parents, :array, "Array of clone parents of all clones except founder", [0,0,1,2,3]
  input :fractions, :array, "Array of clone cell fractions in population", [0.1,0.2,0.3,0.4,0.5,0.6]
  input :mutation_fractions, :array, "Array of fractions of mutations assigned to each clone", [0.1,0.2,0.3,0.4,0.5,0.6]
  task :synthetic_clonal_evolution => :yaml do |somatic_mutations,parents,fractions,mut_fractions|
    num_clones = fractions.length

    parents = [nil] + parents if parents.length < fractions.length

    parents = parents.collect{|p| String === p && p =~ /^\d+$/ ?  p.to_i : p }

    fractions = fractions.collect{|f| f.to_f }
    mut_fractions = mut_fractions.collect{|f| f.to_f }

    total_fractions = Misc.sum(fractions)
    fractions = fractions.collect{|f| f / total_fractions }

    total_mut_fractions = Misc.sum(mut_fractions)
    mut_fractions = mut_fractions.collect{|f| f / total_mut_fractions }

    total_mutations = somatic_mutations.length
    current_mut = 0
    evolution = []

    parents.zip(fractions, mut_fractions).each do |parent,fraction,mut_fraction|

      final_mut = current_mut + (total_mutations * mut_fraction).ceil
      cmuts = somatic_mutations[current_mut..final_mut-1]
      current_mut = final_mut

      evolution << {name: "#{evolution.length+1}", "parent": parent, "fraction": fraction, "mutations": cmuts}
    end

    evolution
  end
end
