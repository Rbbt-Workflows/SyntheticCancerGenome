module SyntheticCancerGenome

  module Clonality
    def self.normalize_fractions(fractions)
      fractions = fractions.collect{|f| f.to_f }
      sum = Misc.sum(fractions)
      fractions.collect{|f| f / sum }
    end

    def self.pick_fraction(clones, fractions, elements)
      fractions = self.normalize_fractions(fractions)
      clones.each_with_index do |info,i|
        list = somatic_mutations[current..final-1]
        current = final

        info{"parent": parent, "fraction": fraction, "mutations": list}
      end
    end
  end
end
