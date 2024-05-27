require_relative 'haploid'
require_relative 'structural_variants'

module SyntheticCancerGenome

  def self.clonal_genotypes(evolution, chr_sizes = nil)

    ancestry = {}
    evolution.each_with_index do |info,i|
      IndiferentHash.setup(info)
      id = info["id"] || "clone-#{i}"

      if parent = info["parent"]
        if String === parent
          parent_info = evolution.select{|info| info["id"] == parent }.first
          parent_info ||= evolution[parent.to_i]
        else
          parent_info = evolution[parent]
        end

        parent_number = evolution.index parent_info

        ancestry[i] ||= ancestry[parent_number] + [parent_number]

      else
        ancestry[i] = []
      end
    end

    #{{{ Remove SV overlaps and out of bounds
    original_SVs = TSV.setup({}, :key_field => "SV ID", :fields => ["Type", "Chromosome", "Start", "End", "Target chromosome", "Target start", "Target end"], :type => :list)
    evolution.each do |info|
      IndiferentHash.setup(info)
      next unless info["SVs"]
      info["SVs"].each_with_index do |values,i|
        values = SyntheticCancerGenome.haploid_SV(values)
        id = Misc.digest(values.compact * ":" + ".#{i}")
        original_SVs[id] = values
      end
    end

    original_SVs = original_SVs.select do |k,values|
      type, chr, start, eend, other_chr, other_start = values
      padding = chr_sizes[:padding] || 0

      size = chr_sizes.values_at(chr, chr.sub(/copy-\d+_/,''), chr.sub(/copy-\d+_chr/,'')).compact.first
      next false if size.nil?
      next false if eend.to_i > size.to_i + padding
      next true unless other_chr && other_chr != ""

      chr = other_chr
      size = chr_sizes.values_at(chr, chr.sub(/copy-\d+_/,''), chr.sub(/copy-\d+_chr/,'')).compact.first
      next false if size.nil?
      next false if other_start.to_i > size.to_i + padding

      true
    end if chr_sizes

    clean_SVs = SyntheticCancerGenome.clean_SV_overlaps(original_SVs)

    evolution.each do |info|
      next unless info["SVs"]
      new = []
      info["SVs"].each_with_index do |values,i|
        values = SyntheticCancerGenome.haploid_SV(values)
        id = Misc.digest(values.compact * ":" + ".#{i}")
        if clean_SVs.include?(id)
          new << values
        else
          Log.low "Remove SV #{id}: #{values.inspect}"
          next
        end
      end 
      info["SVs"] = new
    end
    #}}} Remove SV overlaps


    evolution.each_with_index do |info,i|
      IndiferentHash.setup(info)
      mutations = info["mutations"] || []

      info["haploid_mutations"] = mutations.collect{|mutation| SyntheticCancerGenome.haploid_mutation(mutation) }

      svs = TSV.setup({}, :key_field => "SV ID", :fields => ["Type", "Chromosome", "Start", "End", "Target chromosome", "Target start", "Target end"], :type => :list)

      info["SVs"].each_with_index do |values,i|
        values = SyntheticCancerGenome.haploid_SV(values)
        id = Misc.digest(values.compact * ":" + ".#{i}")
        svs[id] = values
      end if info["SVs"]

      info["haploid_SVs"] = svs
    end

    evolution.each_with_index do |info,i|
      private_SVs = info["haploid_SVs"]
      private_mutations = info["haploid_mutations"]

      parent = ancestry[i].last

      if parent
        parent_mutations = evolution[parent]["all_mutations"]
        parent_SVs = evolution[parent]["all_SVs"]
        all_SVs = SyntheticCancerGenome.clean_SV_overlaps(parent_SVs.annotate(parent_SVs.merge(private_SVs)))
        transposed_SVs = SyntheticCancerGenome.transpose_SVs(parent_SVs, private_SVs)
      else
        parent_mutations = []
        parent_SVs = nil
        all_SVs = private_SVs
      end

      info["all_mutations"] = parent_mutations + private_mutations
      info["parent_SVs"] = parent_SVs
      info["all_SVs"] = all_SVs
    end

    evolution.each_with_index do |info,i|
      private_mutations = info["haploid_mutations"]
      private_SVs = info["haploid_SVs"]

      parent = ancestry[i].last
      if parent
        parent_SVs = evolution[parent]["all_SVs"]

        transposed_SVs = SyntheticCancerGenome.transpose_SVs(parent_SVs, private_SVs)

        transposed_private_mutations_pre = SyntheticCancerGenome.transpose_mutations(parent_SVs, private_mutations).values.collect{|v| choose(v) }
        transposed_private_mutations = SyntheticCancerGenome.transpose_mutations(transposed_SVs, transposed_private_mutations_pre).values.flatten

        transposed_parent_mutations_pre = evolution[parent]["transposed_mutations"]
        transposed_parent_mutations = SyntheticCancerGenome.transpose_mutations(transposed_SVs, transposed_parent_mutations_pre).values.flatten

        transposed_mutations = transposed_parent_mutations + transposed_private_mutations
      else
        private_mutations = info["haploid_mutations"]
        if info["haploid_SVs"]
          transposed_mutations = SyntheticCancerGenome.transpose_mutations(info["haploid_SVs"], private_mutations).values.flatten
        else
          transposed_mutations = private_mutations
        end
      end

      info["transposed_mutations"] = transposed_mutations
    end
  end
end
