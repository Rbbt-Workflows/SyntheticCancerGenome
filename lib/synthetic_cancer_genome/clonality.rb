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
        transposed_SVs = all_SVs
      end

      info["all_mutations"] = parent_mutations + private_mutations
      info["parent_SVs"] = parent_SVs
      info["all_SVs"] = all_SVs
      info["transposed_SVs"] = transposed_SVs
    end

    evolution.each_with_index do |info,i|
      private_mutations = info["haploid_mutations"]
      private_SVs = info["haploid_SVs"]
      parent_SVs = info["parent_SVs"]
      transposed_SVs = info["transposed_SVs"]

      parent = ancestry[i].last
      if parent

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

  def self.clonal_segments(svs, transposed_mutations)
    clone_segments = {}

    collect_fragments, insert_fragments, remove_fragments, duplications = SyntheticCancerGenome.SV_regions(svs)

    chromosome_alterations = {}
    insert_fragments.each do |fragment,chr,start,length,inverse|
      chromosome_alterations[chr] ||= []
      chromosome_alterations[chr] << [:insert, start, fragment, inverse]
    end

    remove_fragments.each do |chr, list|
      chromosome_alterations[chr] ||= []
      list.each do |start,eend,id|
        chromosome_alterations[chr] ||= []
        chromosome_alterations[chr] << [:remove, start, eend]
      end
    end

    chromosome_alterations.each do |chr,list|
      chr_segments = []
      pos = 1
      list.sort_by!{|a| a[1] }.each do |type,start,*other|
        chr_segments << [pos, start-1] * "-" if start > pos
        case type
        when :insert
          fragment, inverse = other
          if inverse
            chr_segments << "(" + collect_fragments[fragment].values_at(0,2,1) * "," + ")"
          else
            chr_segments << "(" + collect_fragments[fragment] * "," + ")"
          end
          pos = start
        when :remove
          pos = other.first+1
        end
      end
      chr_segments << [pos, "END"] * "-"
      clone_segments[chr] = chr_segments
    end

    alt_shift = 0
    last_pre_size = 0
    TSV.traverse transposed_mutations.
      sort_by{|m| m.split(":")[1].to_i }, 
      bar: true do |mutation|
        chr, pos, alt = mutation.split(":")
        pos = pos.to_i

        if alt.start_with?("+")
          clean_alt = alt[1..-1]
          mutation_substitutes = 0
        elsif alt.start_with?("-")
          clean_alt = alt.gsub("-",'')
          mutation_substitutes = alt.length - clean_alt.length
        else
          clean_alt = alt
          mutation_substitutes = 1
        end

        current = clone_segments[chr] ||  [["1", "END"] * "-"]
        new = []
        pre_size = 0
        alt_shift = 0
        current.each do |segment|
          advance = 0

          if segment.start_with?("(")
            source_chr, start, eend = segment[1..-2].split(",")
          elsif segment.include?("-")
            start, eend = segment.split("-")
          else
            change, _sep, segment_substitues = segment.partition("_") 
            segment_substitues = segment_substitues.to_i
            advance = segment_substitues
            size = change.length
          end

          if size.nil?
            start = start.to_i
            advance = eend == "END" ? 0 : eend.to_i - start.to_i + 1
            eend = eend == "END" ? Float::INFINITY : eend.to_i

            inverse = eend < start
            if inverse
              size = start - eend + 1
            else
              size = eend - start + 1
            end
          end

          pre_final_size = pre_size
          if ! ((pos > pre_final_size) && (pos <= pre_final_size + size))
            new << segment
          else
            last_pre_size = pre_final_size

            if segment.start_with?("(")
              offset = pos - (pre_final_size + 1)
              raise "One mutations breaks segment boundary #{mutation} #{segment} #{offset}: #{current}" if offset < 0

              if ! inverse
                new << "(" + [source_chr, start, start + offset - mutation_substitutes ] * "," + ")" if offset > 0
                new << [clean_alt, mutation_substitutes] * "_"
                new << "(" + [source_chr, start + offset + 1, eend] * "," + ")" if eend >= start + offset + 1
              else
                new << "(" + [source_chr, start, start - offset + 1] * "," + ")" if offset > 0
                new << [clean_alt, mutation_substitutes] * "_"
                new << "(" + [source_chr, start - offset - mutation_substitutes, eend] * "," + ")" if eend <= start - offset - 1
              end

            elsif segment.include?("-")
              offset = pos - pre_final_size - 1

              new << [start, start + offset - mutation_substitutes] * "-" if offset >= mutation_substitutes
              new << [clean_alt, mutation_substitutes] * "_"
              if eend == Float::INFINITY
                new << [start + offset+1, "END"] * "-"
              else
                new << [start + offset+1, eend] * "-" if eend >= start + offset + 1
              end
            else
              offset = pos - pre_final_size - 1

              new_change = ""
              new_change << change[0..offset-2]
              new_change << clean_alt
              new_change << change[offset+1..-2]

              segment_substitues += mutation_substitutes - 1

              new << [new_change, segment_substitues] * "_"

            end
          end

          pre_size += advance
        end
        new = new.compact.reject{|s| s.empty? }
        clone_segments[chr] = new
      end

    clone_segments.each{|chr,list| list.each{|s| s.gsub!(/_\d+/,'') } }

    clone_segments
  end

  def self.segments(clonal_genotypes, transposed_mutations = nil)
    segments = []

    clonal_genotypes.each_with_index do |info, i|
      segments << clonal_segments(info[:all_SVs], info[:transposed_mutations])
    end

    segments
  end
end
