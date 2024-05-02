Workflow.require_workflow "COSMIC"

module SyntheticCancerGenome

  helper :chromosome_sizes do |reference_code='hg38'|
    reference = HTS.helpers[:reference_file].call(reference_code)
    reference = GATK.prepare_FASTA reference
    sizes = {}
    reference.replace_extension('dict',true).read.split("\n").each do |line|
      if m = line.match(/SN:(?:chr)?(\d+|Y|X|M)\s+LN:(\d+)/)
        sizes[m[1]] = m[2].to_i
      end
    end

    sizes
  end

  input :build, :select, "Organism build", "hg38", :select_options => %w(hg19 b37 hg39 GRCh38 GRCh37)
  input :histology, :string, "Histology terms separated by column", "bladder:carcinoma"
  task :somatic_mutations_COSMIC => :array do |build,histology|
    histology_terms = histology.split(":")
    build = Organism.GRC_build(build)

    mutations = TSV.traverse COSMIC[build].mutations, 
      :type => :list, 
      :fields => [ "Genomic Mutation", "Primary site", "Site subtype 1", "Site subtype 2", "Site subtype 3", "Primary histology", "Histology subtype 1", "Histology subtype 2", "Histology subtype 3" ],
      :bar => self.progress_bar("Selecting mutations for #{histology}"),
      :into => :stream do |k,values|

        mutation, *hist = values

        next unless mutation =~ /^(?:chr)?[0-9MXY]+:\d+:[ACTG\-\+]+$/i
        next unless (histology_terms - hist).empty?
        next if mutation.split(":").last.include?("?")
        next if mutation.split(":").first.include? "_"

        mutation
      end

      CMD.cmd("sort -u > #{self.tmp_path}", :in => mutations)
      nil
  end

  dep :somatic_mutations_COSMIC
  input :chromosome, :string, "Chromosome to choose from", nil
  input :mutations_per_MB, :float, "Density of mutations in mutations per MB", 100
  task :somatic_mutations => :array do |chr,mutations_per_MB|
    dep = dependencies.first
    build = dep.inputs[:build]

    sizes = chromosome_sizes(build)

    if chr
      chr = chr.sub('chr', '') 
      size = sizes[chr]
      mutations = TSV.traverse step(:genotype_somatic_hg38_COSMIC_histology), :type => :array,
        :bar => self.progress_bar("Choosing chromosome mutations"),
        :into => [] do |mutation|
          next if chr && mutation.split(":").first.sub('chr','') != chr
          mutation
        end
    else
      size = Misc.sum(sizes.values)
      mutations = dep.load
    end

    if mutations_per_MB && mutations_per_MB > 0
      max = size * mutations_per_MB / 1_000_000

      mutations = mutations.shuffle[0..max-1] if mutations.length > max
    end
  end
end
