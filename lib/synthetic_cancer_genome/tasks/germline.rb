Workflow.require_workflow "Genomes1000"

module SyntheticCancerGenome
  HG19_ORGANISM = "Hsa/feb2014"
  HG38_ORGANISM = "Hsa/may2009"

  input :haploid, :boolean, "Use haploid frequencies (half)", false
  task :genotype_germline_hg19_all_chr => :array do |haploid|
    Workflow.require_workflow "Genomes1000"
    TSV.traverse Genomes1000.rsids, :fields => ["Genomic Mutation", "EUR_AF"], :type => :list, :into => :stream, :bar => self.progress_bar("Selecting germline mutations") do |rsid,values|
      mutation, af = values
      af = haploid ? af.to_f / 2 : af.to_f
      next unless rand < af
      chr, pos, alt_l = mutation.split(":")
      alts = alt_l.split(",")
      alt = alts.sample(1).first
      [chr, pos, alt] * ":"
    end
  end

  input :chromosome, :string, "Chromosome to choose from", nil
  dep :genotype_germline_hg19_all_chr
  task :genotype_germline_hg19 => :array do |chr|
    chr = chr.sub('chr', '') if chr
    TSV.traverse step(:genotype_germline_hg19_all_chr), :type => :array, :into => :stream do |mutation|
      next if chr && mutation.split(":").first.sub('chr','') != chr
      mutation
    end
  end

  dep :genotype_germline_hg19
  dep_task :genotype_germline_hg38_lf, Sequence, :lift_over, :positions => :genotype_germline_hg19, :source => HG19_ORGANISM, :target => HG38_ORGANISM

  dep :genotype_germline_hg38_lf
  task :genotype_germline_hg38_lf_chr => :array do 
    chr = recursive_inputs[:chromosome]
    TSV.traverse step(:genotype_germline_hg38_lf), :type => :array, :into => :stream do |line|
      next if ! line.include?(":") || line.split(":").first.include?("_")
      next if chr && line.split(":").first != chr
      line
    end
  end

  dep :genotype_germline_hg38_lf_chr
  dep Sequence, :reference, :positions => :genotype_germline_hg38_lf_chr, :organism => HG38_ORGANISM, :vcf => false, :full_reference_sequence => false
  task :genotype_germline_hg38 => :array do 
    TSV.traverse step(:reference), :into => :stream do |mutation, reference|
      # Make sure we don't take positions that are now non-mutations, as this
      # breaks other tools downstream
      next if mutation.split(":")[2].split(",").include? reference
      mutation
    end
  end

  dep :genotype_germline_hg38, jobname: "father", haploid: true, compute: :produce
  dep :genotype_germline_hg38, jobname: "mother", haploid: true, compute: :produce
  task :diploid_genotype_germline_hg38 => :array do
    Open.open_pipe do |sin|
      TSV.traverse dependencies.first, type: :array, bar: true do |mutation|
        chr, pos, alt = mutation.split(":")
        chr = "chr" + chr unless chr.start_with?("chr")
        chr = "copy-1_" + chr
        sin.puts [chr, pos, alt] * ":"
      end
      TSV.traverse dependencies.last, type: :array, bar: true do |mutation|
        chr, pos, alt = mutation.split(":")
        chr = "chr" + chr unless chr.start_with?("chr")
        chr = "copy-2_" + chr
        sin.puts [chr, pos, alt] * ":"
      end
    end
  end
end
