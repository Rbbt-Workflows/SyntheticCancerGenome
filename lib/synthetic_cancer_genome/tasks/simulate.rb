require_relative 'reference'
require_relative 'structural_variants'
module SyntheticCancerGenome


  input :germline, :array, "Germline mutations"
  input :somatic, :array, "Somatic mutations"
  task :somatic_germline => :array do |germline_mutations,somatic_mutations|
    germline_mutations = germline_mutations.list if Path === germline_mutations
    somatic_mutations = somatic_mutations.list if Path === somatic_mutations
    germline_mutations + somatic_mutations
  end

  dep :somatic_germline
  input :organism, :string, "Organism code"
  input :build, :select, "Organism build", "hg38", :select_options => %w(hg19 b37 hg39 GRCh38 GRCh37)
  dep_task :tumor, NEATGenReads, :NEAT_simulate_DNA, :mutations => :somatic_germline do |jobname,options|
    build = options[:build] || Organism.organism_to_build(options[:organism])

    options[:reference] ||= HTS.helpers[:reference_file].call(build) if build

    options
  end

  input :germline, :array, "Germline mutations"
  input :organism, :string, "Organism code"
  input :build, :select, "Organism build", "hg38", :select_options => %w(hg19 b37 hg39 GRCh38 GRCh37)
  dep_task :normal, NEATGenReads, :NEAT_simulate_DNA, :mutations => :germline do |jobname,options|
    build = options[:build] || Organism.organism_to_build(options[:organism])

    options[:reference] ||= HTS.helpers[:reference_file].call(build) if build

    {inputs: options}
  end

  input :reference, :path, "Haploid reference file before SVs", nil, :nofile => true
  input :mutations, :list, "Final mutations, corrected by SVs", nil, :nofile => true
  input :germline, :array, "Germline mutations before SV correction"
  input :SVs, :tsv, "SVs to apply to reference"
  dep :SV_reference, :reference => :reference, :SVs => :SVs
  dep :SV_germline, :germline => :germline, :SVs => :SVs
  dep :somatic_germline, :somatic => :mutations, :germline => :SV_germline
  input :target_depth, :integer, "Target population depth", 60
  input :fraction, :float, "Clonal fraction", nil, :required => true
  dep_task :clone, NEATGenReads, :NEAT_simulate_DNA, :mutations => :somatic_germline, :organism => nil, :build => nil, 
    :restore_SVs => :SVs, :sample_name => nil, :depth => nil, :haploid_reference => true, :reference => :SV_reference do |id,options|

    options[:depth] = (options[:target_depth].to_f * options[:fraction].to_f).ceil
    {inputs: options}
  end
end
