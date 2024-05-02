Workflow.require_workflow "Genomes1000"

module SyntheticCancerGenome

  input :build, :select, "Organism build", "hg38", :select_options => %w(hg19 b37 hg39 GRCh38 GRCh37)
  input :population, :select, "Population AFs to use when sampling mutations", "EUR", :select_options => %w(EUR)
  task :germline_mutations => :array do |build,population|
    build = Organism.GRC_build(build)
    TSV.traverse Genomes1000[build].mutations, :fields => ["#{population}_AF"], :type => :single, :into => :stream, :bar => self.progress_bar("Selecting germline mutations") do |mutation,af|
      next unless rand < af.to_f
      mutation
    end
  end
end
