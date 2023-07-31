require_relative 'clonality'

module SyntheticCancerGenome

  dep :miniref_ploidy
  dep_task :mini_clonal_tumor, SyntheticCancerGenome, :clonal_tumor, :reference => :miniref_ploidy
end
