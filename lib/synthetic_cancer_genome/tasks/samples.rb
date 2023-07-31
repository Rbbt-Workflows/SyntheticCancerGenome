module SyntheticCancerGenome

  dep :germline_mutations
  dep :miniref_ploidy
  dep_task :minified_normal, SyntheticCancerGenome, :normal, :germline_mutations => :germline_mutations, :reference => :miniref_ploidy

  dep :germline_mutations
  dep :somatic_mutations
  dep :miniref_ploidy
  dep_task :minified_tumor, SyntheticCancerGenome, :tumor, :germline_mutations => :germline_mutations, :somatic_mutations => :somatic_mutations, :reference => :miniref_ploidy
end
