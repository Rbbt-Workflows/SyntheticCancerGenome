require 'synthetic_cancer_genome/structural_variants'

module SyntheticCancerGenome
  input :reference, :path, "Reference file", nil, :nofile => true
  input :SVs, :path, "SVs to apply to reference", nil, :noload => true
  extension 'fa.gz'
  task :SV_reference => :binary do |reference,svs|
    reference = reference.find if Path === reference
    CMD.cmd(:bgzip, "> #{self.tmp_path}", :in => SyntheticCancerGenome.apply_SVs_reference(reference, SyntheticCancerGenome.place_SVs(svs)))
    nil
  end

  input :germline, :array, "Germline mutations to transpose", []
  input :SVs, :path, "SVs to apply to reference", nil, :noload => true
  task :SV_germline => :array do |germline, svs|
    SyntheticCancerGenome.transpose_mutations(svs, germline).values.flatten
  end
end

