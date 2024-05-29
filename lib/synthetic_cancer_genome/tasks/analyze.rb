Workflow.require_workflow "HTS"

module SyntheticCancerGenome

  input :directory, :string, "Directory with the files"
  input :sample, :select, "Sample to align", :tumor, :select_options => %w(tumor normal)
  dep_task :align, HTS, :BAM, 
    :fastq1 => :placeholder, :fastq2 => :placeholder,
    :skip_rescore => true, :skip_duplicates => true, :skip_mark_adapters => true do |jobname,options,dependencies|

      directory, sample = options.values_at :directory, :sample
      sample = sample.to_s

      fastq1 = File.join(directory, sample + '_read1.fq.gz')
      fastq2 = File.join(directory, sample + '_read2.fq.gz')

      jobname += '_' + sample.to_s

      {:inputs => options.merge(:fastq1 => fastq1, :fastq2 => fastq2), :jobname => jobname}
    end

  dep :align, sample: :tumor
  dep :align, sample: :normal
  dep_task :mutect2, HTS, :mutect2, tumor: :placeholder, normal: :placeholder do |jobname,options,dependencies|
    tumor_dep, normal_dep = dependencies.flatten.compact
    {inputs: options.merge(tumor: tumor_dep, normal: normal_dep), jobname: jobname}
  end
end
