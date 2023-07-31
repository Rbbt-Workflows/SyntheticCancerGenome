require 'synthetic_cancer_genome/minify'
require 'synthetic_cancer_genome/haploid'
require 'synthetic_cancer_genome/structural_variants'

module SyntheticCancerGenome

  input :reference, :binary, "Reference file", nil, :nofile => true
  input :organism, :string, "Organism code, no build", "Hsa"
  input :build, :select, "Organism build", "hg38", :select_options => %w(hg19 b37 hg39 GRCh38 GRCh37)
  input :minify_sizes, :tsv, "Sizes of each chromosome's beggining to preserve", nil, :required => true
  input :minify_padding, :integer, "Extra bases to add to reference", 5_000
  input :minify_vcf, :boolean, "Minimize also the vcfs", false
  extension 'fa.gz'
  task :miniref => :binary do |reference,organism,build,minify_sizes,minify_padding,minify_vcf|

    if minify_sizes
      minify_sizes = minify_sizes.to_single if TSV === minify_sizes

      minify_sizes[:padding] = minify_padding

      output = file(build)

      reference_path = Rbbt.share.organisms[organism][build]
      reference_path.produce
      files = reference_path.glob_all("**/*")

      files_info = files.collect{|file| [file, file.sub(reference_path.find, '')] * "<=>" }

      cpus = config :cpus, :miniref, :default => 3
      TSV.traverse files_info, :type => :array, :bar => self.progress_bar("Minifying reference files"), :cpus => cpus do |file_info|
        file ,_sep, subpath = file_info.partition("<=>")

        target = output[subpath].find.remove_extension('gz')
        type = case file
               when /\.vcf(?:\.gz)?$/
                 next unless minify_vcf
                 SyntheticCancerGenome.minify_vcf(file, target, minify_sizes)

               when /\.fa(?:sta)?(?:\.gz)?$/
                 SyntheticCancerGenome.minify_fasta(file, target, minify_sizes)
               else
                 next
               end
        CMD.cmd(:bgzip, target)
        nil
      end
    else
      Open.ln_s reference_path, output
    end

    Open.link output["#{build}.fa.gz"], self.tmp_path
    nil
  end


  dep :miniref
  input :haploid_chromosomes, :array, "Chromosomes to not duplicate", %w(M X Y)
  extension 'fa.gz'
  task :miniref_ploidy => :binary do |haploid_chromosomes|
    haploid_chromosomes = haploid_chromosomes.collect{|c| c.sub('chr', '') }
    sout = Misc.open_pipe do |sin|
      TSV.traverse step(:miniref), :type => :array do |line|
        if line =~ />/
          sin.puts ">copy-1_" + line[1..-1]
        else
          sin.puts line
        end
      end

      skip = false
      TSV.traverse step(:miniref), :type => :array do |line|
        if line =~ />/
          chr = line[1..-1].split(/\s+/).first.sub('chr', '')
          if haploid_chromosomes.include?(chr)
            skip = true
          else
            skip = false
          end
          sin.puts ">copy-2_" + line[1..-1] unless skip
        else
          sin.puts line unless skip
        end
      end
    end

    CMD.cmd(:bgzip,"-c > #{self.tmp_path}", :in => sout)
    nil
  end
end
