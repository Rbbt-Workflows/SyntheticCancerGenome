require_relative 'simulate'
require_relative 'clonality'

module SyntheticCancerGenome

  dep :clonal_tumor, :target_depth => 60
  dep :normal, :depth => 30 do |jobname,options|
    {:inputs => options, :jobname => jobname + "-normal" }
  end
  input :normal_in_tumor_contamination, :float, "Proportion of normal contamination in tumor", 0.1
  input :tumor_in_normal_contamination, :float, "Proportion of tumor contamination in normal", 0.0
  task :tumor_normal => :array do |normal_in_tumor_contamination, tumor_in_normal_contamination|
    Open.mkdir files_dir

    clonal_tumor = step(:clonal_tumor)
    normal = step(:normal)
    tumor_depth = clonal_tumor.recursive_inputs[:target_depth]
    normal_depth = normal.recursive_inputs[:depth]

    if normal_in_tumor_contamination > 0
      ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)

        clones = [normal, clonal_tumor]
        fractions = [normal_in_tumor_contamination, 1.0 - normal_in_tumor_contamination]

        sout = Misc.open_pipe false, false do |sin|

          clones.zip(fractions).each_with_index do |v,ci|
            clone, fraction = v
            fraction = fraction.to_f
            clone_step = Step === clone ? clone : Step.new(clone)

            input = clone_step.files.select{|f| f.include?(".fq.gz") }.sort.collect{|f| clone_step.file(f) }[i]

            skip = nil
            rnd = Random.new 1234
            TSV.traverse input, :type => :array, :bar => self.progress_bar("Processing #{fraction} #{ [clone.short_path, filename] * " => " }") do |line|
              if line =~ /^@.*(clone|normal)/
                rand = rnd.rand
                skip = rand > fraction 
                next if skip
              else
                next if skip
              end

              sin.puts line
            end
          end

          sin.close

        end

        CMD.cmd(:bgzip, "-c > #{output}", :in => sout)
      end
    else
      ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)
        source_file = clonal_tumor.files.select{|f| f.include?(".fq.gz") }.sort.collect{|f| clonal_tumor.file(f) }[i]
        Open.link source_file, output
      end
    end

    if tumor_in_normal_contamination > 0
      ['normal_read1.fq.gz', 'normal_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)

        clones = [tumor, normal]
        fractions = [tumor_in_normal_contamination * (normal_depth.to_f / tumor_depth.to_f), (1 - tumor_in_normal_contamination)]

        sout = Misc.open_pipe false, false do |sin|

          clones.zip(fractions).each_with_index do |v,ci|
            clone, fraction = v
            fraction = fraction.to_f
            clone_step = Step === clone ? clone : Step.new(clone)

            input = clone_step.load.sort[i]

            skip = nil
            rnd = Random.new 1234
            TSV.traverse input, :type => :array, :bar => self.progress_bar("Processing #{fraction} #{ [clone.short_path, filename] * " => " }") do |line|
              if line =~ /^@.*(clone|normal)/
                rand = rnd.rand
                skip = rand > fraction 
                next if skip
              else
                next if skip
              end

              sin.puts line
            end
          end

          sin.close

        end

        CMD.cmd(:bgzip, "-c > #{output}", :in => sout)
      end
    else
      ['normal_read1.fq.gz', 'normal_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)
        source_file = normal.files.select{|f| f.include?(".fq.gz") }.sort.collect{|f| normal.file(f) }[i]
        Open.link source_file, output
      end
    end
    Dir.glob(files_dir + "/*")
  end

  dep :miniref_ploidy
  dep_task :mini_tumor_normal, SyntheticCancerGenome, :tumor_normal, :reference => :miniref_ploidy

  dep :reference_ploidy
  dep_task :full_tumor_normal, SyntheticCancerGenome, :tumor_normal, :reference => :reference_ploidy
end
