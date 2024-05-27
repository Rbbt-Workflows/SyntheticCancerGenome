require_relative 'simulate'
require_relative '../clonality'
module SyntheticCancerGenome

  input :evolution, :text, "Clonal evolution"
  input :minify_sizes, :tsv, "Sizes of each chromosome's beggining to preserve", nil
  input :minify_padding, :integer, "Extra bases to add to reference", 5_000
  task :clonal_genotypes => :array do |evolution,minify_sizes,minify_padding|
    if Step === evolution
      evolution = evolution.load
    elsif Path === evolution
      evolution = evolution.yaml
    elsif String === evolution
      evolution = YAML.load(evolution)
    else
      evolution = evolution
    end

    minify_sizes = TSV.open(minify_sizes) if String === minify_sizes
    minify_sizes = minify_sizes.to_single if TSV === minify_sizes
    minify_sizes[:padding] = minify_padding if minify_sizes

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution, minify_sizes)

    clonal_genotypes.each_with_index do |info,i|
      dir = file("clone_#{i}")
      ['id', 'mutations', 'SVs', 'haploid_mutations', 'haploid_SVs', 'all_mutations', 'parent_SVs', 'all_SVs', 'transposed_mutations'].each do |field|
        value = info[field]
        case value
        when Array 
          if field == 'SVs'
            Open.write(dir[field], value.collect{|v| v * ":" } * "\n" + "\n")
          else
            Open.write(dir[field], value * "\n" + "\n")
          end
        when TSV
          Open.write(dir[field], value.to_s)
        end
      end
    end

    Dir.glob(files_dir + "/clone_*")
  end

  #input :fractions, :array, "Clonal fractions", nil, :required => true
  #dep :clonal_genotypes
  #dep :clone, :SVs => :placeholder, :mutations => :placeholder do |jobname,options,dependencies|
  #  clonal_genotypes = dependencies.flatten.first
  #  jobs = []
  #  fractions = options[:fractions].collect{|f| f.to_f }
  #  sum = Misc.sum(fractions)
  #  fractions = fractions.collect{|f| f / sum }
  #  fractions.each_with_index do |fraction,i|
  #    clone_dir = clonal_genotypes.file("clone_#{i}")
  #    next if fraction == 0
  #    clone_options = {:SVs => clone_dir.all_SVs, :mutations => clone_dir.transposed_mutations, :fraction => fraction, :restore_svs => clone_dir.all_SVs}
  #    clone_name = jobname + "-clone_#{i}"
  #    jobs << {:inputs => options.merge(clone_options), :jobname => clone_name}
  #  end
  #  jobs
  #end
  #input :invert_selection, :boolean, "Invert the selection of reads", false
  #task :clonal_tumor => :array do |fractions,invert_selection|

  #  Open.mkdir files_dir
  #  ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].sort.each_with_index do |filename,i|
  #    output = file(filename)

  #    sout = Misc.open_pipe false, false do |sin|

  #      fractions.each_with_index do |fraction,c|
  #        fraction = 1
  #        fraction = fraction.to_f
  #        
  #        next if fraction == 0.0

  #        clone_step = dependencies[c+1]

  #        input = clone_step.file('output').glob("*.fq.gz").sort[i]

  #        skip = nil
  #        rnd = Random.new 1234
  #        TSV.traverse input, :type => :array, :bar => self.progress_bar("Processing #{fraction}#{invert_selection ? " inverted" : ""} #{ [clone_step.short_path, filename] * " => " }") do |line|
  #          if line =~ /^@.*clone_/
  #            rand = rnd.rand

  #            if invert_selection
  #              skip = rand < 1.0 - fraction 
  #            else
  #              skip = rand > fraction 
  #            end

  #            next if skip
  #          else
  #            next if skip
  #          end

  #          sin.puts line
  #        end
  #      end

  #      sin.close

  #    end

  #    CMD.cmd(:bgzip, "-c > #{output}", :in => sout)
  #  end
  #  Dir.glob(files_dir + "/*")
  #end
  
  input :fractions, :array, "Clonal fractions", nil, :required => true
  dep :clonal_genotypes
  dep :clone, :SVs => :placeholder, :mutations => :placeholder do |jobname,options,dependencies|
    clonal_genotypes = dependencies.flatten.first
    jobs = []
    fractions = options[:fractions].collect{|f| f.to_f }
    sum = Misc.sum(fractions)
    fractions = fractions.collect{|f| f / sum }
    fractions.each_with_index do |fraction,i|
      clone_dir = clonal_genotypes.file("clone_#{i}")
      next if fraction.to_f == 0.0
      clone_options = {:SVs => clone_dir.all_SVs, :mutations => clone_dir.transposed_mutations, :fraction => fraction, :restore_svs => clone_dir.all_SVs}
      clone_name = jobname + "-clone_#{i}"
      jobs << {:inputs => options.merge(clone_options), :jobname => clone_name}
    end
    jobs
  end
  task :clonal_tumor => :array do |fractions,invert_selection|

    Open.mkdir files_dir
    ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].sort.each_with_index do |filename,i|
      output = file(filename)

      sout = Misc.open_pipe false, false do |sin|

        dependencies.each do |clone_step|
          next unless clone_step.task_name.to_sym == :clone
          input = clone_step.file('output').glob("*.fq.gz").sort[i]
          io = CMD.cmd(:gunzip, "'#{input}' -c", :pipe => true)
          Misc.consume_stream(io, false, sin, false)
        end

        sin.close
      end

      CMD.cmd(:bgzip, "-c > #{output}", :in => sout)
    end
    Dir.glob(files_dir + "/*")
  end
end
