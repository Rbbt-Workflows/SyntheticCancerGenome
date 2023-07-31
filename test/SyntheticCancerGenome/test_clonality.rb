require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

require 'scout'

class TestClass < Test::Unit::TestCase

  def datafile(file, type = nil, &block)
    path = Rbbt.data.example[file]

    begin
      path.send(type) 
    rescue Errno::ENOENT
      path.set_extension(type).send(type)
    end



  end

  def test_clonal_genotypes
    evolution = datafile("evolution", :yaml)
    germline = datafile("germline", :list)
    TmpFile.with_path do |genotypes|
      iii SyntheticCancerGenome.clonal_genotypes evolution, germline, genotypes
      iii genotypes.clone_0.all_mutations.list
      iii genotypes.clone_0.all_mutations.list
    end
  end
end

