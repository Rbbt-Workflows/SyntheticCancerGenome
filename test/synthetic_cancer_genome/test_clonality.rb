require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestClass < Test::Unit::TestCase
  def test_clonal_genotypes
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-1_chr1:15:G
    - copy-2_chr1:15:G
  SVs:
    - - DEL
      - copy-1_chr1
      - 10
      - 20
    - - INS
      - copy-2_chr1
      - 10
      - 20
      - copy-2_chr1
      - 21
-
  parent: founder
  mutations:
    - copy-1_chr1:30:C
    - copy-2_chr1:30:C
  SVs:
    - - INS
      - copy-1_chr1
      - 25
      - 35
      - copy-2_chr1
      - 50
-
  parent: 1
  mutations:
  SVs:
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)
    assert_include clonal_genotypes.last[:all_mutations], "copy-2_chr1:15:G"
    assert_include clonal_genotypes.last[:all_mutations], "copy-2_chr1:30:C"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-2_chr1:15:G"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-2_chr1:26:G"
  end
end

