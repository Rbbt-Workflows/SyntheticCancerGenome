require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestClass < Test::Unit::TestCase

  def test_clonal_segments_simple_two
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:-G
  SVs:
    - - INS
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-14"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,11,14)"
  end

  def test_clonal_segments_simple
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:A
  SVs:
    - - INS
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-14"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,11,14)"
  end

  def test_clonal_segments_insertion
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:+GGG
  SVs:
    - - INS
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    - - INS
      - copy-2_chr2
      - 11
      - 20
      - copy-2_chr2
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-15"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,11,15)"
  end

  def test_clonal_segments_deletion
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:--G
  SVs:
    - - INS
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    - - INS
      - copy-2_chr2
      - 11
      - 20
      - copy-2_chr2
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-13"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,11,13)"
  end

  def test_clonal_segments_with_inversion
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:12:+G
  SVs:
    - - DEL
      - copy-1_chr1
      - 11
      - 20
    - - INV
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-12"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,12,11)"
  end


  def test_clonal_segments_double_hit
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:A
    - copy-2_chr1:15:G
  SVs:
    - - INS
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-14"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,11,14)"
  end


  def test_clonal_segments_indels
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:G
    - copy-2_chr2:15:-G
    - copy-2_chr3:15:--G
    - copy-2_chr4:15:+G
    - copy-2_chr5:15:GG
  SVs:
    - - INS
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    - - INS
      - copy-2_chr2
      - 11
      - 20
      - copy-2_chr2
      - 21
    - - INS
      - copy-2_chr3
      - 11
      - 20
      - copy-2_chr3
      - 21
    - - INS
      - copy-2_chr4
      - 11
      - 20
      - copy-2_chr4
      - 21
    - - INS
      - copy-2_chr5
      - 11
      - 20
      - copy-2_chr5
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)

    assert_include segments.last["copy-2_chr1"], "1-14"
    assert_include segments.last["copy-2_chr1"], "(copy-2_chr1,11,14)"
    assert_include segments.last["copy-2_chr2"], "1-14"
    assert_include segments.last["copy-2_chr2"], "(copy-2_chr2,11,14)"
    assert_include segments.last["copy-2_chr3"], "1-13"
    assert_include segments.last["copy-2_chr3"], "(copy-2_chr3,11,13)"
    assert_include segments.last["copy-2_chr4"], "1-15"
    assert_include segments.last["copy-2_chr4"], "(copy-2_chr4,11,15)"
    assert_include segments.last["copy-2_chr5"], "1-14"
    assert_include segments.last["copy-2_chr5"], "(copy-2_chr5,11,14)"
  end

  def test_clonal_indel_multiple
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-2_chr1:15:+A
    - copy-2_chr1:17:--GG
    EOF

    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)

    segments = SyntheticCancerGenome.segments(clonal_genotypes)
    assert_equal ["1-15", "A", "GG", "18-END"], segments.last["copy-2_chr1"]
  end

  def test_clonal_genotypes
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-1_chr1:15:G
    - copy-2_chr1:15:G
    - copy-2_chr1:20:A
    - copy-2_chr1:21:T
  SVs:
    - - DEL
      - copy-1_chr1
      - 11
      - 20
    - - INS
      - copy-2_chr1
      - 11
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
      - 26
      - 35
      - copy-2_chr1
      - 51
    - - INS
      - copy-2_chr1
      - 10
      - 40
      - copy-2_chr1
      - 61
-
  parent: 1
  mutations:
  SVs:
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)
    assert clonal_genotypes.first[:all_mutations].select{|m| m.include?("copy-1_chr1:15:G") }.any?
    refute clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-1_chr1:15:G") }.any?
    assert clonal_genotypes.first[:all_mutations].select{|m| m.include?("copy-2_chr1:15:G") }.any?
    assert clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-2_chr1:15:G") }.any?
    assert clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-2_chr1:25:G") }.any?
    assert clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-2_chr1:30:A") }.any?
    assert clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-2_chr1:31:T") }.any?

    assert_include clonal_genotypes.last[:all_mutations], "copy-2_chr1:15:G"
    assert_include clonal_genotypes.last[:all_mutations], "copy-2_chr1:30:C"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-2_chr1:15:G"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-2_chr1:25:G"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-2_chr1:40:C"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-2_chr1:65:C"
    assert_include clonal_genotypes.last[:transposed_mutations], "copy-1_chr1:20:C"
  end

  def test_clonal_genotypes_with_inversion
    evolution_txt = <<-EOF
---
- 
  id: founder
  mutations:
    - copy-1_chr1:15:G
    - copy-2_chr1:12:G
  SVs:
    - - DEL
      - copy-1_chr1
      - 11
      - 20
    - - INV
      - copy-2_chr1
      - 11
      - 20
      - copy-2_chr1
      - 21
    EOF
    evolution = YAML.load(evolution_txt)

    clonal_genotypes = SyntheticCancerGenome.clonal_genotypes(evolution)
    assert clonal_genotypes.first[:all_mutations].select{|m| m.include?("copy-1_chr1:15:G") }.any?
    refute clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-1_chr1:15:G") }.any?
    assert clonal_genotypes.first[:all_mutations].select{|m| m.include?("copy-2_chr1:12:G") }.any?
    assert clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-2_chr1:12:G") }.any?
    assert clonal_genotypes.first[:transposed_mutations].select{|m| m.include?("copy-2_chr1:29:G") }.any?
  end
end

