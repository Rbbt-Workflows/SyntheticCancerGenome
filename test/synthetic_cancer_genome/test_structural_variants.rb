require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestStructuralVariants < Test::Unit::TestCase
  def test_insertion_in_reference
    reference =<<-EOF
>copy-1_chr1
1234567890
1234567890
1234567890
1234567890
>copy-2_chr1
1234567890
1234567890
1234567890
1234567890
    EOF
    TmpFile.with_file(reference) do |ref|
      svs = [["INS", "copy-2_chr1", 3, 8, "copy-2_chr1", 9, nil]]
      assert_include SyntheticCancerGenome.apply_SVs_reference(ref,svs).read, "1234567834567890"
    end
  end

  def test_insertion
    mutations = ["copy-2_chr1:15:G"]
    svs = [["INS", "copy-2_chr1", 10, 20, "cis", 20, nil]]
    assert_equal 2, SyntheticCancerGenome.transpose_mutations(svs, mutations).values.first.length
  end

  def test_transpose_svs_insert
    svs1 = SyntheticCancerGenome.place_SVs([["INS", "copy-1_chr1", 10, 20, "copy-1_chr1", 21, nil]])
    svs2 = SyntheticCancerGenome.place_SVs([["INS", "copy-2_chr1", 10, 20, "copy-1_chr1", 40, nil]])
    assert_equal ["INS", "copy-2_chr1", "10", "20", "copy-1_chr1", "51"], SyntheticCancerGenome.transpose_SVs(svs1, svs2).values.first.compact.collect{|v| v.to_s }
  end

  def test_transpose_svs_delete
    svs1 = SyntheticCancerGenome.place_SVs([["DEL", "copy-1_chr1", 10, 20]])
    svs2 = SyntheticCancerGenome.place_SVs([["INS", "copy-2_chr1", 10, 20, "copy-1_chr1", 40, nil]])
    assert_equal ["INS", "copy-2_chr1", "10", "20", "copy-1_chr1", "29"], SyntheticCancerGenome.transpose_SVs(svs1, svs2).values.first.compact.collect{|v| v.to_s }
  end

  def test_transpose_svs_delete_twice
    svs1 = SyntheticCancerGenome.place_SVs([["DEL", "copy-1_chr1", 10, 20], ["DEL", "copy-1_chr1", 25, 30]])
    svs2 = SyntheticCancerGenome.place_SVs([["INS", "copy-2_chr1", 10, 20, "copy-1_chr1", 40, nil]])
    assert_equal ["INS", "copy-2_chr1", "10", "20", "copy-1_chr1", "23"], SyntheticCancerGenome.transpose_SVs(svs1, svs2).values.first.compact.collect{|v| v.to_s }
  end
end

