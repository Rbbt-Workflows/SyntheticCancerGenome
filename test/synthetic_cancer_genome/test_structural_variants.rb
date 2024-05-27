require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestStructuralVariants < Test::Unit::TestCase
  def reference
    <<-EOF
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
  end

  def load_ref(ref)
    chr = nil
    seq = {}
    ref.split("\n").each do |line|
      line.strip!
      if line.start_with?(">")
        chr = line[1..-1]
      else
        seq[chr] ||= ""
        seq[chr].concat line
      end
    end

    seq
  end

  def mutation_ref(ref, mutation)
    chr, pos = mutation.split(":")
    seq = load_ref(ref)
    seq[chr][pos.to_i-1]
  end

  def test_insertion_in_reference
    mutation = "copy-2_chr1:14:T"
    assert_equal "4", mutation_ref(reference, mutation)
    TmpFile.with_file(reference) do |ref|
      svs = [["INS", "copy-2_chr1", 3, 8, "copy-2_chr1", 9, nil]]
      new_ref = SyntheticCancerGenome.apply_SVs_reference(ref,svs).read
      assert_include new_ref, "1234567834567890"
      new_mutation = SyntheticCancerGenome.transpose_mutations(svs, [mutation])[mutation].first
      assert_equal "4", mutation_ref(new_ref, new_mutation)
    end
  end

  def test_inversion_in_reference
    mutation = "copy-2_chr1:14:T"
    assert_equal "4", mutation_ref(reference, mutation)
    TmpFile.with_file(reference) do |ref|
      svs = [
        ["INV", "copy-2_chr1", 3, 8, "copy-2_chr1", 9, nil],
        ["INS", "copy-2_chr1", 10, 12, "copy-2_chr1", 13, nil],
      ]
      new_ref = SyntheticCancerGenome.apply_SVs_reference(ref,svs).read
      assert_include new_ref, "1234567887654390"
      new_mutation = SyntheticCancerGenome.transpose_mutations(svs, [mutation])[mutation].first
      assert_equal "4", mutation_ref(new_ref, new_mutation)
    end
  end

  def test_inversion_in_reference_clean
    mutation = "copy-2_chr1:14:T"
    assert_equal "4", mutation_ref(reference, mutation)
    TmpFile.with_file(reference) do |ref|
      svs = [
        ["INV", "copy-2_chr1", 3, 8, "copy-2_chr1", 9, nil],
        ["INV", "copy-2_chr1", 10, 12, "copy-2_chr1", 13, nil],
        ["DEL", "copy-2_chr1", 10, 12, "copy-2_chr1", 13, nil],
      ]
      new_ref = SyntheticCancerGenome.apply_SVs_reference(ref,svs).read
      assert_include new_ref, "12345678876543921"
      new_mutation = SyntheticCancerGenome.transpose_mutations(svs, [mutation])[mutation].first
      assert_equal "4", mutation_ref(new_ref, new_mutation)
    end
  end



  def test_insertion
    mutations = ["copy-2_chr1:15:G"]
    svs = [["INS", "copy-2_chr1", 10, 20, "cis", 20, nil]]
    transposed_mutations = SyntheticCancerGenome.transpose_mutations(svs, mutations).values.first
    assert_equal 2, transposed_mutations.length
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

  def test_trans_multiple
    mutation = "copy-2_chr1:14:T"
    assert_equal "4", mutation_ref(reference, mutation)
    TmpFile.with_file(reference) do |ref|
      svs = [["INV", "copy-1_chr1", 1, 4, "copy-2_chr1", 5, nil]]
      svs += [["DEL", "copy-2_chr1", 6, 7, nil, nil, nil]]
      new_ref = SyntheticCancerGenome.apply_SVs_reference(ref,svs).read
      assert_include new_ref, "12344321589"
      new_mutation = SyntheticCancerGenome.transpose_mutations(svs, [mutation])[mutation].first
      assert_equal "4", mutation_ref(new_ref, new_mutation)
    end
  end

end

