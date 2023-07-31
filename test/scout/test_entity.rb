require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

class TestEntity < Test::Unit::TestCase

  module TestE
    extend Entity
    property :property3 do
      "Property3"
    end
  end

  def test_entity
    TmpFile.with_path do |dir|
      dir.property1.write("Property1")
      dir.property2.write("Property2")

      
      sss 0
      TestE.directory = dir
      assert "Property1", TestE.property(:property1).read
      assert "Property2", TestE.property(:property2).read
      assert "Property3", TestE.property(:property3)
      assert "Property3", TestE.file(:property2).read
    end


  end
end

